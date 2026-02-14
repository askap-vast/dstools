from abc import ABC, abstractmethod
from dataclasses import dataclass

import matplotlib.colors as mcolors
import matplotlib.pyplot as plt
import numpy as np
from matplotlib import cm
from mpl_toolkits.axisartist import GridHelperCurveLinear, SubplotHost
from mpl_toolkits.axisartist.grid_finder import (
    DictFormatter,
    ExtremeFinderSimple,
    FixedLocator,
)
from sklearn.linear_model import LinearRegression

from dstools.dynamic_spectrum import LightCurve, TimeFreqSeries


def get_mask(data, sigma=2):
    mask = data.flux["L"] < sigma * data.flux_err["L"]

    isolated = mask[:-2] & mask[2:]
    mask[1:-1][isolated] = True

    return mask


@dataclass
class PoincareProjection(ABC):
    data: LightCurve
    dphi_deg: int = 30
    dtheta_deg: int = 30

    def __post_init__(self):
        _, P, _, _, _, _ = self._extract_sphere_points(self.data)
        self._set_geometry(P)

    def _extract_sphere_points(self, data: LightCurve):
        # Extract time / polarisation data from lightcurve
        period = data.ds.period if data.ds.fold else 1
        time = data.x.copy() * period
        polangle = np.deg2rad(data.polangle.copy())
        ellipticity = np.deg2rad(data.ellipticity.copy())
        polangle_err = np.deg2rad(data.polangle_err.copy())
        ellipticity_err = np.deg2rad(data.ellipticity_err.copy())

        # Mask low SNR points
        mask = get_mask(data, sigma=3)
        time[mask] = np.nan
        polangle[mask] = np.nan
        ellipticity[mask] = np.nan

        # Poincare sphere coordinates normalised by P
        Sq = np.cos(2 * ellipticity) * np.cos(2 * polangle)
        Su = np.cos(2 * ellipticity) * np.sin(2 * polangle)
        Sv = np.sin(2 * ellipticity)

        P = np.stack([Sq, Su, Sv], axis=1).astype(float)

        return time, P, polangle, ellipticity, polangle_err, ellipticity_err

    def prepare_grid(
        self,
        fig: plt.Figure = None,
        ax: plt.Axes = None,
        xmin: float = -2.2,
        xmax: float = 2.2,
        ymin: float = -2.2,
        ymax: float = 2.2,
        internal_labels: bool = False,
    ):
        phi_ticks = np.arange(0, 360 + self.dphi_deg, self.dphi_deg)
        theta_ticks = np.arange(0, 180 + self.dtheta_deg, self.dtheta_deg)
        phi_ticks = np.deg2rad(phi_ticks)
        theta_ticks = np.deg2rad(theta_ticks)

        grid_helper = GridHelperCurveLinear(
            (self._fwd, self._back),
            extreme_finder=ExtremeFinderSimple(60, 60),
            grid_locator1=FixedLocator(phi_ticks),
            grid_locator2=FixedLocator(theta_ticks),
            tick_formatter1=DictFormatter(
                {p: f"{np.rad2deg(p):.0f}°" for p in phi_ticks}
            ),
            tick_formatter2=DictFormatter(
                {t: f"{90 - np.rad2deg(t):.0f}°" for t in theta_ticks}
            ),
        )

        # Replace supplied axis with curvilinear grid embedded copy
        # or set up default fig / ax if not supplied
        fig = plt.figure(figsize=(6, 6 / 1.6)) if fig is None else fig

        if ax is None:
            ax = SubplotHost(fig, 1, 1, 1, grid_helper=grid_helper)
        else:
            ss = ax.get_subplotspec()
            ax.remove()
            ax = SubplotHost(fig, ss, grid_helper=grid_helper)

        fig.add_subplot(ax)

        # Axes configuration
        ax.set_aspect("equal", adjustable="box")

        # ax.set_xlim(xmin, xmax)
        # ax.set_ylim(ymin, ymax)

        ax.grid(
            True,
            color="0.6",
            alpha=0.6,
            linewidth=0.8,
        )

        if internal_labels:
            # Enable ticklabels on all axes to catch distinct meridians
            ax.axis["bottom"].major_ticklabels.set_visible(True)
            ax.axis["top"].major_ticklabels.set_visible(True)
            ax.axis["left"].major_ticklabels.set_visible(True)
            ax.axis["right"].major_ticklabels.set_visible(True)

            # Add labels for internal parallels that do not cross axes
            xlim, ylim = ax.get_xlim(), ax.get_ylim()

            for th in theta_ticks:
                # Get y coordinate of parallel at maximum y-extent,
                # which occurs at phi = pi / 2
                x0, y0 = self._fwd(np.pi / 2, th)

                if (xlim[0] < x0 < xlim[1]) and (ylim[0] < y0 < ylim[1]):
                    ax.text(
                        0,
                        y0,
                        f" {90 - np.rad2deg(th):.0f}°",
                        va="center",
                        ha="left",
                        fontsize=8,
                    )

        return fig, ax

    def plot(
        self,
        data: LightCurve,
        fig: plt.Figure,
        ax: plt.Axes,
        connect_points: bool = False,
        plot_errors: bool = False,
        fit_linear_model: bool = False,
        cmap: str = "plasma",
        cmin: float = 0,
        cmax: float = 1,
        **scatter_kwargs,
    ):
        time, X, Y, Xerr, Yerr = self.project_data(data)

        if connect_points:
            ax.plot(
                X,
                Y,
                color="k",
                alpha=0.4,
                lw=1,
                zorder=1,
            )

        if plot_errors:
            ax.errorbar(
                X,
                Y,
                xerr=Xerr,
                yerr=Yerr,
                ls="none",
                marker="none",
                alpha=0.15,
                lw=1,
                zorder=1,
            )

        # Combine user plot settings with defaults
        default_kwargs = {"s": 5}
        scatter_kwargs = default_kwargs | scatter_kwargs

        # Colormap normalisation
        cmap = cm.get_cmap(cmap)
        cmap = mcolors.ListedColormap(cmap(np.linspace(cmin, cmax, len(time))))

        sc = ax.scatter(
            X,
            Y,
            c=time,
            ls="-",
            cmap=cmap,
            zorder=2,
            **scatter_kwargs,
        )

        cbar = fig.colorbar(sc, ax=ax, label=f"Time [{self.data.ds.tunit}]")

        if fit_linear_model:
            # Linear regression model of GC
            model = LinearRegression()

            x = X[~np.isnan(X)].reshape(-1, 1)
            y = Y[~np.isnan(Y)]
            model.fit(x, y)

            print(model.intercept_)
            print(model.coef_)

            xmin, xmax = ax.get_xlim()
            x = np.linspace(xmin, xmax, 100)
            y = model.predict(x.reshape(-1, 1))

            ax.plot(
                x,
                y,
                ls="-",
                color="r",
                alpha=0.4,
                zorder=1,
            )

        return sc

    @abstractmethod
    def _set_geometry(self, P: np.ndarray):
        """Set up geometric projection basis."""

    @abstractmethod
    def project_data(self, data: LightCurve):
        """Project from data coordinates into projected coordinate system."""

    @abstractmethod
    def _fwd(self, longitude, colatitude):
        """Convert longitude (PA) and colatitude (ellipticity) into x, y coordinates."""

    @abstractmethod
    def _back(self, x, y):
        """Convert x, y curvilinear coordinates into longitude and colatitude."""


@dataclass
class StereographicProjection(PoincareProjection):
    def _get_stokes(self, polangle, ellipticity):
        Sq = np.cos(2 * ellipticity) * np.cos(2 * polangle)
        Su = np.cos(2 * ellipticity) * np.sin(2 * polangle)
        Sv = np.sin(2 * ellipticity)

        return Sq, Su, Sv

    def _set_geometry(self, P: np.ndarray):
        self.contact_point = np.array([0, 0, -1])
        self.e1 = np.array([1, 0, 0])
        self.e2 = np.array([0, 1, 0])

        return

    def project_data(self, data: LightCurve):
        time, _, pa, ell, pa_err, ell_err = self._extract_sphere_points(data)

        longitude = 2 * pa
        colatitude = np.pi / 2 - 2 * ell

        X, Y = self._fwd(longitude, colatitude)
        Xerr, Yerr = self.project_errors(pa, ell, pa_err, ell_err)

        return time, X, Y, Xerr, Yerr

    def project_errors(self, pa, ell, pa_err, ell_err):
        """1-sigma errors for stereographic X,Y from errors in (pa, ell), radians."""

        c2c, s2c = np.cos(2 * ell), np.sin(2 * ell)
        c2p, s2p = np.cos(2 * pa), np.sin(2 * pa)
        q = c2c * c2p
        u = c2c * s2p
        v = s2c
        p = np.stack([q, u, v], axis=1)
        a1 = p @ self.e1
        a2 = p @ self.e2
        d = p @ self.contact_point
        # Partials of p wrt ψ (pa) and χ (ell)
        dq_dpsi = -2 * c2c * s2p
        du_dpsi = 2 * c2c * c2p
        dv_dpsi = 0.0 * pa
        dq_dchi = -2 * s2c * c2p
        du_dchi = -2 * s2c * s2p
        dv_dchi = 2 * c2c
        dp_dpsi = np.stack([dq_dpsi, du_dpsi, dv_dpsi], axis=1)
        dp_dchi = np.stack([dq_dchi, du_dchi, dv_dchi], axis=1)
        # X = 2 a1 / (1+d), Y = 2 a2 / (1+d)
        da1_dpsi = dp_dpsi @ self.e1
        dd_dpsi = dp_dpsi @ self.contact_point
        da2_dpsi = dp_dpsi @ self.e2
        da1_dchi = dp_dchi @ self.e1
        dd_dchi = dp_dchi @ self.contact_point
        da2_dchi = dp_dchi @ self.e2
        inv = 1.0 / (1.0 + d)
        inv2 = inv * inv
        dX_dpsi = 2.0 * (da1_dpsi * inv - a1 * inv2 * dd_dpsi)
        dX_dchi = 2.0 * (da1_dchi * inv - a1 * inv2 * dd_dchi)
        dY_dpsi = 2.0 * (da2_dpsi * inv - a2 * inv2 * dd_dpsi)
        dY_dchi = 2.0 * (da2_dchi * inv - a2 * inv2 * dd_dchi)
        sx = np.sqrt((dX_dpsi * pa_err) ** 2 + (dX_dchi * ell_err) ** 2)
        sy = np.sqrt((dY_dpsi * pa_err) ** 2 + (dY_dchi * ell_err) ** 2)

        return sx, sy

    def _fwd(self, longitude, colatitude):
        ph, th = longitude, colatitude

        # Build Stokes vector
        p = np.stack(
            [np.sin(th) * np.cos(ph), np.sin(th) * np.sin(ph), np.cos(th)],
            axis=-1,
        )

        d = 1.0 + p @ self.contact_point
        x = 2 * (p @ self.e1) / d
        y = 2 * (p @ self.e2) / d

        return x, y

    def _back(self, x, y):
        rho2 = x**2 + y**2

        v = (
            (1.0 - rho2) * self.contact_point[:, None]
            + 2.0 * x * self.e1[:, None]
            + 2.0 * y * self.e2[:, None]
        ) / (1.0 + rho2)

        # Cartesian -> spherical in the fixed (Q,U,V) frame
        v /= np.linalg.norm(v, axis=0)
        colatitude = np.arccos(np.clip(v[2], -1.0, 1.0))
        longitude = np.mod(np.arctan2(v[1], v[0]), 2 * np.pi)

        return np.vstack((longitude, colatitude))


@dataclass
class GnomonicProjection(PoincareProjection):
    contact_offset: float = 90

    def plot_pol_time(self, data: LightCurve):
        fig, ax = plt.subplots()
        time, P, _, _, _, _ = self._extract_sphere_points(data)

        Q, U, V = P[:, 0], P[:, 1], P[:, 2]

        gamma = 0.5 * np.arctan2(np.sqrt(U**2 + V**2), -Q)

        ax.plot(time, gamma)

        return

    def project_data(self, data: LightCurve):
        time, _, pa, ell, pa_err, ell_err = self._extract_sphere_points(data)

        longitude = 2 * pa
        colatitude = np.pi / 2 - 2 * ell

        X, Y = self._fwd(longitude, colatitude)
        Xerr, Yerr = self.project_errors(pa, ell, pa_err, ell_err)

        return time, X, Y, Xerr, Yerr

    def project_errors(self, pa, ell, pa_err, ell_err):
        c2c, s2c = np.cos(2 * ell), np.sin(2 * ell)
        c2p, s2p = np.cos(2 * pa), np.sin(2 * pa)
        Sq = c2c * c2p
        Su = c2c * s2p
        Sv = s2c
        p = np.stack([Sq, Su, Sv], axis=1)
        a1 = p @ self.e1
        a2 = p @ self.e2
        d = p @ self.contact_point
        # partials of p wrt angles
        dq_dpsi = -2 * c2c * s2p
        du_dpsi = 2 * c2c * c2p
        dv_dpsi = 0.0 * pa
        dq_dchi = -2 * s2c * c2p
        du_dchi = -2 * s2c * s2p
        dv_dchi = 2 * c2c
        dp_dpsi = np.stack([dq_dpsi, du_dpsi, dv_dpsi], axis=1)
        dp_dchi = np.stack([dq_dchi, du_dchi, dv_dchi], axis=1)
        # x = a1/d, y = a2/d
        da1_dpsi = dp_dpsi @ self.e1
        dd_dpsi = dp_dpsi @ self.contact_point
        da2_dpsi = dp_dpsi @ self.e2
        da1_dchi = dp_dchi @ self.e1
        dd_dchi = dp_dchi @ self.contact_point
        da2_dchi = dp_dchi @ self.e2
        dx_dpsi = (da1_dpsi * d - a1 * dd_dpsi) / (d * d)
        dx_dchi = (da1_dchi * d - a1 * dd_dchi) / (d * d)
        dy_dpsi = (da2_dpsi * d - a2 * dd_dpsi) / (d * d)
        dy_dchi = (da2_dchi * d - a2 * dd_dchi) / (d * d)
        sx = np.sqrt((dx_dpsi * pa_err) ** 2 + (dx_dchi * ell_err) ** 2)
        sy = np.sqrt((dy_dpsi * pa_err) ** 2 + (dy_dchi * ell_err) ** 2)

        return sx, sy

    def _set_geometry(self, P: np.ndarray):
        """
        Set projection plane contact point to mean longitude
        at equator and set up basis vectors in the tangent plane
        """

        # Compute mean PA
        mean_Pvec = np.nanmean(P, axis=0)
        mean_q, mean_u = mean_Pvec[0], mean_Pvec[1]
        self.phi_mean = np.arctan2(mean_u, mean_q)

        # Determine meridian limits
        polangle = np.arctan2(P[:, 1], P[:, 0])
        self.phi_min = np.nanmin(polangle)
        self.phi_max = np.nanmax(polangle)

        # Generate small offset from equator to ensure
        # meridians converge to one of the +V/-V poles.
        c = np.cos(np.deg2rad(self.contact_offset))
        s = np.sin(np.deg2rad(self.contact_offset))

        # Set contact point of tangent plane to mean PA and slight offset
        # from equator
        n0 = np.array([np.cos(self.phi_mean) * c, np.sin(self.phi_mean) * c, s])
        n0 /= np.linalg.norm(n0)

        # Create X, Y basis in the tangent plane
        V = np.array([0, 0, 1])
        e2 = V - (V @ n0) * n0
        e2 /= np.linalg.norm(e2)

        e1 = np.cross(e2, n0)

        self.contact_point = n0
        self.e1 = e1
        self.e2 = e2

        return

    def _fwd(self, longitude, colatitude):
        ph, th = longitude, colatitude

        # Build Stokes vector
        p = np.stack(
            [np.sin(th) * np.cos(ph), np.sin(th) * np.sin(ph), np.cos(th)],
            axis=-1,
        )

        # Project onto tangent plane
        d = p @ self.contact_point
        x = (p @ self.e1) / d
        y = (p @ self.e2) / d

        # Mask points from opposite hemisphere
        x = np.where(d > 0, x, np.nan)
        y = np.where(d > 0, y, np.nan)

        return np.vstack((x, y))

    def _back(self, x, y):
        Sv = self.contact_point[:, None] + self.e1[:, None] * x + self.e2[:, None] * y
        Sv /= np.linalg.norm(Sv, axis=0)
        colatitude = np.arccos(np.clip(Sv[2], -1.0, 1.0))
        longitude = np.arctan2(Sv[1], Sv[0]) % (2 * np.pi)

        return np.vstack((longitude, colatitude))


def plot_poincare_sphere(data: TimeFreqSeries):
    """Plot polarisation state of lightcurve / spectrum on the Poincare sphere."""

    fig = plt.figure()
    ax = fig.add_subplot(projection="3d")

    mask = get_mask(data)
    time = data.x.copy()

    polangle = np.deg2rad(data.polangle.copy())
    ellipticity = np.deg2rad(data.ellipticity.copy())
    polangle[mask] = np.nan
    ellipticity[mask] = np.nan

    q = np.cos(2 * polangle) * np.cos(2 * ellipticity)
    u = np.sin(2 * polangle) * np.cos(2 * ellipticity)
    v = np.sin(2 * ellipticity)

    ax.scatter(
        q,
        u,
        v,
        c=time,
        cmap="plasma",
        ls="-",
    )
    ax.plot(
        q,
        u,
        v,
        color="k",
        alpha=0.4,
    )

    # Plot the sphere
    u_sphere, v_sphere = np.mgrid[0 : 2 * np.pi : 100j, 0 : np.pi : 100j]
    x_sphere = np.cos(u_sphere) * np.sin(v_sphere)
    y_sphere = np.sin(u_sphere) * np.sin(v_sphere)
    z_sphere = np.cos(v_sphere)
    ax.plot_surface(
        x_sphere,
        y_sphere,
        z_sphere,
        color="lightgray",
        alpha=0.1,
    )

    # Plot the equator
    theta = np.linspace(0, 2 * np.pi, 200)
    q_eq = np.cos(theta)
    u_eq = np.sin(theta)
    v_eq = np.zeros_like(theta)

    ax.plot(
        q_eq,
        u_eq,
        v_eq,
        color="gray",
        linestyle="--",
        linewidth=1,
    )

    ax.set_xlim(-1, 1)
    ax.set_ylim(-1, 1)
    ax.set_zlim(-1, 1)
    ax.set_box_aspect([1, -1, 1])

    ax.set_xticks([-1, -0.5, 0, 0.5, 1])
    ax.set_yticks([-1, -0.5, 0, 0.5, 1])
    ax.set_zticks([-1, -0.5, 0, 0.5, 1])

    ax.set_xlabel(r"$Q/P$")
    ax.set_ylabel(r"$U/P$")
    ax.set_zlabel(r"$V/P$")

    # Color map for time evolution
    norm = mcolors.Normalize(vmin=np.nanmin(time), vmax=np.nanmax(time))
    sm = cm.ScalarMappable(cmap="plasma", norm=norm)
    sm.set_array([])
    fig.colorbar(sm, ax=ax, label="Phase")

    return
