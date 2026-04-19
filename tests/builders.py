import astropy.units as u
import numpy as np

from dstools.dynamic_spectrum import DynamicSpectrum, LightCurve, Spectrum


def _make_dynamic_spectrum_instance(**overrides):
    ds = object.__new__(DynamicSpectrum)

    defaults = {
        "time": np.linspace(0, 99, 100),
        "freq": np.linspace(500, 1000, 50),
        "tunit": u.s,
        "period": 10.0,
        "period_offset": 0.0,
        "phase_bins": 10,
        "fold_periods": 1,
        "corr_dumptime": 10.1,
        "header": {
            "time_start": "2020-01-01 00:00:00",
            "time_scale": "utc",
            "phasecentre": "00:00:00 +00:00:00",
            "telescope": "ATCA",
            "feeds": "linear",
            "channels": 50,
        },
        "DM": 0,
        "DM_reffreq": None,
        "flag_channels": None,
        "flag_times": None,
        "absolute_times": True,
        "fold": False,
        "trim": False,
        "tavg": 1,
        "favg": 1,
    }

    for key, value in defaults.items():
        setattr(ds, key, value)

    for key, value in overrides.items():
        setattr(ds, key, value)

    return ds


def make_synthetic_dynamic_spectrum(
    absolute_times: bool = True,
    fold: bool = False,
    tmax: int = 6,
    shape: tuple[int, int] = (5, 10),
    **overrides,
):
    tunit = u.hour
    fmin, fmax = 100, 105
    fold_periods = 1
    time = np.linspace(0, tmax, shape[0])
    freq = np.linspace(fmin, fmax, shape[1])

    header = {
        "time_start": "2020-01-01 20:00:00",
        "time_scale": "utc",
    }

    rng = np.random.default_rng()
    data = {
        "I": rng.normal(0, 1, shape) + 1j * rng.normal(0, 1, shape),
        "Q": rng.normal(0, 1, shape) + 1j * rng.normal(0, 1, shape),
        "U": rng.normal(0, 1, shape) + 1j * rng.normal(0, 1, shape),
        "V": rng.normal(0, 1, shape) + 1j * rng.normal(0, 1, shape),
        "L": np.abs(rng.normal(0, 1, shape)),
    }

    label = "Time (s)" if not fold else "Phase (deg)"

    defaults = {
        "data": data,
        "tunit": tunit,
        "fmin": fmin,
        "fmax": fmax,
        "tmin": 0,
        "tmax": tmax,
        "time": time,
        "freq": freq,
        "fold": fold,
        "fold_periods": fold_periods,
        "phase_bins": shape[0],
        "absolute_times": absolute_times,
        "header": header,
        "_timelabel": label,
    }
    defaults.update(overrides)

    return _make_dynamic_spectrum_instance(**defaults)


def make_stokes_data(shape: tuple[int, int] = (5, 10)):
    i_real = np.full(shape, 10.0)
    q_real = np.full(shape, 3.0)
    u_real = np.full(shape, 4.0)
    v_real = np.full(shape, 1.0)

    return {
        "I": i_real + 1j * np.zeros(shape),
        "Q": q_real + 1j * np.zeros(shape),
        "U": u_real + 1j * np.zeros(shape),
        "V": v_real + 1j * np.zeros(shape),
        "L": np.sqrt(q_real**2 + u_real**2) + 1j * np.zeros(shape),
    }


def make_synthetic_dynamic_spectrum_with_stokes_data(
    shape: tuple[int, int] = (5, 10),
    **overrides,
):
    time = np.linspace(0, shape[0] - 1, shape[0])
    freq = np.linspace(100, 100 + shape[1] - 1, shape[1])

    ds = make_synthetic_dynamic_spectrum(
        absolute_times=overrides.pop("absolute_times", True),
        fold=overrides.pop("fold", False),
        shape=shape,
        time=time,
        freq=freq,
        data=make_stokes_data(shape=shape),
        tmin=time[0],
        tmax=time[-1],
        fmin=freq[0],
        fmax=freq[-1],
        **overrides,
    )

    return ds


def make_synthetic_lightcurve(
    absolute_times: bool = True,
    fold: bool = False,
    imag: bool = False,
    pol_sigma: float = 4,
    shape: tuple[int, int] = (5, 10),
    **overrides,
):
    ds = make_synthetic_dynamic_spectrum_with_stokes_data(
        absolute_times=absolute_times,
        fold=fold,
        shape=shape,
        **overrides,
    )

    return LightCurve(ds, imag=imag, pol_sigma=pol_sigma)


def make_synthetic_spectrum(
    absolute_times: bool = True,
    fold: bool = False,
    imag: bool = False,
    pol_sigma: float = 4,
    shape: tuple[int, int] = (5, 10),
    **overrides,
):
    ds = make_synthetic_dynamic_spectrum_with_stokes_data(
        absolute_times=absolute_times,
        fold=fold,
        shape=shape,
        **overrides,
    )

    return Spectrum(ds, imag=imag, pol_sigma=pol_sigma)


def make_pulse_array(period, tres, tsamples: int = 1000, channels: int = 50):
    array = np.zeros((tsamples, channels), dtype=np.complex128)

    sigma = 90
    pulse_halfwidth = sigma // tres
    x = np.tile(np.arange(-pulse_halfwidth, pulse_halfwidth).reshape(-1, 1), channels)
    pulse = 100 * np.exp(-(x**2) / (2 * pulse_halfwidth**2)) + 1j * 0 * x

    for pulse_index in range(10, tsamples - 10, period // tres):
        array[pulse_index - pulse_halfwidth : pulse_index + pulse_halfwidth, :] += pulse

    time = np.arange(0, tres * tsamples, tres)

    return array, time
