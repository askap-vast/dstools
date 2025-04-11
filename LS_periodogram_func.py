from scipy.interpolate import interp1d
from astropy.timeseries import LombScargle
import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import find_peaks, peak_widths
from mpl_toolkits.axes_grid1.inset_locator import inset_axes

def calc_lomb_scargle_raw(t,y):
    """
    Calculate the Lomb-Scargle periodogram of a time series.
    Parameters
    ----------
    t : array-like
        Time values.
    y : array-like
        Flux values.
    Returns
    -------
    freq : array-like
        Frequencies of the periodogram.
    amp : array-like
        Amplitude of the periodogram.
    """
    oversample = 10
    tmax = t.max()
    tmin = t.min()
    df = 1.0 / (tmax - tmin)
    fmin = df
    fmax = 100 # highest frequency (in c/d)
  
    freq = np.arange(fmin, fmax, df / oversample)
    model = LombScargle(t, y)
    sc = model.power(freq, method="fastnifty", normalization="psd")

    fct = np.sqrt(4./len(t))
    amp = np.sqrt(sc) * fct
    return freq, amp

def index_to_xdata(xdata, indices):
    "interpolate the values from signal.peak_widths to xdata"
    ind = np.arange(len(xdata))
    f = interp1d(ind,xdata)
    return f(indices)

def periodogram(df,sig=4, idx = 0,n_peaks=1,opt_period=1.3,opt_e_period=0.05,opt_off=False, inset_sigma=2,save=False):
    """
    Calculate the Lomb-Scargle periodogram of a time series and plot the results.
    Parameters
    ----------
    df : pandas.DataFrame
        DataFrame containing the time series data.
    sig : float
        Significance level for peak detection.
    idx : int
        Index of the peak to highlight in the plot.
    n_peaks : int
        Number of peaks to detect.
    opt_period : float
        Period of the optical telescope in hours for comparison.
    opt_e_period : float
        Uncertainty in the optical period in hours.
    opt_off : bool
        If True, do not plot the optical period.
    inset_sigma : float
        Sigma level for the inset plot.
    save : bool
        If True, save the plot as a PDF file.
    Returns
    -------
    None
    """
    y_i = df['y_i']
    t = df['t_mjd']
    temp,amp = calc_lomb_scargle_raw(t, y_i)
    x = temp/24
    y = amp

    peaks, _ = find_peaks(y,height=sig*np.std(y))
    widths, width_heights, left_ips, right_ips = peak_widths(y, peaks, rel_height=0.5)
    widths = index_to_xdata(x, widths)
    left_ips = index_to_xdata(x, left_ips)
    right_ips = index_to_xdata(x, right_ips)
    diffs=np.abs(right_ips-left_ips)
    hwhm=0.5*diffs

    periods=[]
    e_periods=[]
    for i in range(n_peaks):
        freq=x[peaks][i]
        e_freq = hwhm[i]
        print("\t\tf={:.4f}+/-{:.4f}".format(freq,e_freq))
        period = 1/freq
        e_period = (freq**-2)*e_freq
        periods.append(period)
        e_periods.append(e_period)
        print("\t\tP={:.4f}+/-{:.4f}\n".format(period,e_period))

    fig, ax = plt.subplots(figsize=(10, 6))

    ax.plot(x, y, label='Stokes I LS')
    ax.axhline(np.mean(y)+sig*np.std(y), color='gray', ls='--', label=r'{}$\sigma$'.format(sig))
    ax.plot(x[peaks][idx], y[peaks][idx], "x", label='Peak',c='k')
    ax.hlines(width_heights[idx], left_ips[idx], right_ips[idx], color='r', label='FWHM')

    if opt_off == False:
        ax.axvline(x=1/opt_period, color="darkorange", ls=":", label=r"P$_{\rm{OPT}}=$"+"{}+/-{} hr".format(opt_period,opt_e_period))
        ax.axvspan(xmin=(1/(opt_period-opt_e_period)), xmax=(1/(opt_period+opt_e_period)), alpha=0.2, color='darkorange')
    ax.axvline(x=1/periods[idx], color="darkgreen", ls="--", label=r"P$_{\rm{LS}}=$"+"{:.3f}+/-{:.3f} hr".format(periods[idx],e_periods[idx]))
    ax.axvspan(xmin=(1/(periods[idx]+e_periods[idx])), xmax=(1/(periods[idx]-e_periods[idx])), alpha=0.2, color='darkgreen') 

    ax.set_xlabel('Freq [hour$^{-1}$]',fontsize=20)
    ax.set_ylabel('Power',fontsize=20)
    ax.legend()

    ax_inset = inset_axes(ax, width="30%", height="40%", loc="upper center")
    ax_inset.plot(x, y, label='Stokes I LS')
    ax_inset.plot(x[peaks][idx], y[peaks][idx], "x", c='k')
    ax_inset.hlines(width_heights[idx], left_ips[idx], right_ips[idx], color='r')
    if opt_off == False:
        ax_inset.axvline(x=1/opt_period, color="darkorange", ls=":")
        ax_inset.axvspan(xmin=(1/(opt_period-opt_e_period)), xmax=(1/(opt_period+opt_e_period)), alpha=0.2, color='darkorange')
    ax_inset.axvline(x=1/periods[idx], color="darkgreen", ls="--")
    ax_inset.axvspan(xmin=(1/(periods[idx]+e_periods[idx])), xmax=(1/(periods[idx]-e_periods[idx])), alpha=0.2, color='darkgreen')
    ax_inset.set_xlim((1/(periods[idx]+inset_sigma*e_periods[idx])), (1/(periods[idx]-inset_sigma*e_periods[idx])))
    ax_inset.set_ylim(0.0, 1.10*y.max())
    ax_inset.grid(True)
    ax.indicate_inset_zoom(ax_inset, edgecolor="black")
    ax.set_ylim(0.0, 1.30*y.max())
    ax.set_xlim(x.min(), x.max())
    if save == True:
        fig.savefig('LS_periodogram_{:.2f}.pdf'.format(t[0]),dpi=300)

