"""
Script used for fitting functions, auto correlation functions
"""
import numpy as np
import statsmodels.tsa.stattools as ss
from scipy.signal import savgol_filter
import matplotlib.pyplot as plt
# pylint: disable=E1101


def acf(x):
    """ get the auto correlation function of a signal

    :x: the input signal, 1D array
    :returns: array of auto correlation function

    """
    _acf = ss.acf(x, nlags=len(x), fft=True)
    return _acf


def find_fit_range_acf(x):
    """ Find the fitting range of an auto correlation function

    :x: input auto correlation function, 1D array
    :returns: index of the fitting rage (start, end)

    """
    end = next((i for i, v in enumerate(x) if v < 1e-6), len(x)-1)
    x = x[:end]
    x = x / x.max()
    interval = len(x)/300
    inter = interval if interval > 1 else 1
    x_log = np.log(x[::inter])
    dx1 = np.maximum(0, savgol_filter(x_log, 5, 2, deriv=1))
    dx1sum = [dx1[:i+1].sum() for i in range(len(dx1))]
    end = next((i for i, v in enumerate(dx1sum) if v > 0.1), len(dx1sum)-1)
    x = x[:max(end*inter, 30)]
    interval = len(x)/300
    inter = interval if interval > 1 else 1
    x_log = np.log(x[::inter])
    dx2abs = np.abs(savgol_filter(x_log, 5, 2, deriv=2))
    cutoff = len(x_log)/5
    fit_len = len(x_log)/3
    d2sum = [dx2abs[i:i+fit_len].sum() for i in \
            range(cutoff, len(x)-cutoff-fit_len)]
    start = (np.array(d2sum).argmin() + cutoff) * inter
    end = start + fit_len * inter
    return start, end


def fit_tau_acf(t, x):
    """ Fit the time scale tau by fitting the auto correlation function of trajectory x(t)

    :t: time, 1D array
    :x: 1D trajectory, should be a stationary process
    :returns: tau, float

    """
    _acf = acf(x)
    start, end = find_fit_range_acf(_acf)
    t_fit = t[start:end]
    x_fit = _acf[start:end]
    coefs = np.polyfit(t_fit, np.log(x_fit), 1)
    tau = -1./coefs[0]
    if tau > 1e6 or tau < 0:
        tau = np.nan
    return tau


def find_fit_range_exp_grow(x):
    """ Find the fitting range of an growing exponential function

    :x: input exponential function, 1d array
    :returns: index of the fitting rage (start, end)

    """
    x = (x.max() - x) / x.max()
    end = next((i for i, v in enumerate(x) if v < 1e-3), len(x)-1)
    x = x[:end]
    interval = len(x)/300
    inter = interval if interval > 1 else 1
    x_log = np.log(x[::inter])
    dx2abs = np.abs(savgol_filter(x_log, 5, 2, deriv=2))
    end = next((i for i, v in enumerate(dx2abs) if v > 0.1), len(dx2abs)-1)
    end = end * inter
    fit_len = end/3
    start = end - fit_len
    return start, end


def fit_tau_exp_grow(t, x):
    """ Fit time scale tau from a growing exponential function

    :t: time, 1d array
    :x: x(t) is the growing exponential function, 1d array
    :returns: tau, float

    """
    start, end = find_fit_range_exp(x)
    t_fit = t[start:end]
    x = (x.max() - x) / x.max()
    x_fit = x[start:end]
    coefs = np.polyfit(t_fit, np.log(x_fit), 1)
    tau = -1./coefs[0]
    if tau > 1e6 or tau < 0:
        tau = np.nan
    return tau


def find_fit_range_exp(x):
    """ Find the fitting range of an exponential function

    :x: input exponential function, 1D array
    :returns: index of the fitting rage (start, end)

    """
    x = x / x.max()
    start = 0
    end = next((i for i, v in enumerate(x) if v < 1e-6), len(x)-1)
    x = x[:end]
    interval = len(x)/300
    inter = interval if interval > 1 else 1
    x_log = np.log(x[::inter])
    dx2abs = np.abs(savgol_filter(x_log, 5, 2, deriv=2))
    end = next((i for i, v in enumerate(dx2abs) if v > 0.1), len(dx2abs)-1)
    end = inter * end
    return start, end


def fit_tau_exp(t, x):
    """ Fit time scale tau from a exponential function x(t)

    :t: time, 1D array
    :x: exponential function data, 1D array
    :returns: tau, float

    """
    start, end = find_fit_range_exp(x)
    t_fit = t[start:end]
    x_fit = x[start:end]
    coefs = np.polyfit(t_fit, np.log(x_fit), 1)
    tau = -1./coefs[0]
    if tau > 1e6 or tau < 0:
        tau = np.nan
    return tau


def show_fitting_exp(t, x):
    """ Plot and show the fitting figure

    :t: t of function x(t)
    :x: x of function x(t)
    :returns: TODO

    """
    fig = plt.figure(0)
    ax = fig.add_subplot(111)
    x = x / x.max()
    end = next((i for i, v in enumerate(x) if v < 1e-6), len(x)-1)
    t, x = t[:end], x[:end]
    ax.plot(t, x)
    interval = len(x) / 300
    inter = interval if interval > 1 else 1
    x_log = np.log(x[::inter])
    dx2abs = np.abs(savgol_filter(x_log, 5, 2, deriv=2))
    # ax_twinx = ax.twinx()
    # ax_twinx.plot(t[::inter], dx2abs, 'ro-')
    start, end = find_fit_range_exp(x)
    # start, end = 0, 1000
    ax.fill_betweenx(x, t[start], t[end], alpha=0.1)
    t_fit = t[start:end]
    x_fit = x[start:end]
    coefs = np.polyfit(t_fit, np.log(x_fit), 1)
    x_fit = np.exp(coefs[0]*t+coefs[1])
    ax.plot(t, x_fit)
    ax.set_xlim([0, 5*t[end]])
    ax.set_ylim([0.01, 1])
    # ax.set_yscale('log')
    ax.set_xlabel(r'$t$')
    ax.set_ylabel(r'$z_d$')
    print fit_tau_exp(t, x)
    plt.show()


def show_fitting_exp_grow(t, x):
    """ Plot and show the fitting figure

    :t: t of function x(t)
    :x: x of function x(t)
    :returns: TODO

    """
    fig = plt.figure(0)
    ax = fig.add_subplot(111)
    start, end = find_fit_range_exp_grow(x)
    x = (x.max() - x) / x.max()
    ax.plot(t, x)
    interval = len(x) / 300
    inter = interval if interval > 1 else 1
    x_log = np.log(x[::inter])
    dx2abs = np.abs(savgol_filter(x_log, 5, 2, deriv=2))
    ax_twinx = ax.twinx()
    ax_twinx.plot(t[::inter], dx2abs, 'ro-')
    ax.fill_betweenx(x, t[start], t[end], alpha=0.1)
    t_fit = t[start:end]
    x_fit = x[start:end]
    coefs = np.polyfit(t_fit, np.log(x_fit), 1)
    x_fit = np.exp(coefs[0]*t+coefs[1])
    ax.plot(t, x_fit)
    ax.set_xlim([0, 2*t[end]])
    ax.set_ylim([0.01, 1])
    ax.set_yscale('log')
    ax.set_xlabel(r'$t$')
    ax.set_ylabel(r'$z_d$')
    plt.show()


def main():
    """ The main function
    :returns: TODO

    """
    print "The main function here!"
    # fname = 'stretch_coil_N100/BeadRod_rd_N100_F0.1_T1_19.dat'
    # fname = 'stretch_coil_N100/BeadRod_rd_N100_F0.1_T1_88.dat'
    # fname = 'coil_stretch_N100/BeadRod_rd_N100_F0.1_T1_88.dat'
    # fname = 'data/coil_stretch_N100_F0.1_T1.dat'
    fname = 'data/coil_stretch_N100_F1_T1.dat'
    data = np.loadtxt(fname)
    t, x = data[:, 0], data[:, 1]
    show_fitting_exp_grow(t, x)
    return


if __name__ == "__main__":
    main()
