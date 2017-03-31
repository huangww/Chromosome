# -*- coding: utf-8 -*-
"""
    chromosome.stretch_coil
    ~~~~~~~~~~~~~~~~~~~~~~~

    This script plot the figure of stretch-coil transition of a piined polymer loop

    :copyright: (c) 2017 by YOUR_NAME.
    :license: LICENSE_NAME, see LICENSE for more details.
"""

import numpy as np
import matplotlib.pyplot as plt
import theory as th


def stretch_coil_transition():
    """ Plot the figure of stretch-coil transition
    :returns: TODO

    """
    fig = plt.figure(0, figsize=(5, 4))
    plt.rc('text', usetex=True)
    font = {'family' : 'sans-serif',
            'serif'  : 'Helvetica',
            'weight' : 'normal',
            'size'   : 14}
    plt.rc('lines', lw=2)
    plt.rc('font', **font)
    ax = fig.add_subplot(111)
    fig.subplots_adjust(left=0.16, right=0.95,\
        bottom=0.12, top=0.95, wspace=0.25)
    param_range = [0.1, 0.5, 1]
    for param in param_range:
        fname = 'stretch_coil_N100_F%g_T1.dat' % param
        data = np.loadtxt(fname)
        label_string = r'$F=%g$' % param
        ax.plot(data[:, 0], data[:, 1]/data[:, 1].max(), label=label_string)
    ax.set_yscale('log')
    ax.set_xlim([0, 1000])
    ax.set_ylim([0.01, 1])
    ax.legend()
    ax.set_xlabel(r'$t$')
    ax.set_ylabel(r'$z_d(t) / z_d(0)$')

    ax_inset = plt.axes([0.31, 0.23, 0.3, 0.3])
    ax_inset.plot(data[:, 0], data[:, 1].max()-data[:, 1], 'r-')
    ax_inset.plot(data[:, 0], data[:, 0], 'k--')
    ax_inset.set_xscale('log')
    ax_inset.set_yscale('log')
    ax_inset.set_ylim([0.1, 100])
    ax_inset.set_xlabel(r'$t$', labelpad=0)
    ax_inset.set_ylabel(r'$z_d(0) - z_d(t)$', labelpad=0)

    plt.show()
    return


def coil_stretch_transition():
    """ Plot the figure of coil-stretch transition
    :returns: TODO

    """
    fig = plt.figure(0, figsize=(5, 4))
    plt.rc('text', usetex=True)
    font = {'family' : 'sans-serif',
            'serif'  : 'Helvetica',
            'weight' : 'normal',
            'size'   : 14}
    plt.rc('lines', lw=2)
    plt.rc('font', **font)
    ax = fig.add_subplot(111)
    fig.subplots_adjust(left=0.16, right=0.95,\
        bottom=0.12, top=0.95, wspace=0.25)
    param_range = [0.1, 0.5, 1]
    for param in param_range:
        fname = 'coil_stretch_N100_F%g_T1.dat' % param
        data = np.loadtxt(fname)
        label_string = r'$F=%g$' % param
        ax.plot(data[:, 0], data[:, 1]/data[:, 1].max(), label=label_string)
    # ax.set_yscale('log')
    ax.set_xlim([0, 1000])
    ax.set_ylim([0.01, 1])
    ax.legend()
    ax.set_xlabel(r'$t$')
    ax.set_ylabel(r'$z_d(t) / z_d(\infty)$')

    ax_inset = plt.axes([0.46, 0.24, 0.42, 0.38])
    data = np.loadtxt('slope_coil_stretch_N100_T1.dat')
    # tau_theory = th.relaxation_time_3d(data[:, 0], 1.0, 100)
    # tau_theory = 1. / data[:, 0]
    # factor = tau_theory.max()
    # tau_theory = tau_theory / factor
    # ax_inset.plot(data[:, 0], tau_theory, 'k--')
    ax_inset.plot(data[:, 0], data[:, 1], 'o')
    coefs = np.polyfit(data[:, 0], data[:, 1], 1)
    ax_inset.plot(data[:, 0], coefs[0]*data[:, 0] + coefs[1], 'k--')
    # ax_inset.set_xscale('log')
    # ax_inset.set_yscale('log')
    # ax_inset.set_xlim([0.1, 1])
    # ax_inset.set_ylim([0, 1])
    ax_inset.set_xlabel(r'$F$', labelpad=0)
    ax_inset.set_ylabel(r'$v_{\parallel}$', labelpad=0)

    plt.show()
    return


def main():
    """ The main function
    :returns: TODO

    """
    # stretch_coil_transition()
    coil_stretch_transition()

if __name__ == "__main__":
    main()
