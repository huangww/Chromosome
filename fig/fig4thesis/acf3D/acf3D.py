#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np
import matplotlib.pyplot as plt
import acfFit as af

Teff = [10, 20, 50]

def plotFig():
    fig = plt.figure(0,figsize=(5,4))
    plt.rc('text',usetex = True)
    font = {'family' : 'sans-serif',
            'serif'  : 'Helvetica',
            'weight' : 'normal',
            'size'   : 16 }
    plt.rc('lines', lw=2)
    plt.rc('font', **font)
    fig.subplots_adjust(left=0.15, right =0.95,\
        bottom=0.15, top =0.95, wspace=0.25)
    ax = fig.add_subplot(111)

    fext_list = [0.01, 0.05, 0.1]
    fit_range = [(200, 1000), (200, 600), (100, 300)]
    for fext, fit in zip(fext_list, fit_range):
        fileName = 'acf_N100_F%g_T1.dat' % fext
        acf = np.loadtxt(fileName)
        t = np.linspace(0, len(acf), len(acf))
        labelString = r"$F=%g$" % fext
        ax.plot(t, acf, label=labelString)
        start, end = fit
        coefs = np.polyfit(t[start:end], np.log(acf[start:end]), 1)
        ax.plot(t, np.exp(coefs[0]*t+coefs[1]), 'k--',lw=1)

    ax.legend(fontsize=16)
    ax.set_yscale('log')
    ax.set_xlabel(r'$t$', labelpad=0)
    ax.set_ylabel(r'$\left<z_d(t)z_d(0)\right>$',labelpad=0)
    ax.set_xlim([0, 1e3])
    # ax.set_xticks([0, 1e3, 2e3, 3e3, 4e3])
    # ax.ticklabel_format(axis='x', style='sci',scilimits=(0,1))
    # ax.tick_params(axis='both', which='major', pad=1)
    ax.set_ylim([1e-2, 1])
    plt.show()

def main():
    plotFig()

if __name__ == "__main__":
    main()
        





