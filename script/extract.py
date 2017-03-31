"""
This script treat the trajectory data and extract information such like relaxation time scale
"""
# pylint: disable=E0401
# pylint: disable=E1101
import glob
import os
import numpy as np
import statsmodels.tsa.stattools as ss
# import matplotlib.pyplot as plt
import fitting as af
import theory as th


def get_mean_acf(file_list):
    """ Get the mean auto correlation function from a list of trajectory data

    :file_list: file name list of trajectory data
    :returns: mean_acf, 1d array

    """
    acf_list = []
    for fname in file_list:
        data = np.loadtxt(fname)
        x = data[:, 1]
        acf = ss.acf(x, nlags=len(x), fft=True)
        acf_list.append(acf)
    min_len_acf = min([len(acf) for acf in acf_list])
    mean_acf = np.zeros(min_len_acf)
    for acf in acf_list:
        mean_acf = mean_acf + acf[:min_len_acf]
    mean_acf = mean_acf / len(file_list)
    return mean_acf


def extract_tau_from_data(data):
    """ Extract tau from data of trajectory

    :data: input data, the first column is time
    :returns: a list of tau for each column

    """
    t = data[:, 0]
    tau = []
    for i in range(1, len(data[0])):
        x = data[:, i]
        tau.append(af.fit_tau_acf(t, x))
    return tau


def mean_data(file_list):
    """ Load the data in file list and calculate the mean among them

    :file_list: a list of data file name
    :returns: ndarray of merged data

    """
    file_size = np.array([os.path.getsize(fname) for fname in file_list])
    max_size_index = file_size.argmax()
    data0 = np.loadtxt(file_list[max_size_index])
    data_shape = np.shape(data0)
    data_mean = np.zeros(data_shape)
    count = np.zeros(data_shape)
    for fname in file_list:
        data = np.loadtxt(fname)
        if np.shape(data) == data_shape:
            count = count + 1
            data_mean = data_mean + data
        else:
            for i, x in enumerate(data):
                data_mean[i] = data_mean[i] + x
                count[i] = count[i] + 1
    data_mean = data_mean / count
    return data_mean


def get_tau_from_mean_data():
    """ Extract time scale time from averaged data
    :returns: TODO

    """
    # para_range = np.linspace(0.6, 1, 1)
    para_range = np.logspace(-3, 0)
    slope_data = []
    for para in para_range:
        print para,
        fname = 'coil_stretch_N100/BeadRod_rd_N100_F%g_T1_*.dat'%para
        file_list = [f for f in glob.iglob(fname) if os.path.getsize(f) > 1e5]
        if file_list:
            data = mean_data(file_list)
            x = data[:, 1]/data[:, 1].max()
            end = next((i for i, v in enumerate(x) if abs(v-1) < 1e-2), len(x)-1)
            coefs = np.polyfit(data[:end/2, 0], x[:end/2], 1)
            slope = coefs[0]
            slope_data.append([para, slope])
            print slope
    fname = 'data/slope_coil_stretch_N100_T1.dat'
    np.savetxt(fname, slope_data)
    return


def get_tau_from_data():
    """ Extract time scale tau from the raw trajectory data
    :returns: TODO

    """
    tau_data = []
    # para_range = np.linspace(0.1, 1.0, 10)
    para_range = np.logspace(-3, 0)
    for para in para_range:
        # fname = 'stretch_coil_N100/BeadRod_rd_N100_F%g_T1_*.dat'%para
        fname = 'coil_stretch_N100/BeadRod_rd_N100_F%g_T1_*.dat'%para
        file_list = [f for f in glob.iglob(fname) if os.path.getsize(f) > 1e5]
        if file_list:
            tau_list = []
            for fname in file_list:
                print fname,
                data = np.loadtxt(fname)
                # tau = af.fit_tau_exp(data[:,0], data[:,1])
                zd_theory = th.meanz(100, 1./para)[50]
                x = zd_theory - data[:, 1]
                # x = data[:, 1]
                end = next((i for i, v in enumerate(x) if v < 1e-6), len(x)-1)
                tau = data[end, 0]
                # tau = extract_tau_from_data(data)
                # tau = np.nan if tau > 1e3 else tau
                # tau[tau > 1e3] = np.nan
                # tau[tau < 0] = np.nan
                tau_list.append(tau)
                print tau
            tau = np.nanmean(tau_list, axis=0)
            tau_data.append([para, tau])
            # tau = np.insert(tau, 0, para)
            # tau_data.append(tau.tolist())
    # fname = 'data/tau_stretch_coil_N100_T1.dat'
    fname = 'data/tau_coi_stretch_N100_T1.dat'
    np.savetxt(fname, tau_data)
    return


def main():
    """ The main function
    :returns: TODO

    """
    # get_tau_from_data()
    get_tau_from_mean_data()


if __name__ == "__main__":
    main()
