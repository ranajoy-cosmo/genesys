#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
import os
import sys
from simulation.global_config import global_paths

spectra_folder = os.path.join(global_paths.base_dir, "spectra")

def plot_theoretical(lmax, plot_log=True, plot_list=["TT", "EE", "BB"], lensed=True):

    spectra_fid_r_0001_lensed = np.load(os.path.join(spectra_folder, "r_0001/lensedtot_cls.npy"))[..., 2:lmax+1]
    spectra_fid_r_001_lensed = np.load(os.path.join(spectra_folder, "r_001/lensedtot_cls.npy"))[..., 2:lmax+1]
    spectra_fid_r_01_lensed = np.load(os.path.join(spectra_folder, "r_01/lensedtot_cls.npy"))[..., 2:lmax+1]
    spectra_fid_r_0_lensed = np.load(os.path.join(spectra_folder, "r_0/lensedtot_cls.npy"))[..., 2:lmax+1]
    spectra_fid_r_0001_unlensed = np.load(os.path.join(spectra_folder, "r_0001/unlensed_cls.npy"))[..., 2:lmax+1]
    spectra_fid_r_001_unlensed = np.load(os.path.join(spectra_folder, "r_001/unlensed_cls.npy"))[..., 2:lmax+1]
    spectra_fid_r_01_unlensed = np.load(os.path.join(spectra_folder, "r_01/unlensed_cls.npy"))[..., 2:lmax+1]
    spectra_fid_r_0_unlensed = np.load(os.path.join(spectra_folder, "r_0/unlensed_cls.npy"))[..., 2:lmax+1]

    ell = np.arange(lmax+1)[2:]

    if plot_log:
        plotter = plt.loglog
    else:
        plotter = plt.plot

    #Plotting the TT and EE spectra
    
    if "TT" in plot_list:
        if lensed:
            plotter(ell, np.sqrt(ell*(ell+1)*spectra_fid_r_0_lensed[0]/2/np.pi), color='k')
        else:
            plotter(ell, np.sqrt(ell*(ell+1)*spectra_fid_r_0_unlensed[0]/2/np.pi), color='k')
        plt.annotate("TT", xy=(1.25, 31), size=12)

    if "EE" in plot_list:
        if lensed:
            plotter(ell, np.sqrt(ell*(ell+1)*spectra_fid_r_0_lensed[1]/2/np.pi), color='k')
        else:
            plotter(ell, np.sqrt(ell*(ell+1)*spectra_fid_r_0_unlensed[1]/2/np.pi), color='k')
        plt.annotate("EE", xy=(1.25, 0.2), size=12)


    if "BB" in plot_list:
        plotter(ell, np.sqrt(ell*(ell+1)*spectra_fid_r_01_unlensed[2]/2/np.pi), color='k')

        plotter(ell, np.sqrt(ell*(ell+1)*spectra_fid_r_001_unlensed[2]/2/np.pi), color='k')

        plotter(ell, np.sqrt(ell*(ell+1)*spectra_fid_r_0001_unlensed[2]/2/np.pi), color='k')

        plotter(ell, np.sqrt(ell*(ell+1)*spectra_fid_r_01_lensed[2]/2/np.pi), color='dimgray')

        plotter(ell, np.sqrt(ell*(ell+1)*spectra_fid_r_001_lensed[2]/2/np.pi), color='dimgray')

        plotter(ell, np.sqrt(ell*(ell+1)*spectra_fid_r_0001_lensed[2]/2/np.pi), color='dimgray')

        plotter(ell, np.sqrt(ell*(ell+1)*spectra_fid_r_0_lensed[2]/2/np.pi), color='k')

        plt.annotate("BB", xy=(1.25, 0.003), size=12)

        plt.annotate("r = 0.1", xy=(1.25, 0.062), size=10)

        plt.annotate("r = 0.01", xy=(1.25, 0.019), size=10)

        plt.annotate("r = 0.001", xy=(1.25, 0.0065), size=10)

    plt.show()


def make_decorations(ylim=[1e-4,100], xlim=[1,3000], leg_loc=None, unit="uK^2"):
    plt.ylim(ylim)
    plt.xlim(xlim)
    plt.xlabel('$l$', fontsize=12)
    plt.ylabel('$\sqrt{l(l+1)C_l/2\pi}$ $[\mu K]$', fontsize=12)
    if leg_loc == None:
        leg_loc = "upper left"
    plt.legend(loc=leg_loc , prop={'size':12})


def plot_spectra(cl, lmax=None, label=None, plot_log=True, color=None):
    if plot_log:
        plotter = plt.loglog
    else:
        plotter = plt.plot

    if lmax==None:
        lmax = cl.size - 1
    ell = np.arange(lmax+1)[2:]
    if label:
        plotter(ell, np.sqrt(ell*(ell+1)*cl[2:lmax+1]/2/np.pi), label=label, color=color)
        if color:
            plotter(ell, np.sqrt(ell*(ell+1)*cl[2:lmax+1]/2/np.pi), label=label, color=color)
        else:
            plotter(ell, np.sqrt(ell*(ell+1)*cl[2:lmax+1]/2/np.pi), label=label)
    else:
        if color:
            plotter(ell, np.sqrt(ell*(ell+1)*cl[2:lmax+1]/2/np.pi), color=color)
        else:
            plotter(ell, np.sqrt(ell*(ell+1)*cl[2:lmax+1]/2/np.pi))


def plot_error_margins(cl_mean, cl_std, lmax=None, color=None):

    if lmax==None:
        lmax = cl_mean.size - 1
    ell = np.arange(lmax+1)[2:]

    if color:
        plt.fill_between(ell, np.sqrt(ell*(ell+1)*(cl_mean[2:lmax+1] - cl_std[2:lmax+1])/2/np.pi), np.sqrt(ell*(ell+1)*(cl_mean[2:lmax+1] + cl_std[2:lmax+1])/2/np.pi), facecolor=color, alpha=0.3, linewidth=0)
    else:
        plt.fill_between(ell, np.sqrt(ell*(ell+1)*(cl_mean[2:lmax+1] - cl_std[2:lmax+1])/2/np.pi), np.sqrt(ell*(ell+1)*(cl_mean[2:lmax+1] + cl_std[2:lmax+1])/2/np.pi), alpha=0.3, linewidth=0)


def plot_error_bars(cl, sigma_cl, lmax=None, color='r', alpha=0.5):
    if lmax==None:
        lmax = cl.size - 1
    ell = np.arange(lmax+1)
    upper = np.sqrt(ell*(ell+1)*(cl + sigma_cl)[:lmax+1]/2/np.pi)
    lower = np.sqrt(ell*(ell+1)*(cl - sigma_cl)[:lmax+1]/2/np.pi)
    plt.fill_between(ell, upper, lower, facecolor=color, alpha=alpha, lw=0)


def plot_binned_spectra(cl, ell, label=None, plot_log=True):
    if plot_log:
        plotter = plt.loglog
    else:
        plotter = plt.plot

    if label:
        plotter(ell, np.sqrt(ell*(ell+1)*cl/2/np.pi), label=label)
    else:
        plotter(ell, np.sqrt(ell*(ell+1)*cl/2/np.pi))


def plot_binned_spectra_with_error_bars(cl, ell, error, bin_width, label=None, plot_log=True):
    ax = plt.subplot()
    ax.set_xscale('log')
    ax.set_yscale('log')

    if label:
        ax.errorbar(ell, np.sqrt(ell*(ell+1)*cl/2/np.pi), fmt='o', label=label)
    else:
        ax.errorbar(ell, np.sqrt(ell*(ell+1)*cl/2/np.pi), fmt='o')

def annotate_params():
    plt.annotate(r'$\alpha$ $=$ $65$', xy=(15, 6), size=11)
    plt.annotate(r'$\beta$ $=$ $30$', xy=(15, 4), size=11)
    plt.annotate("$T_{precession}$ $=$ $93$ $mins$", xy=(15, 2.6), size=11)
    plt.annotate("$T_{spin}$ $=$ $600s$", xy=(15, 1.6), size=11)
    plt.annotate("$f_{sampling}$ $=$ $10Hz$", xy=(15, 1.1), size=11)
