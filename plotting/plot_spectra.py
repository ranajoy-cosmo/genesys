#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import os
import sys
from genesys.global_config import global_paths

pi = np.pi

def plot_fiducial(lmax, plot_list=["TT", "EE", "BB", "TE"], lensed=True, fid_spectra_dir="pico_spectra", annotate_size=12, figsize=(20,10), axis_size=12, plot_x_log=True, plot_y_log=True, **kwargs):

    full_spectra_dir = os.path.join(global_paths.spectra_dir, fid_spectra_dir)

    if plot_x_log==True and plot_y_log==True:
        plotter = plt.loglog
    elif plot_x_log==True:
        plotter = plt.semilogx
    elif plot_y_log==True:
        plotter = plt.semilogy
    else:
        plotter = plt.plot

    matplotlib.rc('xtick', labelsize=axis_size)
    matplotlib.rc('ytick', labelsize=axis_size)

    fig = plt.figure(figsize=figsize)

    ax = fig.add_subplot(111)
    
    #Plotting the TT and EE spectra

    if lensed:
        spectra = np.load(os.path.join(full_spectra_dir, "r0.0", "lensedtot_Cls.npy"))[...,2:lmax+1]
    else:
        spectra = np.load(os.path.join(full_spectra_dir, "r0.0", "tot_Cls.npy"))[...,2:lmax+1]

    if lmax > spectra.shape[1] + 1:
        print("Requested lmax of {} greater than lmax of provided spectra. Changing lmax to {}.".format(lmax, spectra.shape[1] - 1))
        lmax = spectra.shape[1] + 1

    ell = np.arange(2,lmax+1)

    if "TT" in plot_list:
        plotter(ell, ell*(ell+1)*spectra[0]/2/pi, 'k', **kwargs)
        plt.annotate("TT", xy=(1.25, 960), size=annotate_size)

    if "EE" in plot_list:
        plotter(ell, ell*(ell+1)*spectra[1]/2/pi, 'k', **kwargs)
        plt.annotate("EE", xy=(1.25, 0.04), size=annotate_size)

    if "TE" in plot_list:
        plotter(ell, ell*(ell+1)*spectra[3]/2/pi, 'k', **kwargs)
        plt.annotate("EE", xy=(1.25, 0.04), size=annotate_size)

    #Plotting the BB fiducial spectra

    if "BB" in plot_list:
        spectra = np.load(os.path.join(full_spectra_dir, "r0.0", "lensedtot_Cls.npy"))[..., 2:lmax+1]
        plotter(ell, ell*(ell+1)*spectra[2]/2/pi, 'k', **kwargs)

        spectra = np.load(os.path.join(full_spectra_dir, "r0.001", "tot_Cls.npy"))[..., 2:lmax+1]
        plotter(ell, ell*(ell+1)*spectra[2]/2/pi, 'k--', **kwargs)

        spectra = np.load(os.path.join(full_spectra_dir, "r0.001", "lensedtot_Cls.npy"))[..., 2:lmax+1]
        plotter(ell, ell*(ell+1)*spectra[2]/2/pi, 'k', **kwargs)

        spectra = np.load(os.path.join(full_spectra_dir, "r0.01", "tot_Cls.npy"))[..., 2:lmax+1]
        plotter(ell, ell*(ell+1)*spectra[2]/2/pi, 'k--', **kwargs)

        spectra = np.load(os.path.join(full_spectra_dir, "r0.01", "lensedtot_Cls.npy"))[..., 2:lmax+1]
        plotter(ell, ell*(ell+1)*spectra[2]/2/pi, 'k', **kwargs)
#
        #  spectra = np.load(os.path.join(full_spectra_dir, "r0.1", "totCls.npy"))[..., 2:lmax+1]
        #  plotter(ell, ell*(ell+1)*spectra[2]/2/pi, 'k--', **kwargs)
#
        #  spectra = np.load(os.path.join(full_spectra_dir, "r0.1", "lensedtotCls.npy"))[..., 2:lmax+1]
        #  plotter(ell, ell*(ell+1)*spectra[2]/2/pi, 'k', **kwargs)

        plt.annotate("BB", xy=(1.25, 0.00000166), size=annotate_size)

        #  plt.annotate("r = 0.1", xy=(1.25, 0.00166), size=annotate_size-2)

        plt.annotate("r = 0.01", xy=(1.25, 0.000166), size=annotate_size-2)

        plt.annotate("r = 0.001", xy=(1.25, 0.0000166), size=annotate_size-2)


    plt.show()

    return fig, ax

def make_decorations(ax, ylim=[2e-7,100], xlim=[1.0,2000], leg_loc=None, label_fontsize=12, legend_fontsize=12, axis_fontsize=12):
    plt.ylim(ylim)
    plt.xlim(xlim)
    
    plt.xlabel('$\ell$', fontsize=label_fontsize)
    plt.ylabel('$\ell(\ell+1)C_{\ell}/2\pi$ $[\mu K^2]$', fontsize=label_fontsize)
    
    if leg_loc == None:
        leg_loc = "upper left"
    plt.legend(loc=leg_loc , prop={'size':legend_fontsize})

    for ticklabel in ax.yaxis.get_ticklabels():
        ticklabel.set_rotation("vertical")

    for item in ([ax.yaxis.label] + [ax.xaxis.label] + ax.get_xticklabels() + ax.get_yticklabels()):
        item.set_fontsize(axis_fontsize)

def plot_Cl(cl, lmax=None, label=None, **kwargs):
    if lmax==None:
        lmax = cl.size + 1
    ell = np.arange(lmax+1)[2:]
    if label:
        plt.plot(ell, ell*(ell+1)*cl[2:lmax+1]/2/pi, label=label, **kwargs)
    else:
        plt.plot(ell, ell*(ell+1)*cl[2:lmax+1]/2/pi, **kwargs)

def plot_Cl_error_bars(cl, sigma_cl, lmax=None, color='r', alpha=0.5):
    if lmax==None:
        lmax = cl.size - 1
    ell = np.arange(lmax+1)
    upper = np.sqrt(ell*(ell+1)*(cl + sigma_cl)[:lmax+1]/2/pi)
    lower = np.sqrt(ell*(ell+1)*(cl - sigma_cl)[:lmax+1]/2/pi)
    plt.fill_between(ell, upper, lower, facecolor=color, alpha=alpha, lw=0)

def plot_Cl_binned(ell, cl, label=None, plot_log=True):
    if plot_log:
        plotter = plt.loglog
    else:
        plotter = plt.plot

    if label:
        plotter(ell, ell*(ell+1)*cl/2/np.pi, label=label)
    else:
        plotter(ell, ell*(ell+1)*cl/2/np.pi)

def plot_Cl_binned_with_error_bars(ell, Cl, std_Cl, bin_width, **kwargs):
    plt.errorbar(ell, ell*(ell+1)*Cl/2/pi, fmt='o', yerr=ell*(ell+1)*std_Cl/2/pi, xerr=bin_width, **kwargs)

def plot_Dl_binned_with_error_bars(ell, Dl, std_Dl, bin_width, **kwargs):
    plt.errorbar(ell, Dl, fmt='o', yerr=std_Dl, xerr=bin_width, **kwargs)

def annotate_params():
    plt.annotate(r'$\alpha$ $=$ $65$', xy=(15, 6), size=11)
    plt.annotate(r'$\beta$ $=$ $30$', xy=(15, 4), size=11)
    plt.annotate("$T_{precession}$ $=$ $93$ $mins$", xy=(15, 2.6), size=11)
    plt.annotate("$T_{spin}$ $=$ $600s$", xy=(15, 1.6), size=11)
    plt.annotate("$f_{sampling}$ $=$ $10Hz$", xy=(15, 1.1), size=11)
