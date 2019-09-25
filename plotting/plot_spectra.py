#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import matplotlib
import os
import sys
from genesys.global_config import global_paths

pi = np.pi

class Spectra_Frame:
    """
    A wrapper class for matplotlib.pyplot for plotting power spectra
    """
    def __init__(self, **kwargs):
        """
        Default constructor that creates the figure and subplot axes objects.
        The user may pass standars matplotlib keyworded arguments such as
        Keyworded arguments passed to pyplot.figure():
            num, figsize, dpi, facecolor, edgecolor,frameon,...
        nrows, ncol, figsize.....
        """
        self.fig, self.axes = plt.subplots(**kwargs, squeeze=False)

    def add_title_to_plot(self, title, fontsize=12):
        """
        Add a title to the main figure
        """
        self.fig.suptitle(title, fontsize=fontsize)

    def add_title_to_subplot(self, subplot_number, title, fontsize):
        """
        Add a title to the main figure
        """
        self.axes[subplot_number].set_title(title, fontsize=fontsize)

    def plot_spectrum(self, spectra_obj, columns, ells, lmax, subplot_num, label, annotate_label, **kwargs):


    def plot_r_variation_fid_spectra(self, subplot_number=(0,0), fid_spectra_dir=None, default_file_name=None, plot_list=None, r_list=None, lensed=True, annotate_size=12, figsize=(20,10), axis_size=12, plot_x_log=True, plot_y_log=True, **kwargs):
        """
        Plot the fiducial spectra with different r_values
        The fiducial spectra are in fid_spectra_dir and in default named:
            unlensed
        Default parameter values if passed as None
            fid_spectra_dir: 
            default_file_name:
            plot_list: ['EE','BB']
            r_list: [0.0, 0.001, 0.01]
        """
        if fid_spectra_dir == None:
            fid_spectra_dir = some_dir
        if default_file_name == None:
            default_file_name = first_rvalue_last
        if plot_list == None:
            plot_list = ['EE','BB']
        if r_list == None:
            r_list = [0.0, 0.001, 0.01]

        full_spectra_dir = os.path.join(global_paths.spectra_dir, fid_spectra_dir)
        r_spectra_dirs_dict = {'0.0': os.path.join(full_spectra_dir, 'r0.0', '0.001': os.path.join(full_spectra_dir, 'r0.001', '0.01': os.path.join(full_spectra_dir, 'r0.01'} 

        lensed_total_name = "lensedtot_Cls.npy"
        total_name = "tot_Cls.npy"

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
            spectra = np.load(os.path.join(r0_0_spectra_dir, lensed_total_name))[...,2:lmax+1]
        else:
            spectra = np.load(os.path.join(r0_0_spectra_dir, total_name))[...,2:lmax+1]

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
            spectra = np.load(os.path.join(r0_0_spectra_dir, lensed_total_name))[..., 2:lmax+1]
            plotter(ell, ell*(ell+1)*spectra[2]/2/pi, 'k', **kwargs)

            spectra = np.load(os.path.join(r0_0_spectra_dir, total_name))[..., 2:lmax+1]
            plotter(ell, ell*(ell+1)*spectra[2]/2/pi, 'k--', **kwargs)

            spectra = np.load(os.path.join(r0_01_spectra_dir, lensed_total_name))[..., 2:lmax+1]
            plotter(ell, ell*(ell+1)*spectra[2]/2/pi, 'k', **kwargs)

            spectra = np.load(os.path.join(r0_01_spectra_dir, total_name))[..., 2:lmax+1]
            plotter(ell, ell*(ell+1)*spectra[2]/2/pi, 'k--', **kwargs)

            spectra = np.load(os.path.join(r0_001_spectra_dir, lensed_total_name))[..., 2:lmax+1]
            plotter(ell, ell*(ell+1)*spectra[2]/2/pi, 'k', **kwargs)

            spectra = np.load(os.path.join(r0_001_spectra_dir, total_name))[..., 2:lmax+1]
            plotter(ell, ell*(ell+1)*spectra[2]/2/pi, 'k--', **kwargs)

            plt.annotate("BB", xy=(1.25, 0.00000166), size=annotate_size)

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

def plot_Cl_binned(ell, cl, label=None):
    if label:
        plt.plot(ell, ell*(ell+1)*cl/2/np.pi, label=label)
    else:
        plt.plot(ell, ell*(ell+1)*cl/2/np.pi)

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
