#!/usr/bin/env python

import numpy as np

def bin_spectra_old(spectra, bins, pol=False):
    n_bins = len(bins) - 1
    ell = np.arange(bins[-1] + 1)
    ell_bins = np.empty(n_bins)
    if pol:
        binned_spectra = np.empty((3, n_bins))
    else:
        binned_spectra = np.empty(n_bins)

    for i in range(n_bins):
        bin_low = bins[i]
        bin_high = bins[i+1]
        if pol:
            for j in range(3):
                binned_spectra[j, i] = np.mean(spectra[j, bin_low:bin_high])
        else:
            binned_spectra[i] = np.mean(spectra[bin_low:bin_high])
        ell_bins[i] = np.mean(ell[bin_low:bin_high])

    return ell_bins, binned_spectra

def bin_spectra(spectra, lmax, pol=False):
    if pol:
        binned_spectra = np.empty((3, lmax+1))
    else:
        binned_spectra = np.empty(lmax+1)

    for ell in range(lmax+1):
        ell_low = int(ell - 0.1*ell)
        ell_high = int(ell + 0.1*ell)
        if pol:
            for j in range(3):
                binned_spectra[j, ell] = np.mean(spectra[j, ell_low:ell_high])
        else:
            binned_spectra[ell] = np.mean(spectra[ell_low:ell_high])

    return binned_spectra

def bin_error_bar(error, bins, pol=False):
    n_bins = len(bins) - 1
    ell = np.arange(bins[-1] + 1)
    bin_width = np.empty(n_bins)
    if pol:
        binned_error = np.empty((3, n_bins))
    else:
        binned_error = np.empty(n_bins)

    for i in range(n_bins - 1):
        bin_low = bins[i]
