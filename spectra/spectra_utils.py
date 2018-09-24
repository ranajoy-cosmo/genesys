"""
This module contains a number of tools for the estimation and manupulation of power spectra
"""
#!/usr/bin/env python

import numpy as np
import healpy as hp
import genesys.numerical.filters as fl
import genesys.maps.map_utils as mu
from genesys.global_config import global_paths
import os

#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
# Power spectrum estimation routines
#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*

def estimate_cl(sky_map, lmax, binary_mask=None, beam_fwhm=0.0, method="anafast"):
    """
    Estimates the TT, EE, BB and TE auto/cross power spectrum.
    The sky map(s) will be masked by the binary mask provided.
    beam_fwhm input is in arcmins. It is converted internally to radians.
    """
    if binary_mask is None:
        binary_mask = mu.get_mask_from_map(sky_map)
    sky_map_masked = mu.mask_map(sky_map, binary_mask=binary_mask, pol=pol)
    f_sky = mu.get_sky_fraction(binary_mask, pol=pol)

    Bl = hp.gauss_beam(fwhm=np.radians(fwhm/60.0), lmax=lmax, pol=pol)

    if pol:
        spectra = hp.anafast(sky_map_masked.filled(), lmax=lmax)[:4]
        # spectra /= f_sky[...,None]*Bl.T**2
        spectra /= f_sky[1]*Bl.T**2
        # spectra /= Bl.T**2
    else:
        spectra = hp.anafast(sky_map_masked.filled(), lmax=lmax)
        spectra /= f_sky*Bl**2

    return spectra


def estimate_alm(sky_map, lmax, binary_mask=None, pol=False):
    if binary_mask == None:
        binary_mask = mu.get_mask_from_map(sky_map)
    sky_map_masked = mu.mask_map(sky_map, binary_mask=binary_mask, pol=pol)

    alm = hp.map2alm(sky_map_masked, lmax=lmax)

    return alm


def deconvolve_spectra(spectra, fwhm, lmax, pol=True):
    Bl = hp.gauss_beam(fwhm=fwhm, lmax=lmax, pol=pol)
    if pol:
        spectra /= Bl.T**2
    else:
        spectra /= Bl**2

    return spectra


def wiener_filter_for_alm(alm, lmax=None, fwhm=0.0, f_sky=1.0, sky_prior=None):
    if lmax is None:
        lmax = hp.Alm.getlmax(len(alm), None)

    if sky_prior is None:
        spectra_th = np.load("spectra_files/planck_params_r_variation_formatted/r0.0/lensedtotCls.npy")[0,:lmax+1]
    else:
        spectra_th = estimate_cl(sky_prior, lmax, fwhm=fwhm, pol=False)

    Bl = hp.gauss_beam(fwhm=fwhm, lmax=lmax, pol=False)
    spectra_ob = hp.alm2cl(alm, lmax_out=lmax)
    spectra_ob /= f_sky*Bl**2

    filter_response = spectra_ob/spectra_th
    filter_response[:2] = 1.0
    filter_response = np.sqrt(filter_response)

    return filter_response


def deconvolve_alm(alms, lmax=None, fwhm_in=0.0, fwhm_out=0.0, f_sky=1.0, pol=False, wiener=True, sky_prior=None):

    fwhm = np.sqrt(fwhm_in**2 - fwhm_out**2)
    if fwhm == 0.0:
        return alms
    factor = (2.0*np.sqrt(2.0*np.log(2.0)))
    sigma = fwhm/factor
    
    retalm = []

    if lmax is None:
        lmax = hp.Alm.getlmax(len(alms[0] if pol else alms), None)
    
    if wiener:
        wiener_filter = wiener_filter_for_alm(alms[0] if pol else alms, lmax, f_sky=f_sky, sky_prior=sky_prior, fwhm=fwhm_in)
        wiener_smooth = fl.filter_butter(wiener_filter, lmax, 100)
        #wiener_smooth[:1500] = 1.0

    if pol:
        for ialm, alm in enumerate(alms):
            ell = np.arange(lmax + 1)
            if ialm >= 1:
                s = 2
            else:
                s = 0
            fact = np.exp(0.5*(ell*(ell + 1) - s**2)*sigma**2)
            if wiener:
                fact /= wiener_smooth
            res = hp.almxfl(alm, fact, inplace=False)
            retalm.append(res)
    else:
        lmax = hp.Alm.getlmax(len(alms), None)
        ell = np.arange(lmax + 1)
        fact = np.exp(0.5*(ell*(ell + 1))*sigma**2)
        if wiener:
            fact /= wiener_smooth
        retalm = hp.almxfl(alms, fact, inplace=False)

    return retalm 

