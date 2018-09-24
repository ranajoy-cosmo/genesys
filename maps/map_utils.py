"""
This module has useful utility functions for manipulating Healpix maps and masks
Healpy mask convention :
    0(False) -> Pixel not seen
    1(True) -> Pixel seen
Masked map mask convention :
    False -> Pixel seen
    True -> Pixel not seen
"""
import numpy as np
import healpy as hp
from genesys.global_config import global_paths
import genesys.noise.noise_utils as nu
import os

#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
#* Masking functions
#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*

def mask_map(sky_map, binary_mask):
    """
    Masks the given sky map with the binary mask.
    At least with healpy version >= 1.11.0, the number of sky maps does not matter.
    For example, passing a single binary mask for a set of I,Q,U sky maps, will apply the same mask for all the three.
    Passing a set of masks in the form of an array will assign the masks to the corresponding sky maps.
    """
    # Checking that the sky map and mask have the same resolution
    assert hp.get_nside(sky_map) == hp.get_nside(binary_mask), "nside of sky map and mask does not match.\nnside of sky map : {}\nnside of mask : {}".format(hp.get_nside(sky_map), hp.get_nside(binary_mask))
        return 0

    sky_map_masked = hp.ma(sky_map) 
    sky_map_masked.mask = np.logical_not(binary_mask)

    return sky_map_masked

def get_mask_from_nan(sky_map):
    """
    Make a mask of the same dimension as the input sky mask.
    A NAN corresponds to an unseen pixel and will be set to 0 in the mask
    """
    mask = np.logical_not(sky_map == np.nan)
    return mask

def get_sky_fraction(binary_mask):
    """
    Returns the fraction of the sky that is valid.
    The function determines the dimension of the binary mask.
    It returns a scalar for a single mask. It returns an array of the same number of masks otherwise.
    """
    dim_mask = hp.maptype(binary_mask)
    n_pix = 12*hp.get_nside(binary_mask)**2
    if dim_mask==0:
        sky_fraction = np.sum(binary_mask.astype(np.float)) / n_pix
    else:
        sky_fraction = np.sum(binary_mask.astype(np.float), axis=1) / n_pix

    return sky_fraction

def ud_grade_mask(mask, nside_out):
    nside_in = hp.get_nside(mask)
    mask_new = hp.ud_grade(mask, nside_out=nside_out)
    if nside_out < nside_in:
        mask_new[mask_new < 1] = 0.0

    return mask_new

def get_equatorial_mask(width, nside, deg=True, pol=True):
    pix = np.arange(12*nside**2)
    theta, phi = hp.pix2ang(nside=nside, ipix=pix)
    if deg:
        width = np.radians(width)

    if pol:
        mask = np.ones((3,12*nside**2))
        mask[...,np.logical_and(theta>np.pi/2-width/2, theta<np.pi/2+width/2)] = 0
    else:
        mask = np.ones(12*nside**2)
        mask[np.logical_and(theta>np.pi/2-width/2, theta<np.pi/2+width/2)] = 0

    return mask

#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*

def deconvolve_map(map_in, fwhm_in=0.0, fwhm_out=0.0, lmax=None, binary_mask=None, pol=False, wiener=True, sky_prior=None):
    
    if fwhm_in == fwhm_out:
        return map_in

    if lmax is None:
        lmax = 3*hp.get_nside(map_in) - 1

    if binary_mask is None:
        binary_mask = get_mask_from_map(map_in)

    f_sky = get_sky_fraction(binary_mask, pol)

    alm_in = su.estimate_alm(map_in, lmax, binary_mask, pol)
    alm_dec = su.deconvolve_alm(alm_in, fwhm_in=fwhm_in, fwhm_out=fwhm_out, f_sky=f_sky, pol=pol, wiener=True, sky_prior=sky_prior)
    map_dec = hp.alm2map(alm_dec, nside=hp.get_nside(map_in), pol=pol)
    
    return map_dec

def fill_empty_pixels(sky_map, max_iter, fail_fill_value=0, pol=True, verbose=False):
    if np.sum(np.isnan(sky_map)) == 0:
        prompt("There are no empty pixels") if verbose
        return

    nside = hp.get_nside(sky_map)

    if pol:
        dim = sky_map.shape[0]
        for i in xrange(max_iter):
            empty_pix = np.where(np.isnan(sky_map[0]))[0]
            theta, phi = hp.pix2ang(nside, empty_pix)
            neighbours = hp.get_all_neighbours(nside, theta, phi).T
            for j in range(dim):
                fill_values = np.nanmean(sky_map[j][neighbours], axis=-1)
                sky_map[j][empty_pix] = fill_values
            if np.sum(np.isnan(sky_map)) == 0:
                break
    else:
        for i in xrange(max_iter):
            empty_pix = np.where(np.isnan(sky_map))[0]
            theta, phi = hp.pix2ang(nside, empty_pix)
            neighbours = hp.get_all_neighbours(nside, theta, phi).T
            fill_values = np.nanmean(sky_map[neighbours], axis=-1)
            sky_map[empty_pix] = fill_values
            if np.sum(np.isnan(sky_map)) == 0:
                break

    num_empty_pix = np.sum(np.isnan(sky_map))
    if num_empty_pix:
        prompt("{} empty pixels remaining after {} iterations. Filling empty pixels with {}\n".format(num_empty_pix, max_iter, fail_fill_value))
        sky_map[np.isnan(sky_map)] = fail_fill_value
