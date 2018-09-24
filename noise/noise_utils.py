"""
This module calculates and converts between the noise parameters for a CMB experiment and also generates the noise power spectra given such parameters.
Calculations are done assuming the observer is at the centre of a 2-Sphere with the last scattering surface being the surface of the sphere and at radius unity away from the observer.
Features include :
    1) Calculate the solid angle of a pixel of a Healpix map from the nside
    2) Time spent by the detector per pixel of the Healpix map for the given nside and time span
    3) Time spent by the detector per arcminute square of the sky for the given time span
    4) Time spent by the detector per unit solid angle of the sky for the given time span
    5) Detector sensitivity to noise rms in each pixel of the Healpix map for the given nside, time span and number of detectors
    6) Detector sensitivity to noise rms in uK.arcmin for the given time span and number of detectors
    7) Detector sensitivity to noise rms in uK.steradian for the given time span and number of detectors
"""

import numpy as np
import healpy as hp
import os
import genesys.numerical.unit_conversion as uc
from genesys.global_config import global_paths

OMEGA_SKY = uc.convert_unit('solid_angle', 1.0, unit_in="full_sky", unit_out="steradian")
OMEGA_DEGREE_SQ = uc.convert_unit('solid_angle', 1.0, unit_in="degree_square", unit_out="steradian")
OMEGA_ARCMIN_SQ = uc.convert_unit('solid_angle', 1.0, unit_in="arcmin_square", unit_out="steradian")

FWHM_FACTOR = 2*np.sqrt(2*np.log(2))                                    #FWHM = FWHM_FACTOR*sigma

fiducial_spectra_folder = os.path.join(global_paths.base_dir, "spectra")
fiducial_spectra_file = os.path.join(global_paths.base_dir, "spectra" , "planck_parameters/r_0.001/lensedtot_cls.npy")              #Fiducial spectra


#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
#* Pixel and time routines
#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*

def pixel_solid_angle(nside, unit_out="arcmin_square", verbose=False):
    """
    Gives the solid angle subtended by a pixel on a Healpix map with a given nside.
    accepted out units : 'arcsec_square', 'arcmin_square', 'degree_square', 'steradian'
    Default units : arcmin_square
    """
    pix_solid_angle_arcmin = hp.nside2pixarea(nside, degrees=False)
    pix_solid_angle = uc.convert_unit('solid_angle', pix_solid_angle_arcmin, unit_in='steradian', unit_out=unit_out)

    if verbose:
        write_string = "nside : {}\n".format(nside)
        write_string += "pixel solid angle : {:.3e} {} ({:.3e} steradian)".format(pix_solid_angle, unit_out, pix_solid_angle_arcmin)
        print(write_string)

    return pix_solid_angle


def time_per_area_unit(t_integration, area_unit_type="arcmin_square", time_unit_in="second", time_unit_out="second", nside=None, verbose=False):
    """
    Time spent by 1 detector per area unit specified.
    The scan is considered to be uniform all over the sky.
    accepted area units : 'pixel', 'arcsec_square', 'arcmin_square', 'degree_square', 'steradian', 'full_sky'
    If area unit is specified as 'pixel', the nside must be specified
    Default input units : seconds
    Default output units : seconds
    """
    if area_unit_type == "pixel":
        assert (nside is not None), "Please enter the nside for unit area type \"pixel\""
        omega_area_unit = {'pixel' : pixel_solid_angle(nside, "steradian")}
    else:
        uc.validate_unit('solid_angle', area_unit_type)
        omega_area_unit = uc.all_unit_dicts['solid_angle']

    time_in_area_unit_tui = t_integration * omega_area_unit[area_unit_type] / OMEGA_SKY
    time_in_area_unit_tuo = uc.convert_unit('time', time_in_area_unit_tui, unit_in=time_unit_in, unit_out=time_unit_out)

    if verbose:
        t_integration_second = uc.convert_unit('time', t_integration, unit_in=time_unit_in, unit_out="second")
        time_in_area_unit_second = uc.convert_unit('time', time_in_area_unit_tui, unit_in=time_unit_in, unit_out="second")
        write_string = "Total integration time : {} {} ({} second)\n".format(t_integration, time_unit_in, t_integration_second)
        write_string += "Area unit : {}\n".format(area_unit_type)
        if area_unit_type == "pixel":
            write_string += "Nside : {}\n".format(nside)
        write_string += "Time in area unit : {:.4e} {} ({:.4e} second)".format(time_in_area_unit_tuo, time_unit_out, time_in_area_unit_second)
        print(write_string)

    return time_in_area_unit_tuo

#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
#* Sensitivity routines and conversions 
#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
area_unit_map = {"arcmin_square" : "arcmin", "degree_square" : "degree", "full_sky" : "sqrt(full_sky)", "steradian" : "radian", "pixel" : "sqrt(pixel)"}

def noise_rms_convert_area_unit(noise_rms_in, area_unit_in="arcmin_square", area_unit_out="arcmin_square", nside_in=None, nside_out=None, verbose=False): 
    """
    Converts the noise rms from uK.unit_in to uK.unit_out.
    This method works on both a scalar quantity as well as on a numpy array.
    For example, a healpix map with each pixel having the rms value can be converted between uK.arcmin to uK rms noise in each pixel.
    """
    time_per_area_unit_in = time_per_area_unit(t_integration=1.0, area_unit_type=area_unit_in, nside=nside_in)
    time_per_area_unit_out = time_per_area_unit(t_integration=1.0, area_unit_type=area_unit_out, nside=nside_out)
    factor = np.sqrt(time_per_area_unit_in / time_per_area_unit_out)
    noise_rms_out = noise_rms_in * factor

    if verbose:
        if np.isscalar(noise_rms_in) == 1:
            write_string = "Noise rms in : {} uK.{}\n".format(noise_rms_in, area_unit_map[area_unit_in])
            write_string += "Noise rms out : {} uK.{}\n".format(noise_rms_out, area_unit_map[area_unit_out])
        else:
            write_string = "Input units : {}\n".format(area_unit_in)
            write_string += "Output units : {}\n".format(area_unit_out)
        if area_unit_in == "pixel":
            write_string += "Nside in : {}\n".format(nside_in)
        if area_unit_out == "pixel":
            write_string += "Nside out : {}\n".format(nside_out)
        print(write_string)

    return noise_rms_out

def detector_sensitivity_to_noise_rms(det_sens, t_integration, area_unit_type="arcmin_square", n_det=1, noise_margin=1.0, detector_yield=1.0, time_unit_in="second", nside=None, verbose=False):
    """
    Converts from the detector sensitivity to the noise in each area unit specified. 
    Assuming an uniform scan over the time of t_integration and by n_det detectors.
    If area unit is specified as 'pixel', the nside must be specified
    det_sens is a numpy array or a single scalar value of unit uK.sqrt(s)
    Array of 3 values returned for each det_sens T, Q, U maps respectively. The noise is sqrt(2) times more for Q and U
    Output unit : uK.sqrt(area_unit_type)
    """

    t_area_unit = time_per_area_unit(t_integration=t_integration, area_unit_type=area_unit_type, nside=nside, time_unit_in=time_unit_in, time_unit_out="second")
    factor = np.array([1.0, np.sqrt(2), np.sqrt(2)])
    noise_rms_in_area_unit_I = det_sens * noise_margin / np.sqrt(t_area_unit) / np.sqrt(n_det*detector_yield)
    if np.isscalar(det_sens):
        noise_rms_in_area_unit = noise_rms_in_area_unit_I * factor
    else:
        noise_rms_in_area_unit = noise_rms_in_area_unit_I[:,None] * factor

    if verbose:
        write_string = "Integration time : {} {}\n".format(t_integration, time_unit_in)
        write_string += "Number of detectors : {}\n".format(n_det)
        write_string += "Detector sensitivity : {} uK.sqrt(s)\n".format(det_sens)
        if area_unit_type == "pixel":
            write_string += "Nside : {}\n".format(nside)
        write_string += "Noise rms : {:.4e} uK.{} in intensity\n".format(noise_rms_in_area_unit[0], area_unit_map[area_unit_type])
        write_string += "          : {:.4e} uK.{} in polarisation".format(noise_rms_in_area_unit[1], area_unit_map[area_unit_type])
        print(write_string)

    return noise_rms_in_area_unit

def noise_rms_to_detector_sensitivity(noise_rms, t_integration, area_unit_type="arcmin_square", nside=None, noise_margin=1.0, detector_yield=1.0, n_det=1, time_unit_in="second", comp="I", verbose=False):
    """
    Converts from the noise in each area unit of a Healpix map to detector sensitivity.
    Assuming an uniform scan over the time of t_integration and by n_det detectors.
    noise_arcmin is a numpy array or a single scalar value of unit uK.area_unit_type.
    det_sens is a numpy array or a scalar value of unit uK.sqrt(s)
    """
    factor = {'I' : 1.0, 'P' : np.sqrt(2)}
    t_area_unit = time_per_area_unit(t_integration=t_integration, area_unit_type=area_unit_type, nside=nside, time_unit_in=time_unit_in, time_unit_out="second")
    det_sens = noise_rms * np.sqrt(t_area_unit) * np.sqrt(n_det*detector_yield) / factor[comp] / noise_margin

    if verbose:
        write_string = "Integration time : {} {}\n".format(t_integration, time_unit_in)
        write_string += "Number of detectors : {}\n".format(n_det)
        write_string += "Noise rms : {} uK.{} in {}\n".format(noise_rms, area_unit_map[area_unit_type], comp)
        if area_unit_type == "pixel":
            write_string += "Nside : {}\n".format(nside)
        write_string += "Detector sensitivity : {} uK.sqrt(s)".format(det_sens)
        print(write_string)

    return det_sens

#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
#* Hitmap, rms maps and sensitivity maps
#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*

def noise_rms_to_uniform_sensitivity_maps(noise_rms, nside, area_unit_type_in="arcmin_square", area_unit_type_out="arcmin_square", nside_in=None, output_type="I+P"): 
    """
    This produces a Healpix sky map of the noise rms.
    The Healpix map will be of a particular nside but the noise rms units is specified explicitly. For example, it can be uK.arcmin but projected on the Healpix pixelised map.
    The output will be a set of 1, 2 or 3 maps for I, (Q,U) or (I,Q,U). The noise rms for the polarisation components will we sqrt(2) times greater than I.
    Currently, this routine is for a single set of (I,Q,U) maps.
    This routine should not be modified for use of multiple bands with different smoothing scales as the noise does not remain white anymore.
    """
    noise_rms_out = noise_rms_convert_area_unit(noise_rms_in=noise_rms, area_unit_in=area_unit_type_in, area_unit_out=area_unit_type_out, nside_in=nside_in, nside_out=nside)
    if output_type == "I+P":
        assert (len(noise_rms) == 3), "Please check the input noise rms values. It does not conform with the output type of Intensity + Polarisation"
        sensitivity_map = np.ones((3,12*nside**2)) * noise_rms_out[:,None]
    elif output_type == "I":
        assert (np.isscalar(noise_rms)), "Please check the input noise rms values. It does not conform with the output type of Intensity alone"
        sensitivity_map = np.ones(12*nside**2) * noise_rms_out
    elif output_type == "P":
        assert (len(noise_rms) == 2), "Please check the input noise rms values. It does not conform with the output type of Polarisation alone"
        sensitivity_map = np.ones((2,12*nside**2)) * noise_rms_out[1]

    return sensitivity_map

def detector_sensitivity_to_uniform_sensitivity_maps(det_sens, nside, t_integration, n_det, area_unit_type="pixel", time_unit_in="second", output_type="I+P"): 
    """
    This produces a Healpix sky map of the noise rms.
    The Healpix map will be of a particular nside but the noise rms units is specified explicitly. For example, it can be uK.arcmin but projected on the Healpix pixelised map.
    The output will be a set of 1, 2 or 3 maps for I, (Q,U) or (I,Q,U). The noise rms for the polarisation components will we sqrt(2) times greater than I.
    Currently, this routine is for a single set of maps.
    This routine should not be modified for use of multiple bands with different smoothing scales as the noise does not remain white anymore.
    """
    noise_rms = detector_sensitivity_to_noise_rms(det_sens=det_sens, t_integration=t_integration, area_unit_type=area_unit_type, n_det=n_det, time_unit_in=time_unit_in, nside=nside)
    if output_type == "I+P":
        sensitivity_map = noise_rms_to_uniform_sensitivity_maps(noise_rms, nside, area_unit_type_in=area_unit_type, area_unit_type_out=area_unit_type, nside_in=nside, output_type=output_type)
    elif output_type == "I":
        sensitivity_map = noise_rms_to_uniform_sensitivity_maps(noise_rms[0], nside, area_unit_type_in=area_unit_type, area_unit_type_out=area_unit_type, nside_in=nside, output_type=output_type)
    elif output_type == "P":
        sensitivity_map = noise_rms_to_uniform_sensitivity_maps(noise_rms[:2], nside, area_unit_type_in=area_unit_type, area_unit_type_out=area_unit_type, nside_in=nside, output_type=output_type)

    return sensitivity_map

def hitmap_to_sensitivity_maps(hitmap, det_sens, sampling_rate, area_unit_type="arcmin_square", spread="real"):
    """
    Produces sensitivity maps given the hitmap, detector sensitivity and sampling rate.
    These are just the square root of the diagonal terms of the covariance matrix.
    The scan is considered to be ideal, that is, all pixels have been observed at sufficient angles such that the off-diagonal terms are vanishing and the ratio of the polarisation sensitivity to that of the intensity tends to 2.
    The noise rms of each individual measurement is given by sensitivity * sampling_rate and has the units of uK.sqrt(s).
    The area_unit_type will determine the units of the sensitivity map.
    If spread == real, the actual hitmap is used.
    If spread == uniform, the hitmap is made uniform and gives us the equivalent sensitivity map for an uniform scan
    """
    nside = hp.get_nside(hitmap)
    npix = 12*nside**2
    if spread == "uniform":
        total_hits = np.sum(hitmap)
        hitmap_local = np.full(npix, total_hits / npix)
    else:
        hitmap_local = hitmap
    area_factor = noise_rms_convert_area_unit(1.0, area_unit_in='pixel', area_unit_out=area_unit_type, nside_in=nside, nside_out=nside)**2
    noise_variance_sample = (det_sens * np.sqrt(sampling_rate))**2
    sensitivity_map_II = np.sqrt( (area_factor * noise_variance_sample / hitmap_local) * np.ones(hitmap.size) ) 

    return np.row_stack((sensitivity_map_II, np.sqrt(2)*sensitivity_map_II, np.sqrt(2)*sensitivity_map_II))

#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
#* Noise power spectra calculations
#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*

def noise_rms_to_noise_mean_Cl(noise_rms_in, lmax, beam_fwhm=0.0, area_unit_type_in="arcmin_square", nside_in=None, f_sky=1.0, input_type="I+P", convolve=False):
    """
    Given the noise rms for a uniform scan of the sky, this gives us the theoretical mean of the noise spectra for TT, EE, BB and TE.
    The input noise rms must be given as either "I", "P", or "I+P".
    Providing the rms for "I+P" means there will be two input values, not necessarily being in the ration of sqrt(2).
    If either, "I" or "P" is provided, the other will be assumed to be in the ratio of sqrt(2).
    Set convolve=True if you want the spectrum to be additionally smoothed by the beam.
    """
    if area_unit_type_in == "pixel":
        assert (nside_in is not None), "Please enter the nside for unit area type \"pixel\""
    if input_type == "I":
        noise_rms_local = np.array([1.0, np.sqrt(2)]) * noise_rms_in
    elif input_type == "P":
        noise_rms_local = np.array([1.0/np.sqrt(2), 1.0]) * noise_rms_in
    else:
        noise_rms_local = noise_rms_in
    w_inv = noise_rms_convert_area_unit(noise_rms_local, area_unit_in=area_unit_type_in, area_unit_out="steradian", nside_in=nside_in)**2
    Bl_squared = hp.gauss_beam(fwhm=np.radians(beam_fwhm/60.0), lmax=lmax, pol=True)**2
    if convolve:
        Bl_squared = 1.0 / Bl_squared
    noise_mean_spectrum = np.empty((4, lmax+1))
    noise_mean_spectrum[0] = np.sqrt(1.0/f_sky) * w_inv[0] / Bl_squared[...,0]
    noise_mean_spectrum[1] = np.sqrt(1.0/f_sky) * w_inv[1] / Bl_squared[...,1]
    noise_mean_spectrum[2] = np.sqrt(1.0/f_sky) * w_inv[1] / Bl_squared[...,2]
    noise_mean_spectrum[3] = np.zeros(lmax+1)
    noise_mean_spectrum[...,:2] = 0.0

    return noise_mean_spectrum

def detector_sensitivity_to_noise_mean_Cl(det_sens, t_integration, n_det, lmax, beam_fwhm=0.0, f_sky=1.0, noise_margin=1.0, detector_yield=1.0, time_unit_in="second", convolve=False):
    """
    Given the sensitivity of the instrument, observation time and number of detectors, and assuming that the scan was uniform, tis gives us the theoretical mean of the noise spectra for TT, EE, and BB.
    """
    noise_rms = detector_sensitivity_to_noise_rms(det_sens=det_sens*noise_margin, t_integration=t_integration, area_unit_type="arcmin_square", n_det=n_det*detector_yield, time_unit_in=time_unit_in)
    noise_mean_spectrum = noise_rms_to_noise_mean_Cl(noise_rms_in=noise_rms[:2], lmax=lmax, beam_fwhm=beam_fwhm, area_unit_type_in="arcmin_square", f_sky=f_sky, input_type="I+P", convolve=convolve)

    return noise_mean_spectrum


def cosmic_variance_Cl(spectra, lmax, f_sky=1.0):
    """
    This gives us the variance on the CMB power spectra due to cosmic variance alone.
    """
    ell = np.arange(lmax + 1)

    cosmic_variance = (2.0 / (2.0*ell + 1) / f_sky) * spectra**2
    cosmic_variance[...,:2] = 0.0

    return cosmic_variance

def noise_rms_to_noise_variance_Cl(noise_rms_in, lmax, beam_fwhm=0.0, area_unit_type_in="arcmin_square", nside_in=None, f_sky=1.0, input_type="I+P", convolve=False):
    """
    Given the noise rms for a uniform scan of the sky, this gives us the theoretical mean of the noise spectra for TT, EE, BB and TE.
    The input noise rms must be given as either "I", "P", or "I+P".
    Providing the rms for "I+P" means there will be two input values, not necessarily being in the ration of sqrt(2).
    If either, "I" or "P" is provided, the other will be assumed to be in the ratio of sqrt(2).
    Set convolve=True if you want the spectrum to be additionally smoothed by the beam.
    """
    if area_unit_type_in == "pixel":
        assert (nside_in is not None), "Please enter the nside for unit area type \"pixel\""
    if input_type == "I":
        noise_rms_local = np.array([1.0, np.sqrt(2)]) * noise_rms_in
    elif input_type == "P":
        noise_rms_local = np.array([1.0/np.sqrt(2), 1.0]) * noise_rms_in
    else:
        noise_rms_local = noise_rms_in
    w_inv = noise_rms_convert_area_unit(noise_rms_local, area_unit_in=area_unit_type_in, area_unit_out="steradian", nside_in=nside_in)**2
    Bl_squared = hp.gauss_beam(fwhm=np.radians(beam_fwhm/60.0), lmax=lmax, pol=True)**2
    if convolve:
        Bl_squared = 1.0 / Bl_squared
    ell = np.arange(lmax+1)
    noise_variance_spectrum = np.empty((4, lmax+1))
    noise_variance_spectrum[0] = (2.0 / (2.0*ell + 1) / f_sky) * (w_inv[0] / Bl_squared[...,0])**2
    noise_variance_spectrum[1] = (2.0 / (2.0*ell + 1) / f_sky) * (w_inv[1] / Bl_squared[...,1])**2
    noise_variance_spectrum[2] = (2.0 / (2.0*ell + 1) / f_sky) * (w_inv[1] / Bl_squared[...,2])**2
    noise_variance_spectrum[3] = np.zeros(lmax+1)
    noise_variance_spectrum[...,:2] = 0.0

    return noise_variance_spectrum

def detector_sensitivity_to_noise_variance_Cl(det_sens, t_integration, n_det, lmax, beam_fwhm=0.0, f_sky=1.0, noise_margin=1.0, detector_yield=1.0, time_unit_in="second", convolve=False):
    """
    Given the sensitivity of the instrument, observation time and number of detectors, and assuming that the scan was uniform, tis gives us the theoretical mean of the noise spectra for TT, EE, BB and TE.
    """
    noise_rms = detector_sensitivity_to_noise_rms(det_sens=det_sens*noise_margin, t_integration=t_integration, area_unit_type="arcmin_square", n_det=n_det*detector_yield, time_unit_in=time_unit_in)
    noise_variance_spectrum = noise_rms_to_noise_variance_Cl(noise_rms_in=noise_rms[:2], lmax=lmax, beam_fwhm=beam_fwhm, area_unit_type_in="arcmin_square", f_sky=f_sky, input_type="I+P", convolve=convolve)

    return noise_variance_spectrum

def noise_rms_to_variance_Cl(noise_rms_in, spectra, lmax, beam_fwhm=0.0, area_unit_type_in="arcmin_square", nside_in=None, f_sky=1.0, input_type="I+P"):
    noise_variance_spectrum = noise_rms_to_noise_variance_Cl(noise_rms_in, lmax, beam_fwhm=beam_fwhm, area_unit_type_in=area_unit_type, nside_in=nside_in, f_sky=f_sky, input_type=input_type)
    cosmic_variance = cosmic_variance_Cl(spectra=spectra, lmax=lmax, f_sky=f_sky)

    return noise_variance_spectrum + cosmic_variance

def detector_sensitivity_to_variance_Cl(det_sens, spectra, t_integration, n_det, lmax, beam_fwhm=0.0, f_sky=1.0, time_unit_in="second"):
    noise_variance_spectrum = detector_sensitivity_to_noise_variance_Cl(det_sens=det_sens, t_integration=t_integration, n_det=n_det, lmax=lmax, beam_fwhm=beam_fwhm, f_sky=f_sky, time_unit_in=time_unit_in)
    cosmic_variance = cosmic_variance_Cl(spectra=spectra, lmax=lmax, f_sky=f_sky)

    return noise_variance_spectrum + cosmic_variance

#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
#* Random noise maps
#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*

def noise_rms_to_uniform_random_noise_maps(noise_rms, nside, lmax=None, beam_fwhm=0.0, area_unit_type_in="arcmin_square", nside_in=None, input_type="I+P", convolve=False, seed=None):
    if lmax == None:
        lmax = np.int(2.65*nside)

    spectra_noise_mean = noise_rms_to_noise_mean_Cl(noise_rms_in=noise_rms, lmax=lmax, beam_fwhm=beam_fwhm, f_sky=1.0, area_unit_type_in=area_unit_type_in, nside_in=nside_in, convolve=convolve, input_type=input_type) 
    np.random.seed(seed)
    random_noise_map = hp.synfast(cls=spectra_noise_mean, nside=nside, lmax=lmax, pol=True, new=True)

    return random_noise_map

def detector_sensitivity_to_uniform_random_noise_maps(det_sens, t_integration, n_det, nside, lmax=None, beam_fwhm=0.0, time_unit_in="second", convolve=False, seed=None):
    if lmax == None:
        lmax = np.int(2.65*nside)

    spectra_noise_mean = detector_sensitivity_to_noise_mean_Cl(det_sens=det_sens, t_integration=t_integration, n_det=n_det, lmax=lmax, beam_fwhm=beam_fwhm, f_sky=1.0, time_unit_in=time_unit_in, convolve=convolve) 
    np.random.seed(seed)
    random_noise_map = hp.synfast(cls=spectra_noise_mean, nside=nside, lmax=lmax, pol=True, new=True)

    return random_noise_map
