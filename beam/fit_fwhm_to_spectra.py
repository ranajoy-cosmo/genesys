#!/usr/bin/env python

import numpy as np
import healpy as hp
from lmfit import minimize, Parameters, fit_report
import simulation.power_spectrum.noise_power_spectra as nps
import simulation.power_spectrum.spectra_tools as st
from simulation.global_config import global_paths
import sys
import os

fid_spectra_folder = os.path.join(global_paths.base_dir, "spectra")
fid_spectra_file = os.path.join(fid_spectra_folder, "r_0001/lensedtot_cls.npy")              #Fiducial spectra

def residual(params, ell, out_spectra, var_spectra, fid_spectra):
    factor = 2*np.sqrt(2*np.log(2))
    beam_sigma = np.radians(params['fwhm']/60.0)/factor
    model = fid_spectra*np.exp(-ell*(ell+1)*beam_sigma**2)

    return (out_spectra - model)/np.sqrt(var_spectra)

def estimate_fwhm(spectra, lmax):

    ell = np.arange(lmax+1)[2:]
    spectra_var = nps.cosmic_variance_Cl(lmax)[0,2:]
    spectra_fid = np.load(fid_spectra_file)[0,2:lmax+1]

    params = Parameters()
    params.add('fwhm', value=5.0)

    out = minimize(residual, params, args=(ell, spectra, spectra_var, spectra_fid))

    return out


if __name__=="__main__":
    generate_new = sys.argv[1]

    if generate_new == "new":
        map_file = sys.argv[2]
        lmax = int(sys.argv[3])
        sky_map = hp.read_map(map_file)
        spectra = st.estimate_cl(sky_map, lmax=lmax, pol=False)[2:lmax+1]
        np.save("temp_spectra", spectra)
    else:
        spectra_file = sys.argv[2]
        lmax = int(sys.argv[3])
        spectra = np.load(spectra_file)[2:lmax+1]

    out = estimate_fwhm(spectra, lmax)
    print(fit_report(out))
