#!/usr/bin/env python

import numpy as np
import healpy as hp
import os

def make_map_from_spectra(spectra, nside, lmax, resol, seed=1234, pixwin=False, verbose=False):
    np.random.seed(seed)
    cmb_map = hp.synfast(spectra, nside=nside, lmax=lmax, fwhm=np.radians(resol/60.0), new=True, pol=True, pixwin=pixwin, verbose=verbose)

    return cmb_map

if __name__=="__main__":
    r_list = [0.0, 0.001, 0.01, 0.02, 0.03, 0.04, 0.05, 0.1] 
    spectra_dir = "../spectra/spectra_files"
    spectra_tag = "planck_params_r_variation_formatted"

    out_dir = os.path.join("map_files/CMB_Fiducial",spectra_tag)
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)
    nside_resol_lmax_list = ((16, 600.0, 47), (128, 60.0, 383))

    for r in r_list:
        print("Doing r = {}".format(r))
        spectra_file = os.path.join(spectra_dir, spectra_tag, 'r'+str(r), 'lensedtotCls.npy')
        if os.path.isfile(spectra_file):
            print("{} : EXISTS".format(spectra_file))
        else:
            print("{} : DOES NOT EXIST".format(spectra_file))
        spectra = np.load(spectra_file)
        print(spectra.shape)
        for nside, resol, lmax in nside_resol_lmax_list:
            print("nside : {}\nresol : {}\nlmax : {}\n".format(nside, resol, lmax))
            cmb_map = make_map_from_spectra(spectra, nside, lmax, resol)
            out_dir_r = os.path.join(out_dir, 'r'+str(r))
            if not os.path.exists(out_dir_r):
                os.makedirs(out_dir_r)
            hp.write_map(os.path.join(out_dir_r, 'sky_map_ns_'+str(nside)+'_'+str(int(resol))+'_arcmins.fits'), cmb_map)
