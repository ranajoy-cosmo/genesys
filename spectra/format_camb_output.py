#!/usr/bin/env python 

import numpy as np
import os
import sys

def format_spectra(spectra, lmax):
    ell = np.arange(2,lmax+1) 
    norm_factor = ell*(ell + 1)/(2*np.pi)
    spectra_formatted = np.zeros((4, lmax+1))
    spectra_formatted[...,2:] = spectra/norm_factor
    return spectra_formatted

if __name__=="__main__":
    r_list = [0.0, 0.001, 0.01]
    lmax = int(sys.argv[1])
    spectra_tag = sys.argv[2]
    spectra_folder = os.path.join("spectra_files", spectra_tag)

    for r in r_list:

        lensed_in_file = os.path.join(spectra_folder, 'tau0.055_r'+str(r)+'_lensedtotCls.dat')
        unlensed_in_file = os.path.join(spectra_folder, 'tau0.055_r'+str(r)+'_totCls.dat')
        # if os.path.isfile(unlensed_in_file):
            # print("{} : EXISTS".format(unlensed_in_file))
        # else:
            # print("{} : DOES NOT EXIST".format(unlensed_in_file))

        read_spectra_lensed = np.loadtxt(lensed_in_file)[:lmax-1,1:5].T
        read_spectra_unlensed = np.loadtxt(unlensed_in_file)[:lmax-1,1:5].T

        spectra_lensed = format_spectra(read_spectra_lensed, lmax) 
        spectra_unlensed = format_spectra(read_spectra_unlensed, lmax) 


        spectra_out_folder = spectra_folder+"_formatted"
        if not os.path.exists(spectra_out_folder):
            os.makedirs(spectra_out_folder)
        r_dir = os.path.join(spectra_out_folder, 'r'+str(r))
        if not os.path.exists(r_dir):
            os.makedirs(r_dir)
        lensed_out_file = os.path.join(r_dir, 'lensedtotCls.npy')
        unlensed_out_file = os.path.join(r_dir, 'totCls.npy')
        np.save(lensed_out_file, spectra_lensed)
        np.save(unlensed_out_file, spectra_unlensed)
