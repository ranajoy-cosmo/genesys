#!/usr/bin/env python 

import numpy as np
import os
import sys

"""
The output from CAMB are really the Dl and not the CL.
This means, the output is of the form ell*(ell+1)*Cl / 2\pi

The first index of the camb output spectra arrays are of the form:
Scalar spectra: The scalar contribution
_scalCls: [ell,TT,EE,TE,\phi\phi,T\phi]
Tensor spectra: The tensor contribution
_tensCls: [ell,TT,EE,BB,TE]
Lensed scalar spectra: The lensed scalar spectra
_lensedCls: [ell,TT,EE]
Total spectra: Scalar + Tensor
_totCls: [ell,TT,EE,BB,TE]
Lensed total spectra: Lensed spectra of Scalar + Tensor
_lensedtotCls: [ell,TT,EE,BB,TE]

The multipoles range from 2 to lmax
"""

def format_spectra(spectra, lmax):
    ell = np.arange(2,lmax+1) 
    norm_factor = ell*(ell + 1)/(2*np.pi)
    spectra_formatted = np.zeros((4, lmax+1))
    spectra_formatted[...,2:] = spectra[:lmax-1,1:5].T / norm_factor
    return spectra_formatted

if __name__=="__main__":
    lmax = 5000 
    spectra_tag = "pico_spectra"
    camb_output_folder = os.path.join("camb_output", spectra_tag)
    formatted_output_folder = os.path.join("spectra_files", spectra_tag)
    
    #  tau_list = [0.055]
    r_list = [0.0, 0.0001, 0.001, 0.003, 0.01]

    #  for tau in tau_list:
    for r in r_list:
        #  spectra_folder = os.path.join(camb_output_folder, 'tau'+str(tau)+'_r'+str(r))
        #  unlensed_in_file = os.path.join(spectra_folder, 'tau'+str(tau)+'_r'+str(r)+'_totCls.dat')
        #  lensed_in_file = os.path.join(spectra_folder, 'tau'+str(tau)+'_r'+str(r)+'_lensedtotCls.dat')
        spectra_folder = os.path.join(camb_output_folder, 'r'+str(r))
        scal_in_file = os.path.join(spectra_folder, 'spectra_scalCls.dat')
        lensed_in_file = os.path.join(spectra_folder, 'spectra_lensedCls.dat')
        tens_in_file = os.path.join(spectra_folder, 'spectra_tensCls.dat')
        tot_in_file = os.path.join(spectra_folder, 'spectra_totCls.dat')
        lensedtot_in_file = os.path.join(spectra_folder, 'spectra_lensedtotCls.dat')
        #  if os.path.isfile(unlensed_in_file):
            #  print("{} : EXISTS".format(unlensed_in_file))
        #  else:
            #  print("{} : DOES NOT EXIST".format(unlensed_in_file))

        read_spectra_scal = np.loadtxt(scal_in_file)
        read_spectra_lensed = np.loadtxt(lensed_in_file)
        read_spectra_tens = np.loadtxt(tens_in_file)
        read_spectra_tot = np.loadtxt(tot_in_file)
        read_spectra_lensedtot = np.loadtxt(lensedtot_in_file)

        spectra_scal = format_spectra(read_spectra_scal, lmax)
        spectra_lensed = format_spectra(read_spectra_lensed, lmax)
        spectra_tens = format_spectra(read_spectra_tens, lmax)
        spectra_tot = format_spectra(read_spectra_tot, lmax)
        spectra_lensedtot = format_spectra(read_spectra_lensedtot, lmax)

        if not os.path.exists(formatted_output_folder):
            os.makedirs(formatted_output_folder)
        #  spectra_output_folder = os.path.join(formatted_output_folder, 'tau'+str(tau)+'_r'+str(r))
        spectra_output_folder = os.path.join(formatted_output_folder, 'r'+str(r))
        if not os.path.exists(spectra_output_folder):
            os.makedirs(spectra_output_folder)
        scal_out_file = os.path.join(spectra_output_folder, 'scal_Cls.npy')
        lensed_out_file = os.path.join(spectra_output_folder, 'lensed_Cls.npy')
        tens_out_file = os.path.join(spectra_output_folder, 'tens_Cls.npy')
        tot_out_file = os.path.join(spectra_output_folder, 'tot_Cls.npy')
        lensedtot_out_file = os.path.join(spectra_output_folder, 'lensedtot_Cls.npy')

        np.save(scal_out_file, spectra_scal)
        np.save(lensed_out_file, spectra_lensed)
        np.save(tens_out_file, spectra_tens)
        np.save(tot_out_file, spectra_tot)
        np.save(lensedtot_out_file, spectra_lensedtot)
