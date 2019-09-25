#!/usr/bin/env python

import sys
import os
import numpy as np
import camb

if __name__=="__main__":
    # It is assumed that the camb parameter file is in the camb_params folder
    param_file_name = sys.argv[1]
    # It is assumed the output folder will be a child directory of spectra_files
    output_folder_name = sys.argv[2]
    # Prefix folder will be a sub-folder of output_folder
    prefix_folder_name = sys.argv[3]
    # Output prefix is prepended to each output file name
    output_prefix = sys.argv[4]

    output_folder = os.path.join(os.getcwd(), "spectra_files", output_folder_name)
    if not os.path.exists(output_folder):
        os.makedirs(output_folder)
    if prefix_folder_name != None:
        output_folder = os.path.join(output_folder, prefix_folder_name)
        if not os.path.exists(output_folder):
            os.makedirs(output_folder)

    param_file_full_name = os.path.join(os.getcwd(), "camb_params", param_file_name)
    param = camb.read_ini(param_file_full_name)

    camb_results = camb.get_results(param)
    Dl_power = camb_results.get_cmb_power_spectra(param, CMB_unit='muK')

    for name in list(Dl_power.keys()):
        ell = np.arange(Dl_power[name].shape[0])[2:]
        Cl_power = 2*np.pi * Dl_power[name][2:] / ell[...,None] / (ell[...,None]+1)
        if output_prefix == "None":
            np.save(os.path.join(output_folder, name), Cl_power)
        else:
            np.save(os.path.join(output_folder, output_prefix + '_' + name), Cl_power)
