#!/usr/bin/env python

import os
import sys
import argparse
from genesys import global_paths

rel_output_dir = os.path.join(global_paths['base_dir'], 'output')
rel_data_dir = os.path.join(global_paths['base_dir'], 'data')
rel_maps_dir = os.path.join(global_paths['base_dir'], 'maps', 'map_files')
rel_spectra_dir = os.path.join(global_paths['base_dir'], 'spectra', 'spectrum_files')

def check_and_make_dir(dir_path):
    if not os.path.exists(dir_path):
        os.makedirs(dir_path)

def make_links():
    check_and_make_dir(global_paths['storage_dir'])
    check_and_make_dir(global_paths['output_dir'])
    os.symlink(global_paths['output_dir'], rel_output_dir)
    check_and_make_dir(global_paths['data_dir'])
    os.symlink(global_paths['data_dir'], rel_data_dir)
    check_and_make_dir(global_paths['maps_dir'])
    os.symlink(global_paths['maps_dir'], rel_maps_dir)
    check_and_make_dir(global_paths['spectra_dir'])
    os.symlink(global_paths['spectra_dir'], rel_spectra_dir)

def remove_links():
    os.unlink(rel_output_dir)
    os.unlink(rel_data_dir)
    os.unlink(rel_maps_dir)
    os.unlink(rel_spectra_dir)

if __name__=="__main__":
    parser = argparse.ArgumentParser(description="Make or severe symlinks to storage directories")
    parser.add_argument('action', type=str, choices=['link', 'unlink'], help='Action to take. \'link\' OR \'unlink\'')
    in_args = parser.parse_args()
    if in_args.action == "link":
        make_links()
    else:
        remove_links()
