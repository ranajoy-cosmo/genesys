"""
This module contains the standardised I/O routines for the genesys environment
The I/O is done using parallel HDF5
"""
import numpy as np
import healpy as hp
import h5py
import os
import math
from genesys.global_config import global_paths
from genesys.utilities import prompt

#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
#* Path naming routines
#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
def get_path_to_parent_data_set_dir(config):
    """
    Top level directory containing:
        Simulated timestream data HDF5 file: sim_data.dhf5
        Reconstructed maps folders: rec_<rec_tag>
        Run output folder
    Convention:
        <global_paths.output_dir>/<config.data_set_tag>
    """
    return os.path.join(global_paths.output_dir, config.data_set_tag)

def get_path_to_simulation_data_file(config):
    """
    Simulated timestream data HDF5 file
    Convention:
        <global_paths.output_dir>/<config.data_set_tag>/sim_data.hdf5
    """
    return os.path.join(get_path_to_parent_data_set_dir(config), sim_data + '.hdf5')

def get_path_to_reconstructed_maps_directory(config):
    """
    Path to a specific set of reconstructed maps
    Convention:
        <global_paths.output_dir>/<config.data_set_tag>/rec_<rec_tag>
    """
    return os.path.join(global_paths.output_dir, "rec_" + config.rec_tag)

def get_path_to_reconstructed_map_file(config, map_type):
    """
    Path to a specific reconstructed map.
    Convention:
        Sky signal map: sky_map.fits
        Covariance map: covariance_map.fits
        Inverse-covariance map: inverse_covariance_map.fit
        Hit-map: hit_map.fits
    """
    return os.path.join(get_path_to_reconstructed_maps_directory(config), map_type + ".fits") 

#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
# HDF5 objects
#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*

def get_simulation_data_file_object(config, mode='a'):
    sim_data_file_name = get_path_to_simulation_data_file(config)
    return h5py.File(sim_data_file_name, mode=mode)

def get_band_group(sim_df_obj, band_name, mode='a'):

#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
# Band and Detector naming convention
def get_band_id
def get_detector_id(band_name, detector_name):
    return band_name + '--' + detector_name
def get_detector_a_name(detector_name):
    return detector_name + 'a'
def get_detector_b_name(detector_name):
    return detector_name + 'b'
#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*

def combine_return_and_write_stream(return_field, write_field):
    return dict.fromkeys(list(set().union(return_field + write_field)))

def get_return_stream(return_field):
    return dict.fromkeys(return_field)


#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
#* File existence
#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
def check_if_files_exist(segment_dir, file_list):
    exists = True
    for file_name in file_list:
        if not os.path.exists(os.path.join(segment_dir, file_name)):
            exists = False
            break

    return exists
#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*

def make_directory(dir_name, exist_ok=False):
    if not os.path.exists(dir_name):
        os.makedirs(dir_name, exist_ok=exist_ok)

def make_top_data_directories(config, dir_list, verbose=False):
    """
    This function will create the top level directories under which the simulated data will be written to.
    This function is, under normal circumstances, only called by rank 0, and a barrier is set so that all necessary directories are generated before simulation begins.
    Known issues: Rarely, parallel and independent runs of the code might conflict leading to a directory not being written. Report if this occurs.
    """
    # The main sim directory under which all simulated timestream data and corresponding maps will be written to.
    if "sim_dir" in dir_list:
        sim_dir = get_path_to_sim_dir(config)
        if verbose:
            prompt("Simulation directory: {}\n".format(sim_dir))
        make_directory(sim_dir, exist_ok=True)

    # The scan directory under which all simulated timestream data will be written to.
    if "sim_dir" in dir_list:
        scan_dir = get_path_to_scan_dir(config)
        if verbose:
            prompt("Scan directory: {}\n".format(scan_dir))
        make_directory(scan_dir, exist_ok=True)

    # The recon directory under which the reconstructed map and map products will be stored
    if "recon_dir" in dir_list:
        recon_dir = get_path_to_recon_dir(config)
        map_segment_dir = os.path.join(recon_dir, "sky_map_segments")
        hitmap_segment_dir = os.path.join(recon_dir, "hitmap_segments")
        cov_matrix_segment_dir = os.path.join(recon_dir, "covariance_matrix_segments")
        inv_cov_matrix_segment_dir = os.path.join(recon_dir, "inverse_covariance_matrix_segments")
        if verbose:
            prompt("Recon directory: {}\n".format(recon_dir))
            #  prompt("Map-segment directory: {}\n".format(map_segment_dir))
            #  prompt("Hitmap-segment directory: {}\n".format(hitmap_segment_dir))
            #  prompt("Covariance-matrix-segment directory: {}\n".format(cov_matrix_segment_dir))
            #  prompt("Inverse-covariance-matrix-segment directory: {}\n".format(inv_cov_matrix_segment_dir))
        make_directory(recon_dir, exist_ok=True)
        #  make_directory(map_segment_dir, exist_ok=True)
        #  make_directory(hitmap_segment_dir, exist_ok=True)
        #  make_directory(cov_matrix_segment_dir, exist_ok=True)
        #  make_directory(inv_cov_matrix_segment_dir, exist_ok=True)

    prompt("\n")
        
def make_detector_directory(config):
    detector_path = get_path_to_detector_dir(config)
    make_directory(detector_path, exist_ok=True)

def make_segment_directory(config, segment):
    segment_path = get_path_to_segment_dir(config, segment)
    make_directory(segment_path, exist_ok=True)

#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*

def write_segment_data(config, segment, data, data_name):
    segment_path = get_path_to_segment_dir(config, segment)
    if config.special_tag == "":
        data_special_name = data_name
    else:
        data_special_name = config.special_tag + '_' + data_name
    np.save(os.path.join(segment_path, data_special_name), data)

def read_segment_data(config, segment, data_name):
    segment_path = get_path_to_segment_dir(config, segment)
    if config.special_tag == "":
        data_special_name = data_name
    else:
        data_special_name = config.special_tag + '_' + data_name
    data = np.load(os.path.join(segment_path, data_special_name + '.npy'))
    return data

#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*

def write_sky_map(config, sky_map, map_name):
    recon_dir = get_path_to_recon_dir(config)
    map_file_name = os.path.join(recon_dir, map_name + '.fits')
    hp.write_map(map_file_name, sky_map)
