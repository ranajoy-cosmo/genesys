import numpy as np
import healpy as hp
import h5py
import os
import math
from .. import Genesys_Class
from genesys.global_config import global_paths
from genesys.utilities import prompt

#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
#* Path naming routines
#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*

class Path(Genesys_Class):
    """
    Structure of the simulation data
    The parent directory containing all simulation and analysis runs are in <global_paths.output_dir>
    <sim_dir> : Parent directory for a single simulation and analysis run
        - <recon_maps_dir> : Could be several different directories depending on how the maps are made
        - <tod_dir> : One single for this set of simulations
            - <band_dir> : For each frequency channel. Also contains arrays and metadata that is common to all the detectors in the band.
                - metadata_file
                - <detector_dir> : Individual detector in the focal plane. Also contains metadata that is common to all the data chunks
                    - segment_files : HDF5 files containing all the segmented TOD data

    The name tags are provided in the config file passed to the simulation run.
    sim_tag
    tod_tag
    recon_tag
    band_name
    detector_name
    segment_name
    """
    def get_path_to_sim_dir(self, config):
        return os.path.join(self.global_paths['data_dir'], config.sim_tag)

    def get_path_to_tod_dir(self, config):
        return os.path.join(self.get_path_to_sim_dir(config), config.tod_tag)
     
    def get_path_to_recon_maps_dir(self, config):
        return os.path.join(self.get_path_to_sim_dir(config), config.recon_tag)

    def get_path_to_band_dir(config):
        return os.path.join(get_path_to_tod_dir(config), 

    def get_path_to_detector_dir(config):
        return os.path.join(global_paths.output_dir, config.sim_tag, config.scan_tag, config.detector_id)

    def get_path_to_segment_dir(config, segment):
        """
        The segment name is given as an integer. It is converted by this routine to the appropriate 4 digit string preceded by 0s.
        By convention, 1 is added to the segment name. So, the 0th segment will be named '0001'.
        """
        segment_name = get_segment_name(segment)
        return os.path.join(global_paths.output_dir, config.sim_tag, config.scan_tag, config.detector_id, segment_name)

#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
# Detector naming convention
def get_detector_id(band_name, detector_name):
    return band_name + '--' + detector_name
def get_detector_a_name(detector_name):
    return detector_name + 'a'
def get_detector_b_name(detector_name):
    return detector_name + 'b'
# Segment naming convention
def get_segment_name(segment):
    return str(segment+1).zfill(4)
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
