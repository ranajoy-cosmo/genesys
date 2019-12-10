import numpy as np
import healpy as hp
import h5py
import os
import math
from genesys import Genesys_Class
from genesys import add_method

#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
#* Path naming routines
#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*

class Data_IO(Genesys_Class):
    """
    Structure of the simulation data
    The parent directory containing all simulation and analysis runs are in <global_paths.output_dir>
    <sim_dir> : Parent directory for a single simulation and analysis run
        - <recon_maps_dir> : Could be several different directories depending on how the maps are made
        - <tod_dir> : One single for this set of simulations
            - <channel_dir> : For each frequency channel. Also contains arrays and metadata that is common to all the detectors in the band.
                - metadata_file
                - <detector_dir> : Individual detector in the focal plane. Also contains metadata that is common to all the data chunks
                    - segment_files : HDF5 files containing all the segmented TOD data
    Segment naming convention:
        - Segment 1 starts with time stamp 0.0
        - There are 7 charachters in the segment name, with empty places filled with 0. Ex. segment 42 will be named 0000042
    """
    def __init__(self, sim_params=None):
        self.paths = {}
        if sim_params is not None:
            self.paths['sim_dir'] = os.path.join(self.global_paths['output_dir'], sim_params['sim_tag'])
            self.paths['tod_dir'] = os.path.join(self.paths['sim_dir'], sim_params['tod_tag'])
            if 'recon_tag' in sim_params:
                self.paths['recon_dir'] = os.path.join(self.global_paths['sim_dir'], sim_params['recon_tag'])

    def set_path_for_channel(self, channel_name):
        self.paths['channel_dir'] = os.path.join(self.paths['tod_dir'], channel_name)

    def set_path_for_detector(self, detector_name):
        self.paths['detector_dir'] = os.path.join(self.paths['channel_dir'], detector_name)

    def get_path_to_segment_file(self, segment, fill_size=7):
        segment_name = self.get_segment_name(segment, fill_size=fill_size)
        return os.path.join(self.paths['detector_dir'], segment_name + '.hdf5')

    #*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
    # SEGMENT NAMING CONVENTION
    #*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
    def get_segment_name(self, segment, fill_size=7):
        """
        * COPY OF DESCRIPTION IN CLASS DOC *
        Segment naming convention:
            - Segment 1 starts with time stamp 0.0
            - There are 7 charachters in the segment name, with empty places filled with 0. Ex. Segment 42 will be named 0000042
        """
        return str(segment).zfill(fill_size)

    #*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
    #* FILE EXISTENCE
    #*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
    def check_if_files_exist(segment_dir, file_list):
        exists = True
        for file_name in file_list:
            if not os.path.exists(os.path.join(segment_dir, file_name)):
                exists = False
                break

        return exists

    #*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
    #* MAKING NEW DATA DIRECTORIES
    #*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*

    def make_top_data_directories(self, dir_list, verbose=False):
        """
        This function will create the top level directories under which the simulated data will be written to.
        This function is, under normal circumstances, only called by rank 0, and a barrier is set so that all necessary directories are generated before simulation begins.
        Known issues: Rarely, parallel and independent runs of the code might conflict leading to a directory not being written. Report if this occurs.
        """
        if "sim_dir" in dir_list:
            if verbose:
                self.prompt("Simulation directory: {}".format(iself.paths['sim_dir']))
            self.make_directory(self.paths['sim_dir'], exist_ok=True)

        if "tod_dir" in dir_list:
            if verbose:
                self.prompt("TOD directory: {}".format(iself.paths['tod_dir']))
            self.make_directory(self.paths['tod_dir'], exist_ok=True)

        if "recon_dir" in dir_list:
            map_segment_dir = os.path.join(self.paths['recon_dir'], "sky_map_segments")
            hitmap_segment_dir = os.path.join(self.paths['recon_dir'], "hitmap_segments")
            cov_matrix_segment_dir = os.path.join(self.paths['recon_dir'], "covariance_matrix_segments")
            inv_cov_matrix_segment_dir = os.path.join(self.paths['recon_dir'], "inverse_covariance_matrix_segments")
            if verbose:
                prompt("Recon directory: {}\n".format(self.paths['recon_dir']))
                prompt("Map-segment directory: {}\n".format(map_segment_dir))
                prompt("Hitmap-segment directory: {}\n".format(hitmap_segment_dir))
                prompt("Covariance-matrix-segment directory: {}\n".format(cov_matrix_segment_dir))
                prompt("Inverse-covariance-matrix-segment directory: {}\n".format(inv_cov_matrix_segment_dir))
            self.make_directory(self.paths['recon_dir'], exist_ok=True)
            self.make_directory(map_segment_dir, exist_ok=True)
            self.make_directory(hitmap_segment_dir, exist_ok=True)
            self.make_directory(cov_matrix_segment_dir, exist_ok=True)
            self.make_directory(inv_cov_matrix_segment_dir, exist_ok=True)

    def make_channel_directory(self, verbose=False):
        self.make_directory(self.paths['channel_dir'], exist_ok=True)

    def make_detector_directory(self, verbose=False):
        self.make_directory(self.paths['detector_dir'], exist_ok=True)

    #*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*


@add_method(Data_IO)
def make_directory(dir_name, exist_ok=False):
    """
    Make a new directory if it does not exist already
    If exists_ok is False, an OSErroris raised if the directory already exists
    """
    if not os.path.exists(dir_name):
        os.makedirs(dir_name, exist_ok=exist_ok)
