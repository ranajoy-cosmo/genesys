import h5py
import os
import math
import numpy as np
import healpy as hp
from genesys.gen_io import huffman
from genesys import Genesys_Class

#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
#* Path naming routines
#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*

class GenIO(Genesys_Class):
    """
    The parent directory containing all simulation, map-making and analysis runs are in <global_paths.output_dir>
    The simulation data is written in a HDF5 format whereas maps are written in the HEALPix fits format
    The simulation data I/O is implemented using parallell HDF5 when mpi is available
    Structure of the simulation/map-making/analysis directory:
    <sim_dir> : Parent directory for a single simulation set and map-making/analysis runs
        - <recon_maps_dir> : Could be several different directories depending on how the maps are made
            - <recon_map.fits> :
            - <cov_map.fits> :
            - <hit_map.fits> :
        - <tod_dir> : One single set for each simulation
            - <channel_name>_<segment_name>.h5
    The structure of the HDF5 file
    <channel_name>_<data_block_name>.h5
        - common : (GROUP) Contains common data values for the data block and all detectors 
            - fsamp : (DSET) Sampling frequency
            - segl : (DSET) Segment length in seconds
            - num_segments : (DSET) Number of segments in the data block
            - det : (DSET) List of detector names
            - polang : (DSET) Dictionary {<det_name>: polang_value, ...}. Relative orientation of polarisation axis.
            - coords : (DSET) Coordinate system
            - vsun : (DSET) Solar orbital velocity.
        - segment : (GROUP) The TODs from the different detector are here
            - common : (GROUP) Contains common data values for the segment and all detectors
                - time_0 : (DSET) The first time step of the segment
            - <detector_name> : (GROUP) The unique TOD for the particular detector
                - signal : (DSET) The signal observed by the detector
                - theta : (DSET) The pointing theta
                - phi : (DSET) The pointing phi
                - psi : (DSET) Discretised polarisation angle values.
                - noise : (DSET) The noise independently
                - scalars : (GROUP) Signal properties. [gain, sigma0, fknee, alpha]
                    - gain : (DSET)
                    - sigma0 : (DSET)
                    - fKnee : (DSET)
                    - alpha : (DSET)
    Data block and segment naming convention:
        - Data block 1 starts with time stamp 0.0
        - There are 5 charachters in the data block name, with empty places filled with 0. ex. data block 42 will be named 00042
        - Segments are numbered consecutively, having 7 charachters with empty places filled with 0. ex. data block 5 with 10 segments per block will have segments from 0000041-0000050
    """
    def __init__(self, sim_params, verbosity=0, huffman=False):
        self.paths = {}
        self.verbosity = verbosity
        self.copy_params(sim_params)
        self.paths['sim_dir'] = os.path.join(self.global_paths['output_dir'], sim_params['sim_tag'])
        if huffman:
            self.paths['tod_dir'] = os.path.join(self.paths['sim_dir'], 'tod_huff')
        else:
            self.paths['tod_dir'] = os.path.join(self.paths['sim_dir'], 'tod')
        if 'recon_tag' is not None:
            self.paths['recon_dir'] = os.path.join(self.paths['sim_dir'], sim_params['recon_tag'])

    #*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
    # File and segment naming convention
    #*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*

    def get_data_block_name(self, data_block):
        """
        See data block naming convention in the class docstring
        """
        fill_size = 5
        return str(data_block).zfill(fill_size)

    def get_segment_name(self, segment):
        """
        See segment naming convention in the class docstring
        """
        fill_size = 7
        return str(segment).zfill(fill_size)

    def get_tod_file_path(self, channel_name, data_block):
        """
        For description of tod naming convention read the class docstring
        """
        data_block_name = self.get_data_block_name(data_block)
        tod_file_name = channel_name + '_' + data_block_name + '.h5'
        tod_file_path = os.path.join(self.paths['tod_dir'], tod_file_name)
        return tod_file_path

    def get_map_file_path(self, map_name):
        """
        For decription of types of maps read the class docstring.
        """
        map_file_name = map_name + '.fits'
        map_file_path = os.path.join(self.paths['recon_dir'], map_file_name)
        return map_file_path

    #*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
    # Open/Close TOD file
    #*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*

    def open_tod_file(self, channel_name, data_block, mode):
        tod_file_path = self.get_tod_file_path(channel_name, data_block)
        self.f_tod = h5py.File(tod_file_path, mode)

    def close_tod_file(self):
        self.f_tod.close()

    #*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
    # Reading/Writing TOD file 
    #*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*

    def write_channel_common(self, channel_commons, attributes=None):
        """
        channel_common_data is a dictionary of the vlaues to write to file
        """
        #  det = np.string_(', '.join(channel_commons['det_list']))
        prefix = '/common/'
        for item in channel_commons.keys():
            self.f_tod.create_dataset(prefix + item, data=channel_commons[item])
        if attributes != None:
            for item in attributes.keys():
                for attr_item in attributes[item]:
                    self.f_tod[prefix + item].attrs[attr_item] = attributes[item][attr_item]

    def write_segment_common(self, segment, segment_commons, attributes=None):
        """
        channel_common_data is a dictionary of the vlaues to write to file
        """
        prefix = os.path.join(self.get_segment_name(segment), 'common', '')
        for item in segment_commons.keys():
            self.f_tod.create_dataset(prefix + item, data=segment_commons[item])
        if attributes != None:
            for item in attributes.keys():
                for attr_item in attributes[item]:
                    self.f_tod[prefix + item].attrs[attr_item] = attributes[item][attr_item]

    def write_detector_data(self, segment, detector_name, tod, scalars=None, tod_write_field=['signal', 'theta', 'phi', 'psi']):
        prefix = os.path.join(self.get_segment_name(segment), detector_name, '')
        for item in tod_write_field:
            self.f_tod.create_dataset(prefix + item, data=tod[item])
        if scalars != None:
            prefix = os.path.join(self.get_segment_name(segment), detector_name, 'scalars', '')
            for item in scalars.keys():
                self.f_tod.create_dataset(prefix +item, data=scalars[item])

    def read_channel_common(self, data_fields):
        prefix = '/common/'
        channel_commons = {}
        for item in data_fields:
            channel_commons[item] = self.f_tod[prefix + item][()] 
        return channel_commons

    def read_segment_common(self, segment, data_fields):
        prefix = os.path.join(self.get_segment_name(segment), 'common', '')
        segment_commons = {}
        for item in data_fields:
            segment_commons[item] = self.f_tod[prefix + item][()] 
        return segment_commons

    def read_detector_data(self, segment, detector_name, tod_fields):
        prefix = os.path.join(self.get_segment_name(segment), detector_name, '')
        tod = {}
        for item in tod_fields:
            tod[item] = self.f_tod[prefix + item][:]
        return tod

    def read_detector_scalars(self, segment, detector_name, scalar_fields):
        prefix = os.path.join(self.get_segment_name(segment), detector_name, 'scalars', '')
        scalars = {}
        for item in scalar_fields:
            scalars[item] = self.f_tod[prefix + item][()] 
        return scalars

    def write_map(self, sky_map, map_name):
        map_file_path = self.get_map_file_path(map_name)
        hp.write_map(map_file_path, sky_map)

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
    #* Making new data directories
    #*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*

    def make_directory(self, dir_name, exist_ok=False):
        """
        Make a new directory if it does not exist already
        If exists_ok is False, an OSErroris raised if the directory already exists
        """
        if not os.path.exists(dir_name):
            os.makedirs(dir_name, exist_ok=exist_ok)

    def make_top_data_directories(self, dir_list):
        """
        This function will create the top level directories under which the simulated data will be written to.
        This function is, under normal circumstances, only called by rank 0, and a barrier is set so that all necessary directories are generated before simulation begins.
        """
        if self.verbosity > 1:
            self.prompt(f"Simulation directory: {self.paths['sim_dir']}")
        self.make_directory(self.paths['sim_dir'], exist_ok=True)

        if "tod_dir" in dir_list:
            if self.verbosity > 1:
                self.prompt(f"TOD directory: {self.paths['tod_dir']}")
            self.make_directory(self.paths['tod_dir'], exist_ok=True)

        if "recon_dir" in dir_list:
            if self.verbosity > 1:
                self.prompt(f"Recon directory: {self.paths['recon_dir']}")
            self.make_directory(self.paths['recon_dir'], exist_ok=True)

    #*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
    #* Huffman compression
    #*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
    """
    The structure of the Huffman encoded HDF5 file
    <channel_name>_<segment_name>.h5
        - common : (GROUP) Contains common data values for the segment and all detectors 
            - fsamp : (DSET) Sampling frequency
            - nside : (DSET) NSide of HEALPix map used for Huffman compressing pointing
            - npsi : (DSET) Discretisation over 2pi for Huffman compressing psi
            - det : (DSET) List of detector names
            - polang : (DSET) Dictionary {<det_name>: polang_value, ...}. Relative orientation of polarisation axis.
        - segment : (GROUP) The TODs from the different detector are here
            - common : (GROUP) The common TODs for the entire detector set
                - time : (DSET) The time steps. (*Huffman)
                    - unit : (ATTR) Time unit and convention
                - vsun : (DSET) Solar orbital velocity. (*Huffman)
                    - info : (ATTR) axes
                    - coords : (ATTR) Coordinate system
                - hufftree : (DSET)
                - huffsymb : (DSET)
            - <detector_name> : (GROUP) The unique TOD for the particular detector
                - signal : (DSET) The signal observed by the detector
                - flag : (DSET) Flagging of invalid data by 0. (*Huffman)
                - pixels : (DSET) Discretised pointing at a given HEALPix NSide. (*Huffman)
                - psi : (DSET) Discretised polarisation angle values. (*Huffman)
                - scalars : Array of gain and noise properties. [gain, sigma0, fknee, alpha]
                    - legend : Description of the scalar elements as a ',' separated string.
    """
