import h5py
import os
import math
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
    The structure of the Huffman encoded HDF5 file
    <channel_name>_<segment_name>.h5
        - common : (GROUP) Contains common data values for the segment and all detectors 
            - fsamp : (DSET) Sampling frequency
            - nside : (DSET) NSide of HEALPix map used for Huffman compressing pointing
            - npsi : (DSET) Discretisation over 2pi for Huffman compressing psi
            - npsi_hwp : (*Optional) (DSET) Discretisation over 2pi for Huffman compressing psi_hwp
            - det : (DSET) List of detector names
            - polang : (DSET) Dictionary {<det_name>: polang_value, ...}. Relative orientation of polarisation axis.
        - tod : (GROUP) The TODs from the different detector are here
            - common : (GROUP) The common TODs for the entire detector set
                - time : (DSET) The time steps. (*Huffman)
                    - unit : (ATTR) Time unit and convention
                - vsun : (DSET) Solar orbital velocity. (*Huffman)
                    - info : (ATTR) axes
                    - coords : (ATTR) Coordinate system
                - psi_hwp : (*Optional) (DSET) Discretised HWP-polarisation angle values. (*Huffman)
                - hufftree : (DSET)
                - huffsymb : (DSET)
            - <detector_name> : (GROUP) The unique TOD for the particular detector
                - signal : (DSET) The signal observed by the detector
                - flag : (DSET) Flagging of invalid data by 0. (*Huffman)
                - pixels : (DSET) Discretised pointing at a given HEALPix NSide. (*Huffman)
                - psi : (DSET) Discretised polarisation angle values. (*Huffman)
                - scalars : Array of gain and noise properties. [gain, sigma0, fknee, alpha]
                    - legend : Description of the scalar elements as a ',' separated string.
    Segment naming convention:
        - segment 1 starts with time stamp 0.0
        - there are 7 charachters in the segment name, with empty places filled with 0. ex. segment 42 will be named 0000042
    """
    def __init__(self, sim_params, verbosity=0):
        self.paths = {}
        self.verbosity = verbosity
        self.copy_params(sim_params, ['sim_tag', 'tod_tag', 'recon_tag', 'special_tag', 'data_products', 'huffman'])
        self.paths['sim_dir'] = os.path.join(self.global_paths['output_dir'], sim_params['sim_tag'])
        self.paths['tod_dir'] = os.path.join(self.paths['sim_dir'], 'tod')
        if 'recon_tag' is not None:
            self.paths['recon_dir'] = os.path.join(self.paths['sim_dir'], sim_params['recon_tag'])

    #*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
    # File and segment naming convention
    #*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*

    def get_segment_name(self, segment):
        """
        See segment naming convention in the class docstring
        """
        fill_size = 7
        return str(segment).zfill(fill_size)

    def get_tod_file_path(self, channel_name, segment):
        """
        For description of tod naming convention read the class docstring
        """
        segment_name = self.get_segment_name(segment)
        tod_file_name = channel_name + '_' + segment_name + '.h5'
        tod_file_path = os.path.join(self.paths['tod_dir'], tod_file_name)
        return tod_file_path

    def get_map_file_path(self, map_type):
        """
        For decription of types of maps read the class docstring.
        """
        valid_map_types = ['recon_map', 'cov_map', 'hit_map']
        assert map_type in valid_map_types, f"'{map_type}' not among valid map types. Valid types: {', '.join(valid_map_types)}"
        map_file_name = map_type + '.fits'
        map_file_path = os.path.join(self.params['recon_dir'], map_file_name)
        return map_file_name

    #*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
    # Reading/writing TOD
    #*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*

    def open_tod_file(self, channel_name, segment, comm=None):
        tod_file_path = self.get_tod_file_path(channel_name, segment)
        if comm == None:
            self.f_tod = h5py.File(file_name, 'w')
        else:
            self.f_tod = h5py.File(file_name, 'w', driver='mpio', comm=comm)

    def close_tod_file(self):
        self.f_tod.close()

    def write_channel_common(self, channel_common_data):
        """
        channel_common_data is a dictionary of the vlaues to write to file
        """
        prefix = '/'

    def read_tod_from_file(self, tod, segment, data_fields, comm=None):
        file_name = self.get_path_to_segment_file(segment)
        if comm == None:
            f = h5py.File(file_name, 'r')
        else:
            f = h5py.File(file_name, 'r', driver='mpio', comm=comm)
        for item in data_fields:
            tod.__dict__[item] = f[item][:]
        f.close()

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
