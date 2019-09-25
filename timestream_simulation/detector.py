import os
import importlib
import sys
import numpy as np
import healpy as hp
from termcolor import colored
from itertools import compress
from .pointing.pointing import Pointing
from ..noise.noise import Noise 
from ..utilities import Generic_Class, prompt
from ..utilities import unit_conversion as uc
from ..global_config import global_paths
from ..data_io import data_io as dio

"""
This module contains the Bolo class. The Bolo class is a generator of instances of individual bolometers with their independent properties given in the config file. The member functions of the class generate timestream signals for the particular bolometer configuration and 
"""

#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
# The Detector Class 
#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*

class Detector:

    def __init__(self, band_name, detector_name, config, sim_run=True):
        # Loading the config files to self
        self.config = Generic_Class()
        self.band_config = Generic_Class()
        self.band_config = importlib.import_module("genesys.instrument." + config.focal_plane_config).focal_plane.bands[band_name]
        detector_config = self.band_config.detectors[detector_name]
        self.config.__dict__.update(config.__dict__)
        self.config.__dict__.update(detector_config.__dict__)
        # Setting up the derived parameters and the local directories
        self.setup_derived_params()
        self.noise = Noise(self.config)
        self.pointing = Pointing(self.config)
        if sim_run:
            dio.make_detector_directory(self.config)
            self.get_sky_map()

    #*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
    # Reading the timestream data for the detector that was already simulated. 
    #*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*

    def read_timestream_data(self, segment, return_field, noise_only=False):
        time_streams_dict = {}

        if "signal" in return_field:
            if noise_only:
                time_streams_dict["signal"] = dio.read_segment_data(self.config, segment[0], "noise")
            else:
                time_streams_dict["signal"] = dio.read_segment_data(self.config, segment[0], "signal")
            return_field.remove("signal")
        for return_item in return_field:
            time_streams_dict[return_item] = dio.read_segment_data(self.config, segment[0], return_item)

        return time_streams_dict

    #*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
    # Simulating the time-ordered signal data for a given bolo with any beam
    #*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
    def simulate_timestream_data(self, segment, return_field=[]):
        # Initialising the pointing object for the present segment
        self.pointing.initialise_for_segment(segment)

        # Making the directories where the segment data will be stored
        dio.make_segment_directory(self.config, segment[0])

        time_streams_dict = {}
        write_field = self.config.ts_data_products

        #Generating the pointing and orientation vectors
        pointing_vec = self.pointing.get_vec_obv()
        if "pointing_vec" in write_field:
            dio.write_segment_data(self.config, segment[0], pointing_vec, "pointing_vec") 
        hitpix = hp.vec2pix(self.config.nside_in, pointing_vec[...,0], pointing_vec[...,1], pointing_vec[...,2])

        pol_ang = self.pointing.get_pol_ang(pointing_vec)
        if "pol_ang" in write_field:
            dio.write_segment_data(self.config, segment[0], pol_ang, "pol_ang") 
        if "pointing_vec" not in return_field:
            del pointing_vec
        else:
            time_streams_dict['pointing_vec'] = pointing_vec

        # Generating the signal (without noise)
        signal = self.data_model(hitpix, pol_ang)
        if "pol_ang" not in return_field:
            del pol_ang
        else:
            time_streams_dict['pol_ang'] = pol_ang

        del hitpix

        if self.config.noise_type != "none":
            self.noise.initialise_for_segment(segment)
            noise = self.noise.simulate_ts_noise()
            signal += noise 
            if "noise" in write_field:
                dio.write_segment_data(self.config, segment[0], noise, "noise") 
            if "noise" not in return_field:
                del noise
            else:
                time_streams_dict['noise'] = noise


        if "signal" in write_field:
            dio.write_segment_data(self.config, segment[0], signal, "signal") 
        if "signal" not in return_field:
            del signal
        else:
            time_streams_dict['signal'] = signal

        if "time" in write_field + return_field:
            time = self.pointing.get_t_steps()
            if "time" in write_field:
                dio.write_segment_data(self.config, segment[0], time, "time") 
            if "time" not in return_field:
                del time
            else:
                time_streams_dict['time'] = time

        return time_streams_dict

    def data_model(self, hitpix, pol_ang):
        if self.config.polarisation_modulation == "instrumental_scanning":
            pol_phase_factor = 2
        else:
            pol_phase_factor = 4

        if self.config.sim_pol_type == "noise_only":
            signal = np.zeros(self.nsamples - 2*self.pointing.pad) 
        if self.config.sim_pol_type == "I":
            signal = self.sky_map[hitpix]
        if self.config.sim_pol_type == "QU":
            signal = self.sky_map[0][hitpix]*np.cos(pol_phase_factor*pol_ang) + self.sky_map[1][hitpix]*np.sin(pol_phase_factor*pol_ang)
        if self.config.sim_pol_type == "_QU":
            signal = self.sky_map[1][hitpix]*np.cos(pol_phase_factor*pol_ang) + self.sky_map[2][hitpix]*np.sin(pol_phase_factor*pol_ang)
        if self.config.sim_pol_type == "IQU":
            signal = self.sky_map[0][hitpix] + self.sky_map[1][hitpix]*np.cos(pol_phase_factor*pol_ang) + self.sky_map[2][hitpix]*np.sin(pol_phase_factor*pol_ang)
        return signal

#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*

    #Calculating additional scan parameters
    def setup_derived_params(self):
        #Cross scan resolution
        self.config.theta_cross = 360.0*60.0*np.sin(np.radians(self.config.alpha))*self.config.t_spin/self.config.t_prec
        #Co scan resolution
        total_beta_pos = np.radians(self.config.beta) + np.radians(self.config.focal_plane_pos[1]/60.0) + np.radians(self.config.offset[1]/60.0/60.0)
        self.config.theta_co = 360*60*np.sin(total_beta_pos)/self.config.sampling_rate/self.config.t_spin
        self.config.scan_resolution = self.config.theta_co


    def get_sky_map(self):
        if self.config.sim_pol_type == "noise_only":
            self.sky_map = None
        elif self.config.sim_pol_type == "T":
            self.sky_map = hp.read_map(self.config.input_sky_map, verbose=False)
        elif self.config.sim_pol_type == "QU":
            self.sky_map = np.array(hp.read_map(self.config.input_sky_map, field=(0,1), verbose=False))
        elif self.config.sim_pol_type == "_QU":
            self.sky_map = np.array(hp.read_map(self.config.input_sky_map, field=(1,2), verbose=False))
        else:
            self.sky_map = np.array(hp.read_map(self.config.input_sky_map, field=(0,1,2), verbose=False))

        self.config.nside_in = hp.get_nside(self.sky_map)

    #*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*

    def segment_start_prompt(self, rank, segment):
        t_start = segment[1]
        t_stop = segment[2]
        t_duration = t_stop - t_start
        display_string = colored("#* Rank {} doing Detector {} and Segment {}:\n".format(rank, self.config.detector_id, segment[0]), color="blue")
        display_string += "\nSTART TIME: {}s, {}d".format(t_start, uc.convert_unit('time', t_start, 'second', 'day'))
        display_string += "\nSTOP TIME: {}s, {}d".format(t_stop, uc.convert_unit('time', t_stop, 'second', 'day'))
        display_string += "\nSCAN DURATION: {}s, {}d".format(t_duration, uc.convert_unit('time', t_duration, 'second', 'day'))
        display_string += colored("#*#*#*\n\n", color="blue")

    def display_params(self):
        display_string = colored("#* Detector {} parameters\n".format(self.config.detector_id), color="green")
        display_string += "IDEAL BAND PARAMETERS ***\n"
        display_string += "\tBAND NAME: {}\n".format(self.band_config.band_id)
        display_string += "\tBAND CENTRAL FREQUENCY: {}\n".format(self.band_config.band_frequency)
        display_string += "\tBAND_WIDTH: {}\n".format(self.band_config.band_width)
        display_string += "\tBEAM FWHM: {}\n".format(self.band_config.beam_fwhm)
        display_string += "\tDETECTOR NET: {}\n".format(self.band_config.det_NET)
        display_string += "\tDETECTOR YIELD: {}\n".format(self.band_config.detector_yield)
        display_string += "\tNOISE MARGIN: {}\n".format(self.band_config.noise_margin)
        display_string += "***\n"
        display_string += "POSITION ON FOCAL PLANE: {} arcmins\n".format(self.config.focal_plane_pos)
        display_string += "OFFSET FROM POSITION: {} arcsecs\n".format(self.config.offset)
        display_string += "SAMPLING RATE: {} Hz\n".format(self.config.sampling_rate)
        display_string += "THETA CO: {} arcmin\n".format(self.config.theta_co)
        display_string += "THETA CROSS: {} arcmin\n".format(self.config.theta_cross)
        display_string += "SCAN RESOLUTION: {} arcmin\n".format(self.config.scan_resolution)
        display_string += "BEAM FWHM MAJOR: {}\n".format(self.config.beam_fwhm_major)
        display_string += "BEAM FWHM MINOR: {}\n".format(self.config.beam_fwhm_minor)
        if self.config.polarisation_modulation == "instrumental_scanning":
            display_string += "INITIAL POLARISATION ANGLE: {}\n".format(self.config.pol_phase_ini)
        display_string += "INPUT SKY MAP: {}\n".format(self.config.input_sky_map)
        if self.config.noise_type != "none":
            display_string += "WHITE NOISE RMS: {}\n".format(self.config.white_noise_sigma)
            if self.config.noise_type == "1_over_f":
                display_string += "f-knee: {} mHZ\n".format(self.config.f_knee)
                display_string += "1-over-f alpha: {}\n".format(self.config.one_over_f_alpha)
        display_string += colored("#*#*#*\n\n", color="green")
        prompt(display_string, sys.stdout)

    #*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*


#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
# This section conatains code for checking parameters. Will fix it later and add
    #  def make_check(self):
        #  check_flag = 1
#
        #  # Checking main config variables
        #  main_var_list = ['sim_tag', 'scan_tag', 'special_tag', 'simulate_data', 'sim_pol_type', 'coordinate_system', 'tod_type', 'ts_data_products']
        #  main_var_list += ['scan_strategy_config', 'focal_plane_config']
        #  main_var_list += ['noise_type', 'beam_type', 'write_beam']
        #  main_var_list += ['t_segment', 'band_detector_segment_dict']
        #  main_var_exists = list(map(lambda var_name: var_name in dir(self.config), main_var_list))
        #  if False in main_var_exists:
            #  check_flag = 0
            #  prompt(colored("WARNING: ", color="red") + "{} is not defined.\n".format(list(compress(main_var_list, np.logical_not(main_var_exists)))))
#
        #  if self.config.simulate_data:
            #  sim_pol_type_list = ['IQU', 'QU', 'I', '_QU', 'noise_only']
            #  if self.config.sim_pol_type not in sim_pol_type_list:
                #  check_flag = 0
                #  prompt(colored("WARNING: ", color="red") + "sim_pol_type is defines as \"{}\". It should be one of {}.\n".format(self.config.sim_pol_type, sim_pol_type_list))
            #  if self.config.sim_pol_type != "noise_only" and "nside_in" not in dir(self.config):
                #  check_flag = 0
                #  prompt(colored("WARNING: ", color="red") + "nside_in is not defined\n")
#
            #  if self.config.tod_type not in ['signal', 'gradient']:
                #  ckeck_flag = 0
                #  prompt(colored("WARNING: ", color="red") + "tod_type is defined as \"{}\". It should be either \"signal\" or \"gradient\".\n".format(self.config.sim_pol_type))
            #  if self.config.tod_type == "gradient":
                #  grad_type_list = ['zero_order', 'grad_co', 'grad_cross', 'grad_coXgrad_co', 'grad_crossXgrad_cross', 'grad_coXgrad_cross']
                #  if not set(self.config.grad_type).issubset(set(grad_type_list)):
                    #  check_flag = 0
                    #  prompt(colored("WARNING: ", color="red") + "grad_type is set to {}. Please check as it is not one of the defined types\n".format(self.config.grad_type))
#
        #  if self.config.coordinate_system not in ['ecliptic', 'galactic']:
            #  check_flag = 0
            #  prompt(colored("WARNING: ", color="red") + "coordinate_system is defined as \"{}\". It should be either \"ecliptic\" or \"galactic\".\n".format(self.config.coordinate_system))
        #  #-----------------------------------------------------------
#
        #  # Checking that the scan strategy variables are in place
#
        #  scan_strategy_var_list = ['t_year', 't_prec', 't_spin', 'sampling_rate', 'oversampling_rate', 'alpha', 'beta', 'scan_strategy_module']
        #  scan_strategy_var_exists = list(map(lambda var_name: var_name in dir(self.config), scan_strategy_var_list))
        #  if False in scan_strategy_var_exists:
            #  check_flag = 0
            #  prompt(colored("WARNING: ", color="red") + "{} is not defined for the scan strategy.\n".format(list(compress(scan_strategy_var_list, np.logical_not(scan_strategy_var_exists)))))
        #  #-----------------------------------------------------------
#
        #  # Checking that the I/O variables are in place
        #  ts_data_product_list = ['signal', 'pointing_vec', 'pol_ang', 'noise', 'mask']
        #  if not set(self.config.ts_data_products).issubset(set(ts_data_product_list)):
            #  check_flag = 0
            #  prompt(colored("WARNING: ", color="red") + "ts_data_products are {}. Please check as it is not one of the defined types\n".format(self.config.ts_data_products))
        #  #-----------------------------------------------------------
#
        #  # Noise and beam settings
        #  beam_type_list = ['pencil', 'full_simulated', 'from_file']
        #  if not self.config.beam_type in beam_type_list:
            #  check_flag = 0
            #  prompt(colored("WARNING: ", color="red") + "beam_type is {}. Please check as it is not one of the defined types\n".format(self.config.beam_type))
        #
        #  if self.config.beam_type == "full_simulated":
            #  if "beam_cutoff" not in dir(self.config):
                #  check_flag = 0
                #  prompt(colored("WARNING: ", color="red") + "beam_cutoff not defined")
#
        #  noise_type_list = ['white', '1_over_f', 'none']
        #  if self.config.noise_type not in noise_type_list:
            #  check_flag = 0
            #  prompt(colored("WARNING: ", color="red") + "noise_type is {}. Please check as it is not one of the defined types\n".format(self.config.noise_type))
        #  #-----------------------------------------------------------
#
        #  self.check_flag = check_flag
#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
