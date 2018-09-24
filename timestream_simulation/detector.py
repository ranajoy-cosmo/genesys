import numpy as np
import healpy as hp
import os
import importlib
import sys
from genesys.beam.beam import Beam
from genesys.timestream_simulation.pointing.pointing import Pointing
from genesys.noise.noise import Noise 
from genesys.utilities import Generic_Class, prompt
from termcolor import colored
from genesys.global_config import global_paths
from itertools import compress
import genesys.data_io.data_io as dio

"""
This module contains the Bolo class. The Bolo class is a generator of instances of individual bolometers with their independent properties given in the config file. The member functions of the class generate timestream signals for the particular bolometer configuration and 
"""

#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
# The Bolo Class 
#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*

class Detector:

    def __init__(self, detector_name, config):
        self.config = Generic_Class()
        detector_config = importlib.import_module("genesys.instrument." + config.focal_plane_config).focal_plane.detectors[detector_name]
        self.config.__dict__.update(config.__dict__)
        self.config.__dict__.update(detector_config.__dict__)
        self.set_dirs()
        self.check_flag = self.make_check()
        self.setup_derived_params()
        if config.simulate_data:
            #  self.beam = Beam(self.config)
            self.noise_class = Noise(self.config)
            self.get_sky_map()
            #  self.beam.write_beam(self.bolo_dir)



#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
# Reading the timestream data for the bolo that was already simulated. 
#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*

    def read_timestream(self, segment, return_field, noise_only=False):
        segment_dir = self.get_segment_dir(segment)
        t_stream = {}

        if "signal" in return_field:
            if noise_only:
                t_stream["signal"] = np.load(os.path.join(segment_dir, "noise.npy"))
            else:
                t_stream["signal"] = np.load(os.path.join(segment_dir, "signal.npy"))
            return_field.remove("signal")
        for return_item in return_field:
            t_stream[return_item] = np.load(os.path.join(segment_dir, return_item + ".npy"))

        return t_stream 

#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
# Simulating the time-ordered signal data for a given bolo with any beam
#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
    def simulate_timestream_signal(self, segment, return_field=[]):
        t_return_stream = dict.fromkeys(return_field)
        write_field = self.config.timestream_data_products
        if write_field:
            self.make_write_dir(segment)

        #Generating the pointing and orientation vectors for the central line
        prompt("0.0\n")
        self.pointing = Pointing(self.config, segment, self.beam.del_beta.size)
        
        pointing_vec = self.pointing.get_vec_obv(0.0)
        hitpix = hp.vec2pix(self.config.nside_in, pointing_vec[...,0], pointing_vec[...,1], pointing_vec[...,2])

        #Generating the polarisation angle of the detector on the sky
        pol_ang = self.pointing.get_pol_ang(pointing_vec)
        if self.config.sim_pol_type == "T" or self.config.sim_pol_type == "noise_only":
            cos2=None
            sin2=None
        else:
            cos2 = np.cos(2*pol_ang)
            sin2 = np.sin(2*pol_ang)

        if "pointing_vec" in write_field:
            self.write_timestream_data(pointing_vec, "pointing_vec", segment)
        if "pointing_vec" in return_field:
            if self.config.beam_type == "pencil":
                t_return_stream["pointing_vec"] = pointing_vec
            else:
                t_return_stream["pointing_vec"] = pointing_vec[self.pointing.pad:-self.pointing.pad][::self.config.oversampling_rate]

        del pointing_vec

        if "pol_ang" in write_field:
            self.write_timestream_data(pol_ang, "pol_ang", segment)
        if "pol_ang" in return_field:
            if self.config.beam_type == "pencil":
                t_return_stream["pol_ang"] = pol_ang
            else:
                t_return_stream["pol_ang"] = pol_ang[self.pointing.pad:-self.pointing.pad][::self.config.oversampling_rate]

        del pol_ang

        if self.config.beam_type == "pencil":
            signal = self.generate_signal_pencil_beam(hitpix, cos2, sin2)
        elif self.config.beam_type in ["full_simulated", "from_file"]:
            beam_kernel_row = self.beam.get_beam_row(0.0)                       #The input argument is the beam offset from the centre
            signal = self.generate_signal_full_beam(hitpix, beam_kernel_row, cos2, sin2)
            #Iterating over the beam map and integrating
            for del_beta in self.beam.del_beta:
                if del_beta == 0.0:
                    continue
                prompt("{}\n".format(del_beta))
                beam_kernel_row = self.beam.get_beam_row(del_beta)
                pointing_vec = self.pointing.get_vec_obv(del_beta)
                hitpix = hp.vec2pix(self.config.nside_in, pointing_vec[...,0], pointing_vec[...,1], pointing_vec[...,2])
                signal += self.generate_signal_full_beam(hitpix, beam_kernel_row, cos2, sin2)
            beam_sum = np.sum(self.beam.beam_kernel[0])
            signal /= beam_sum
        else:
            prompt("Beam type not recognised", sys.stdout)
            sys.exit()
                

        if self.config.noise_type != "none":
            noise = self.noise_class.simulate_timestream_noise_from_parameters()
            if "noise" in write_field:
                self.write_timestream_data(noise, "noise", segment)
            if "noise" in return_field:
                t_return_stream["noise"] = signal

            signal[::self.config.oversampling_rate] += noise 

        if "signal" in write_field:
            self.write_timestream_data(signal, "signal", segment)
        if "signal" in return_field:
            t_return_stream["signal"] = signal[::self.config.oversampling_rate]

        if return_field:
            return t_return_stream

#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
# Simulating the time-ordered gradient data for a given bolo with any beam
#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
    
    def simulate_timestream_template(self, segment):
        self.make_write_dir(segment)

        if self.config.template_type == "tm_bandpass":
            promtp("Patience. Not ready yet\n")
        if self.config.template_type == "tm_gradient":
            signal = dict.fromkeys(["top", "central", "bottom"])

            if list(set(["tm_grad_co", "tm_grad_co_co", "tm_grad_cross_cross"]) & set(self.config.tm_gradient_type)):
                signal["central"] = self.simulate_timestream_signal(segment, return_field=["signal"])["signal"]

            if list(set(["tm_grad_cross", "tm_grad_cross_cross", "tm_grad_co_cross"]) & set(self.config.tm_gradient_type)):
                self.beam.shift_beam_kernel("top")
                signal["top"] = self.simulate_timestream_signal(segment, return_field=["signal"])["signal"]
                self.beam.restore_beam_kernel()
                self.beam.shift_beam_kernel("bottom")
                signal["bottom"] = self.simulate_timestream_signal(segment, return_field=["signal"])["signal"]

            for gradient_type in self.config.tm_gradient_type:
                gradient = self.generate_signal_gradient(gradient_type, signal)
                self.write_timestream_data(gradient, gradient_type, segment)

    #Generating timestream signal, no noise. Option of Intensity only, Polarisation only or all components
    #Although this promotes code duplicity, the signal generation subroutine is done separately for pencil beams, full beams and gradient calculation to not make the single subroutine too congested and difficult to understand and modify
    def generate_signal_pencil_beam(self, hit_pix, cos2, sin2):
        if self.config.sim_pol_type == "noise_only":
            signal = np.zeros(self.nsamples - 2*self.pointing.pad) 

        elif self.config.sim_pol_type == "T":
            signal = 0.5*self.sky_map[hit_pix]

        elif self.config.sim_pol_type in ["QU", "_QU"]:
            signal = 0.5*(self.sky_map[0][hit_pix]*cos2 + self.sky_map[1][hit_pix]*sin2) 

        else:
            signal = 0.5*(self.sky_map[0][hit_pix] + self.sky_map[1][hit_pix]*cos2 + self.sky_map[2][hit_pix]*sin2) 

        return signal

    def generate_signal_full_beam(self, hit_pix, beam_kernel_row, cos2, sin2):
        if self.config.sim_pol_type == "noise_only":
            signal = np.zeros(self.nsamples - 2*self.pointing.pad) 

        elif self.config.sim_pol_type == "T":
            signal = np.convolve(0.5*self.sky_map[hit_pix], beam_kernel_row[0], mode='valid')

        elif self.config.sim_pol_type in ["QU", "_QU"]:
            signal = np.convolve(0.5*(self.sky_map[0][hit_pix]*cos2 + self.sky_map[1][hit_pix]*sin2), beam_kernel_row[1], mode='valid')
            signal += np.convolve(0.5*(-1.0*self.sky_map[0][hit_pix]*sin2 + self.sky_map[1][hit_pix]*cos2), beam_kernel_row[2], mode='valid')

        else:
            signal = np.convolve(0.5*self.sky_map[0][hit_pix], beam_kernel_row[0], mode='valid')
            signal += np.convolve(0.5*(self.sky_map[1][hit_pix]*cos2 + self.sky_map[2][hit_pix]*sin2), beam_kernel_row[1], mode='valid')
            signal += np.convolve(0.5*(-1.0*self.sky_map[1][hit_pix]*sin2 + self.sky_map[2][hit_pix]*cos2), beam_kernel_row[2], mode='valid')

        return signal

    def generate_signal_gradient(self, gradient_type, signal):
        if gradient_type == "tm_grad_co":
            gradient = (np.roll(signal["central"], -1) - np.roll(signal["central"], 1)) / (2*self.config.scan_resolution)
        elif gradient_type == "tm_grad_cross":
            gradient = (signal["top"] - signal["bottom"]) / (2*self.config.scan_resolution) 
        elif gradient_type == "tm_grad_co_co":
            gradient = (np.roll(signal["central"], -1) - 2*signal["central"] + np.roll(signal["central"], 1)) / self.config.scan_resolution**2
        elif gradient_type == "tm_grad_cross_cross":
            gradient = (signal["top"] - 2*signal["central"] + signal["bottom"]) / self.config.scan_resolution**2
        elif gradient_type == "tm_grad_co_cross":
            gradient = (np.roll(signal["top"], -1) - np.roll(signal["bottom"], -1) - np.roll(signal["top"], 1) + np.roll(signal["bottom"], 1)) / (4*self.config.scan_resolution**2)
        else:
            prompt("Gradient type : {}, does not match with the available types.\n".format(gradient_type))

        return gradient[1:-1]


    #Calculating additional scan parameters
    def setup_derived_params(self):
        #Cross scan resolution
        self.config.theta_cross = 360.0*60.0*np.sin(np.radians(self.config.alpha))*self.config.t_spin/self.config.t_prec
        #Co scan resolution
        self.config.theta_co = 360*60*np.sin(np.radians(self.config.beta))/self.config.sampling_rate/self.config.t_spin

        self.config.scan_resolution = self.config.theta_co/self.config.oversampling_rate


    def get_sky_map(self):
        if self.config.sim_pol_type == "noise_only":
            self.sky_map = None
        elif self.config.sim_pol_type == "T":
            self.sky_map = hp.read_map(self.config.input_map, verbose=False)
        elif self.config.sim_pol_type == "QU":
            self.sky_map = np.array(hp.read_map(self.config.input_map, field=(0,1), verbose=False))
        elif self.config.sim_pol_type == "_QU":
            self.sky_map = np.array(hp.read_map(self.config.input_map, field=(1,2), verbose=False))
        else:
            self.sky_map = np.array(hp.read_map(self.config.input_map, field=(0,1,2), verbose=False))

        if not self.config.sim_pol_type == "noise_only":
            map_nside = hp.get_nside(self.sky_map)
            assert map_nside == self.config.nside_in, "NSIDE of config does not match NSIDE of map"
            ## TODO
            #  if self.config.fill_empty_pixels:
                #  fill_empty_pixels(self.sky_map, 10, 0, self.config.sim_pol_type in ['QU','_QU','TQU'])

        self.npix = hp.nside2npix(self.config.nside_in)


    def set_dirs(self):
        self.sim_dir = dio.get_path_to_sim_dir(self.config.sim_tag)
        self.scan_dir = dio.get_path_to_scan_dir(self.config.sim_tag, self.config.scan_tag)
        self.detector_dir = dio.get_path_to_det_dir(self.config.sim_tag, self.config.scan_dir, self.config.detector_name)


    def make_write_dir(self, segment):
        if not os.path.exists(self.detector_dir):
            os.makedirs(self.detector_dir, exist_ok=True)

        segment_dir = dio.get_path_to_segment_dir(self.config.sim_tag, self.config.scan_tag, self.detector_name, segment)
        if not os.path.exists(segment_dir):
            os.makedirs(segment_dir, exist_ok=True)


    def write_timestream_data(self, ts_data, data_name, segment):
        write_dir = self.get_segment_dir(segment)
        if self.config.tod_type == "signal":
            if data_name == "signal":
                np.save(os.path.join(write_dir, self.config.special_tag + data_name), ts_data[::self.config.oversampling_rate])
            elif data_name == "noise":
                np.save(os.path.join(write_dir, self.config.special_tag + data_name), ts_data)
            else:
                if self.config.beam_type == "pencil":
                    np.save(os.path.join(write_dir, self.config.special_tag + data_name), ts_data)
                else:
                    np.save(os.path.join(write_dir, self.config.special_tag + data_name), ts_data[self.pointing.pad:-self.pointing.pad][::self.config.oversampling_rate])
        else:
            np.save(os.path.join(write_dir, self.config.special_tag + data_name), ts_data[::self.config.oversampling_rate])

    #*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
    #*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*

    def make_check(self):
        check_flag = 1

        # Checking main config variables
        main_var_list = ['sim_tag', 'scan_tag', 'sim_pol_type', 'coordinate_system', 'tod_type', 'oversampling_rate', 'notes', 'focal_plane_config', 'noise_type', 't_segment', 'detector_segment_dict', 'ts_data_products', 'beam_type', 'write_beam']
        main_var_exists = list(map(lambda var_name: var_name in dir(self.config), main_var_list))
        if False in main_var_exists:
            check_flag = 0 
            prompt(colored("WARNING: ", color="red") + "{} is not defined.\n".format(list(compress(main_var_list, np.logical_not(main_var_exists)))))

        sim_pol_type_list = ['IQU', 'QU', 'I', '_QU', 'noise_only']
        if self.config.sim_pol_type not in sim_pol_type_list:
            check_flag = 0
            prompt(colored("WARNING: ", color="red") + "sim_pol_ype is defines as \"{}\". It should be one of {}.\n".format(self.config.sim_pol_type, sim_pol_type_list))
        if self.config.sim_pol_type != "noise_only" and "nside_in" not in dir(self.config):
            check_flag = 0
            prompt(colored("WARNING: ", color="red") + "nside_in is not defined\n")

        if self.config.tod_type not in ['signal', 'gradient']:
            ckeck_flag = 0 
            prompt(colored("WARNING: ", color="red") + "tod_type is defined as \"{}\". It should be either \"signal\" or \"gradient\".\n".format(self.config.sim_pol_type))
        if self.config.tod_type == "gradient":
            grad_type_list = ['zero_order', 'grad_co', 'grad_cross', 'grad_coXgrad_co', 'grad_crossXgrad_cross', 'grad_coXgrad_cross']
            if not set(self.config.grad_type).issubset(set(grad_type_list)):
                check_flag = 0
                prompt(colored("WARNING: ", color="red") + "grad_type is set to {}. Please check as it is not one of the defined types\n".format(self.config.grad_type))

        if self.config.coordinate_system not in ['ecliptic', 'galactic']:
            check_flag = 0
            prompt(colored("WARNING: ", color="red") + "coordinate_system is defined as \"{}\". It should be either \"ecliptic\" or \"galactic\".\n".format(self.config.coordinate_system))
        #-----------------------------------------------------------

        # Checking that the scan strategy variables are in place

        scan_strategy_var_list = ['t_year', 't_prec', 't_spin', 'sampling_rate', 'alpha', 'beta', 'scan_strategy_name', 'scan_strategy_note']
        scan_strategy_var_exists = list(map(lambda var_name: var_name in dir(self.config), scan_strategy_var_list))
        if False in scan_strategy_var_exists:
            check_flag = 0
            prompt(colored("WARNING: ", color="red") + "{} is not defined for the scan strategy.\n".format(list(compress(scan_strategy_var_list, np.logical_not(scan_strategy_var_exists)))))
        #-----------------------------------------------------------

        # Checking that the I/O variables are in place
        ts_data_product_list = ['signal', 'pointing_vec', 'pol_ang', 'noise', 'mask']
        if not set(self.config.ts_data_products).issubset(set(ts_data_product_list)):
            check_flag = 0
            prompt(colored("WARNING: ", color="red") + "ts_data_products are {}. Please check as it is not one of the defined types\n".format(self.config.ts_data_products))
        #-----------------------------------------------------------

        # Noise and beam settings
        beam_type_list = ['pencil', 'full_simulated', 'from_file']
        if not self.config.beam_type in beam_type_list:
            check_flag = 0
            prompt(colored("WARNING: ", color="red") + "beam_type is {}. Please check as it is not one of the defined types\n".format(self.config.beam_type))
        
        if self.config.beam_type == "full_simulated":
            if "beam_cutoff" not in dir(self.config):
                check_flag = 0
                prompt(colored("WARNING: ", color="red") + "beam_cutoff not defined")

        noise_type_list = ['white', '1_over_f', 'none']
        if self.config.noise_type not in noise_type_list:
            check_flag = 0
            prompt(colored("WARNING: ", color="red") + "noise_type is {}. Please check as it is not one of the defined types\n".format(self.config.noise_type))
        #-----------------------------------------------------------

        return check_flag


    def display_params(self):
        display_string = colored("#* Detector {} parameters\n".format(self.config.detector_name), color="green")
        display_string += "Alpha: {} degrees\n".format(self.config.alpha)
        display_string += "Beta: {} degrees\n".format(self.config.beta)
        display_string += "Sampling rate: {} Hz\n".format(self.config.sampling_rate)
        display_string += "Theta co: {} arcmin\n".format(self.config.theta_co)
        display_string += "Theta cross: {} arcmin\n".format(self.config.theta_cross)
        display_string += "Oversampling rate: {}\n".format(self.config.oversampling_rate)
        display_string += "Scan resolution for beam integration: {} arcmin\n".format(self.config.scan_resolution)
        if self.config.sim_pol_type != "noise_only":
            display_string += "Pixel size for NSIDE = %d: {} arcmin\n".format(self.config.nside_in, hp.nside2resol(self.config.nside_in, arcmin=True))
        n_steps = int(self.config.t_segment*self.config.sampling_rate)*self.config.oversampling_rate 
        display_string += "#Samples per segment: {}\n".format(n_steps)
        prompt(display_string, sys.stdout)
