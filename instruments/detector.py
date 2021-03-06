import numpy as np
import healpy as hp
import copy
from termcolor import colored
from genesys import Genesys_Class
from genesys.pointing import Pointing
from genesys.noise import Noise
from .default_units import unit_dict
 
class Detector(Genesys_Class):
    def __init__(self, channel_obj, detector_name):
        self.params = {}
        self.params['detector_name'] = detector_name
        self.params['channel_name'] = channel_obj.params['channel_name']
        self.params.update(channel_obj.params['detectors'][detector_name])
        self.params['scan_strategy'] = {}
        self.params['scan_strategy'].update(channel_obj.params['scan_strategy'])
        self.params['HWP'] = {}
        self.params['HWP'].update(channel_obj.params['HWP'])
        get_param_if_None(self.params['noise'], channel_obj.params['noise'], ['white_noise_rms', 'f_knee', 'noise_alpha']) 
        get_param_if_None(self.params, channel_obj.params, ['sampling_rate', 'pol_modulation', 'central_frequency'])

    #*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
    # This is where all the action happens
    # Scanning the map with the simulated pointing
    #*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
        
    def observe_sky(self, pnt_obj, sky_map_obj, tod, segment_length, sat_vel, noise_type, tod_write_field):
        sky_map = sky_map_obj.sky_map

        self.axis_sight = pnt_obj.get_detector_axis(self.params['pos'])
        v_pointing = pnt_obj.get_pointing(self.axis_sight)
        psi = pnt_obj.get_polariser_phase(v_pointing, self.params['pol_phase_ini'])

        if np.sum(sat_vel**2) != 0.0:
            orb_dip = pnt_obj.get_orbital_dipole(v_pointing, sat_vel, self.params['central_frequency']) 

        theta, phi = hp.vec2ang(v_pointing)
        del v_pointing

        if self.params['pol_modulation'] == 'passive':
            pol_phase = 2*psi
        else:
            pol_phase = 2*psi + 4*tod['psi_hwp']

        hit_pix = hp.ang2pix(sky_map_obj.nside, theta, phi)

        signal = sky_map[0][hit_pix] + sky_map[1][hit_pix]*np.cos(pol_phase) + sky_map[2][hit_pix]*np.sin(pol_phase)

        if np.sum(sat_vel**2) != 0.0:
            signal += orb_dip

        if noise_type != 'no_noise':
            noise_params = copy.deepcopy(self.params['noise'])
            noise_params.update({'sampling_rate': self.params['sampling_rate'], 'noise_type': noise_type})
            noise_obj = Noise(noise_params)
            noise = noise_obj.simulate_noise(segment_length)
            signal += noise
            if "noise" in tod_write_field:
                tod['noise'] = noise

        tod['theta'] = theta
        tod['phi'] = phi
        tod['psi'] = psi
        tod['signal'] = signal
        if "orb_dip" in tod_write_field:
            tod['orb_dip'] = orb_dip
                
    #*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
    # Paramtere display routines
    #*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*

    def display_hwp(self):
        hwp_param = self.params['HWP']
        self.prompt("HALF WAVE PLATE:")
        for item in list(hwp_param.keys()):
            self.prompt(f"\t{item:<27}{hwp_param[item]} {unit_dict[item]}")

    def display_noise(self):
        noise_param = self.params['noise']
        self.prompt("noise:")
        for item in list(noise_param.keys()):
            self.prompt(f"\t{item:<27}{noise_param[item]} {unit_dict[item]}")

    def display_scan_strategy(self):
        ss_params = self.params['scan_strategy']
        self.prompt("SCAN STRATEGY:")
        for item in ss_params.keys():
            self.prompt(f"\t{item:<27}{ss_params[item]} {unit_dict[item]}")

    def info(self):
        self.prompt(colored(f"#{self.params['channel_name'] + ' -- ' + self.params['detector_name']:^29}#", 'green')) 
        for item in self.params.keys():
            if item == 'noise':
                self.display_noise()
            elif item == 'HWP' and self.params['pol_modulation'] != 'scan':
                self.display_hwp()
            elif item in ['detector_name', 'channel_name', 'scan_strategy']:
                pass
            else:
                self.prompt(f"{item:<35}{self.params[item]}")
        self.prompt(colored(15 * '#*' + '#', color="green"))

                
#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
# Helper routines not part of detector class
#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*

def get_param_if_None(dict_1, dict_2, items):
    for item in items:
        if item not in dict_1 or dict_1[item] == None:
            dict_1[item] = dict_2[item]
