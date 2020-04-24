import numpy as np
import healpy as hp
from termcolor import colored
from genesys import Genesys_Class
from genesys.pointing import Pointing
from genesys.noise import Noise
 
class Detector(Genesys_Class):
    def __init__(self, channel_obj, detector_name):
        self.params = {}
        self.params['detector_name'] = detector_name
        self.params['channel_name'] = channel_obj.params['channel_name']
        self.params.update(channel_obj.params['detectors'][detector_name])
        self.params['scan_strategy'] = {}
        self.params['scan_strategy'].update(channel_obj.params['scan_strategy'])
        if 'noise' not in self.params:
            self.params['noise'] = {}
        get_param_if_None(self.params['noise'], channel_obj.params['noise'], ['noise_type', 'white_noise_rms']) 
        if self.params['noise']['noise_type'] == '1_over_f':
            get_param_if_None(self.params['noise'], channel_obj.params['noise'], ['f_knee', 'noise_alpha'])
        # Other params
        get_param_if_None(self.params, channel_obj.params, ['sampling_rate', 'pol_modulation', 'input_map_file'])
        if self.params['pol_modulation'] != 'scan':
            self.params['HWP'] = {}
            self.params['HWP'].update(channel_obj.params['HWP'])

    #*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
    # This is where all the action happens
    # Scanning the map with the simulated pointing
    #*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
        
    def observe_sky(self, pointing, sky_map, tod, pol_type):
        axis_detector_sight = pointing.get_detector_axis(self.params['pos'])
        theta, phi, psi = pointing.get_pointing_and_pol(axis_detector_sight)

        if self.params['pol_modulation'] != 'scan':
            pol_phase = 2*psi + 4*tod['hwp_psi']
        else:
            pol_phase = 2*psi
        # TO DO: orbital_dipole = self.pointing.get_orbital_dipole()

        hit_pix = hp.ang2pix(sky_map.nside, theta, phi)
        if pol_type == 'noise':
            signal = np.zeros(tod['t_steps'].size)
        if pol_type == 'I':
            signal = sky_map.sky_map[hit_pix]
        if pol_type == 'QU':
            signal = sky_map.sky_map[0][hit_pix]*np.cos(pol_phase) + sky_map.sky_map[1][hit_pix]*np.sin(pol_phase)
        if pol_type == '_QU':
            signal = sky_map.sky_map[1][hit_pix]*np.cos(pol_phase) + sky_map.sky_map[2][hit_pix]*np.sin(pol_phase)
        if pol_type == 'IQU':
            signal = sky_map.sky_map[0][hit_pix] + sky_map.sky_map[1][hit_pix]*np.cos(pol_phase) + sky_map.sky_map[2][hit_pix]*np.sin(pol_phase)

        noise_params = self.params['noise']
        noise_params['sampling_rate'] = self.params['sampling_rate']
        noise_obj = Noise(noise_params)
        noise = noise_obj.simulate_noise(tod['t_steps'].size)
        signal += noise

        tod['theta'] = theta
        tod['phi'] = phi
        tod['psi'] = psi
        tod['signal'] = signal
        tod['noise'] = noise
                
    #*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
    # Paramtere display routines
    #*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*

    def display_hwp(self):
        hwp_param = self.params['HWP']
        self.prompt("HALF WAVE PLATE:")
        for item in list(hwp_param.keys()):
            self.prompt(f"\t{item:<35}{hwp_param[item]}")

    def display_noise(self):
        noise_param = self.params['noise']
        self.prompt("noise:")
        for item in list(noise_param.keys()):
            self.prompt(f"\t{item:<27}{noise_param[item]}")

    def display_scan_strategy(self):
        unit_dict = {'alpha': 'degrees', 'beta': 'degrees', 't_precession': 'seconds', 't_spin': 'seconds', 'duration': 'years'}
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
