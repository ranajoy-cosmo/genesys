import numpy as np
import healpy as hp
from termcolor import colored
from genesys import Genesys_Class
from genesys.pointing import Pointing
from genesys.maps import Sky_Map
 
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
        # OTHER PARAMS
        get_param_if_None(self.params, channel_obj.params, ['sampling_rate', 'pol_modulation', 'input_map_file'])
        if self.params['pol_modulation'] != 'scan':
            self.params['HWP'] = {}
            self.params['HWP'].update(channel_obj.params['HWP'])

    def load_map(self, pol_type='IQU'):
        """
        pol_type:
            I -> Only I read. field=(0)
            IQU -> All I,Q,U read. field=(0,1,2)
            _QU -> All I,Q,U read but I set to 0. field=(0,1,2) 
        """
        input_map_path = self.params['input_map_file']
        if pol_type == 'I':
            self.sky_map = Sky_Map(self.params['input_map_file'], field=(0))
        if pol_type == 'IQU':
            self.sky_map = Sky_Map(self.params['input_map_file'], field=(0,1,2))
        if pol_type == '_QU':
            self.sky_map = Sky_Map(self.params['input_map_file'], field=(0,1,2))
            self.sky_map[0] *= 0.0

    #*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
    # THIS IS WHERE ALL THE ACTION HAPPENS
    # SCANNING THE MAP WITH THE SIMULATED POINTING
    #*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*

    def prepare_for_observation(self, t_stream, segment_length, segment):
        t_stream.gen_t_steps(segment_length, segment, self.params['sampling_rate'])
        self.initialise_pointing()

    def initialise_pointing(self):
        pointing_params = {}
        pointing_params.update(self.params['scan_strategy'])
        pointing_params['pos'] = self.params['pos']
        pointing_params['HWP_spin_rate'] = self.params['HWP']['spin_rate']
        self.pointing = Pointing(pointing_params)
        
    def observe_sky(self, t_stream, coordinate_system, pol_type):
        theta, phi, psi = self.pointing.get_pointing_and_pol_angle(t_stream.t_steps, coordinate_system)
        if self.params['pol_modulation'] != 'scan':
            hwp_psi = self.pointing.get_hwp_angle(t_stream.t_steps)
            pol_phase = 2*psi + 4*hwp_psi
        else:
            pol_phase = 2*psi
        # TO DO: orbital_dipole = self.pointing.get_orbital_dipole()
        hit_pix = hp.ang2pix(self.sky_map.nside, theta, phi)
        if pol_type == 'noise':
            signal = np.zeros(t_stream.t_steps.size)
        if pol_type == 'I':
            signal = self.sky_map.sky_map[hit_pix]
        if pol_type == 'QU':
            signal = self.sky_map.sky_map[0][hit_pix]*np.cos(pol_phase) + self.sky_map.sky_map[1][hit_pix]*np.sin(pol_phase)
        if pol_type == '_QU':
            signal = self.sky_map.sky_map[1][hit_pix]*np.cos(pol_phase) + self.sky_map.sky_map[2][hit_pix]*np.sin(pol_phase)
        if pol_type == 'IQU':
            signal = self.sky_map.sky_map[0][hit_pix] + self.sky_map.sky_map[1][hit_pix]*np.cos(pol_phase) + self.sky_map.sky_map[2][hit_pix]*np.sin(pol_phase)

        t_stream.theta = theta
        t_stream.phi = phi
        t_stream.psi = psi
        t_stream.hwp_psi = hwp_psi
        t_stream.signal = signal
                
    #*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
    # PARAMTERE DISPLAY ROUTINES
    #*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*

    def display_hwp(self):
        hwp_param = self.params['HWP']
        self.prompt("HALF WAVE PLATE:")
        for item in list(hwp_param.keys()):
            self.prompt(f"\t{item}: {hwp_param[item]}")

    def display_noise(self):
        noise_param = self.params['noise']
        self.prompt("Noise:")
        for item in list(noise_param.keys()):
            self.prompt(f"\t{item}: {noise_param[item]}")

    def display_scan_strategy(self):
        unit_dict = {'alpha': 'degrees', 'beta': 'degrees', 't_precession': 'seconds', 't_spin': 'seconds', 'duration': 'years'}
        ss_params = self.params['scan_strategy']
        self.prompt("SCAN STRATEGY:")
        for item in ss_params.keys():
            self.prompt(f"\t{item}: {ss_params[item]} {unit_dict[item]}")

    def info(self):
        self.prompt(colored("#*#*#* ", color="green") + f"{self.params['channel_name']} -- {self.params['detector_name']}"+ colored(" #*#*#* ", color="green")) 
        for item in self.params.keys():
            if item == 'noise':
                self.display_noise()
            elif item == 'HWP' and self.params['pol_modulation'] != 'scan':
                self.display_hwp()
            elif item in ['detector_name', 'channel_name', 'scan_strategy']:
                pass
            else:
                self.prompt(f"{item}: {self.params[item]}")
        self.prompt(colored("#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*", color="green"))

                
#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
# HELPER ROUTINES NOT PART OF Detector CLASS
#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*

def get_param_if_None(dict_1, dict_2, items):
    for item in items:
        if item not in dict_1 or dict_1[item] == None:
            dict_1[item] = dict_2[item]
