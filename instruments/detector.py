import os
import copy
from genesys import Genesys_Class
from genesys.pointing import Pointing
#  from genesys.maps import Sky_Map
#  from genesys.timestream.timestream import TStream
 
class Detector(Genesys_Class):
    # Instrument ALREADY INHERITS FROM Genesys_Class
    def __init__(self, channel_obj, detector_name, other=None):
        if other != None:
            self.copy_attributes(other)
        else:
            self.params = {}
            self.params.update(channel_obj.params['detectors'][detector_name])
            self.params['name'] = detector_name
            self.params['channel_name'] = channel_obj.params['commons']['name']
            self.params['scan_strategy'] = copy.deepcopy(channel_obj.params['commons']['scan_strategy'])
            if 'noise' not in self.params:
                self.params['noise'] = {}
            get_param_if_None(self.params['noise'], channel_obj.params['commons']['noise'], ['noise_type', 'white_noise_rms']) 
            if self.params['noise']['noise_type'] == '1_over_f':
                get_param_if_None(self.params['noise'], channel_obj.params['commons']['noise'], ['f_knee', 'noise_alpha'])
            # OTHER PARAMS
            get_param_if_None(self.params, channel_obj.params['commons'], ['sampling_rate', 'pol_modulation'])
            if self.params['pol_modulation'] != 'scan':
                self.params['HWP'] = channel_obj.params['commons']['HWP']

    def initialise_pointing(self):
        pointing_params = {}
        pointing_params.update(self.params.scan_strategy)
        pointing_params['pos'] = self.params['pos']
        self.pointing = Pointing(self.params)

    def initialise_timestream_segment(self, sim_config, segment):
        self.ts = TStream(sim_config, segment)

    def load_map(self, pol_type='IQU'):
        """
        pol_type:
            I -> Only I read. field=(0)
            IQU -> All I,Q,U read. field=(0,1,2)
            _QU -> All I,Q,U read but I set to 0. field=(0,1,2) 
        """
        input_map_path = self.params['input_map_path']
        if pol_type == 'I':
            self.sky_map = Sky_Map(self.params['input_map_path'], field=(0))
        if pol_type == 'IQU':
            self.sky_map = Sky_Map(self.params['input_map_path'], field=(0,1,2))
        if pol_type == '_QU':
            self.sky_map = Sky_Map(self.params['input_map_path'], field=(0,1,2))
            self.sky_map[0] *= 0.0
                
    #*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
    # PARAMTERE DISPLAY ROUTINES
    #*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*

    def display_hwp(self):
        hwp_param = self.params['HWP']
        print("HALF WAVE PLATE:")
        for item in list(hwp_param.keys()):
            print(f"\t{item}: {hwp_param[item]}")

    def display_noise(self):
        noise_param = self.params['noise']
        print("Noise:")
        for item in list(noise_param.keys()):
            print(f"\t{item}: {noise_param[item]}")

    def display_scan_strategy(self):
        unit_dict = {'alpha': 'degrees', 'beta': 'degrees', 't_precession': 'seconds', 't_spin': 'seconds', 'duration': 'years'}
        ss_params = self.params['scan_strategy']
        print("SCAN STRATEGY:")
        for item in ss_params.keys():
            print(f"\t{item}: {ss_params[item]} {unit_dict[item]}")

    def info(self):
        print(f"Channel {self.params['channel_name']} and Detector {self.params['name']}") 
        for item in self.params.keys():
            if item == 'noise':
                self.display_noise()
            elif item == 'HWP' and self.params['pol_modulation'] != 'scan':
                self.display_hwp()
            elif item == 'scan_strategy':
                self.display_scan_strategy()
            else:
                print(f"{item}: {self.params[item]}")
    #*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
    # Path naming conventions
    #*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*

    def path_to_detector_param_file(self, instrument_name, detector_file_name):
        return os.path.join(self.path_to_instrument_dir(instrument_name), detector_file_name)

def get_param_if_None(dict_1, dict_2, items):
    for item in items:
        if item not in dict_1 or dict_1[item] == None:
            dict_1[item] = dict_2[item]
