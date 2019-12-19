import os
import copy
from genesys import Genesys_Class
from genesys.instruments.detector import Detector

class Channel(Genesys_Class):
    """
    Class containing the specifications of a particular frequency channel
    """
    def __init__(self, inst_obj, channel_name):
        """
        Configure the channel params
        """
        self.params = {}
        self.params['commons'] = {}
        self.params['commons']['name'] = channel_name
        # IDEAL CHANNEL _PARAMETERS
        self.params['commons'].update(inst_obj.params['channels'][channel_name])
        # HWP AND SCAN STRATEGY
        pol_mod_name = self.params['commons']['pol_modulation']
        if pol_mod_name != 'scan':
            self.params['commons']['HWP'] = copy.deepcopy(inst_obj.params['HWP'][pol_mod_name])
            self.params['commons']['HWP']['HWP_label'] = pol_mod_name
        self.params['commons']['scan_strategy'] = copy.deepcopy(inst_obj.params['scan_strategy'])
        # ALL DETECTOR PARAMETERS FROM CHANNEL FILE
        channel_param_file_name = self.params['commons']['channel_param_file']
        channel_param_file = os.path.join(inst_obj.path_to_instrument_dir(), channel_param_file_name)
        channel_params = self.load_param_file(file_path=channel_param_file)
        self.params['detectors'] = channel_params
        self.detector_list = list(channel_params.keys())

    #*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
    # PARAMTERE DISPLAY ROUTINES
    #*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
    def display_channel_commons(self):
        unit_dict = {'name': '', 'central_frequency': 'GHz', 'bandwidth': 'GHz', 'beam_fwhm': 'arcmins', 'sampling_rate': 'Hz', 'pixel_size': 'mm', 'num_detectors': '', 'noise_type': '', 'white_noise_rms': 'uk.sqrt(s)', 'detector_yield': '', 'noise_margin': '', 'pol_modulation': '', 'channel_param_file': '', 'rotation': '', 'rpm': '', 'HWP_label': ''}
        unit_dict.update({'alpha': 'degrees', 'beta': 'degrees', 't_precession': 'seconds', 't_spin': 'seconds', 'duration': 'years'})
        print("IDEAL CHANNEL PARAMETERS")
        for item in list(self.params['commons']):
            if item == 'noise':
                print("\tnoise/detector:")
                for noise_item in list(self.params['commons']['noise']):
                    print(f"\t\t{noise_item}: {self.params['commons']['noise'][noise_item]} {unit_dict[noise_item]}")
            elif item == 'HWP':
                print("\tHWP:")
                for hwp_item in list(self.params['commons']['HWP']):
                    print(f"\t\t{hwp_item}: {self.params['commons']['HWP'][hwp_item]} {unit_dict[hwp_item]}")
            elif item == 'scan_strategy':
                print("\tscan_strategy:")
                for ss_item in list(self.params['commons']['scan_strategy']):
                    print(f"\t\t{ss_item}: {self.params['commons']['scan_strategy'][ss_item]} {unit_dict[ss_item]}")
            else:
                print(f"\t{item}: {self.params['commons'][item]} {unit_dict[item]}")

    def display_detectors(self):
        unit_dict = {'name': '', 'input_map_path': '', 'pol_phase_ini': 'degrees', 'central_frequency': 'GHz', 'bandwidth': 'GHz', 'beam_fwhm': 'arcmins', 'sampling_rate': 'Hz', 'pixel_size': 'mm', 'num_detectors': '', 'noise_type': '', 'white_noise_rms': 'uk.sqrt(s)', 'f_knee': 'mHz', 'noise_alpha': '', 'detector_yield': '', 'noise_margin': '', 'pol_modulation': '', 'channel_param_file': '', 'pos': 'arcmins'}
        print("DETECTORS")
        for detector in self.detector_list:
            print(f"\t{detector}:")
            for item in list(self.params['detectors'][detector]):
                if item != 'noise':
                    print(f"\t\t{item}: {self.params['detectors'][detector][item]} {unit_dict[item]}")
                else:
                    print("\t\tnoise:")
                    for noise_item in list(self.params['detectors'][detector]['noise']):
                        print(f"\t\t\t{noise_item}: {self.params['detectors'][detector]['noise'][noise_item]} {unit_dict[noise_item]}")

    def info(self):
        self.display_channel_commons()
        print('\n')
        self.display_detectors()

    #*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
    # DETECTOR OBJECT 
    #*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
    
    def get_detector_object(self, detector_name):
        pass
