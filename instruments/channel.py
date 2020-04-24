import os
from termcolor import colored
from genesys import Genesys_Class
from genesys.instruments.detector import Detector

class Channel(Genesys_Class):
    """
    Class containing the parameters of a particular frequency channel
    """
    def __init__(self, inst_obj, channel_name):
        """
        Configure the channel params.
        The channel object can only be created from the instrument object
        """
        self.params = {}
        self.params['channel_name'] = channel_name
        self.params.update(inst_obj.params['channels'][channel_name])
        # HWP and scan strategy
        pol_mod_name = self.params['pol_modulation']
        if pol_mod_name != 'scan':
            self.params['HWP'] = {}
            self.params['HWP'].update(inst_obj.params['HWP'][pol_mod_name])
            self.params['HWP']['HWP_name'] = pol_mod_name
        self.params['scan_strategy'] = {}
        self.params['scan_strategy'].update(inst_obj.params['scan_strategy'])
        # All detector parameters from channel file
        channel_param_file_path = os.path.join(inst_obj.instrument_dir, inst_obj.params['channels'][channel_name]['channel_param_file'])
        self.params['detectors'] = self.load_param_file(file_path=channel_param_file_path)

    #*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
    # Paramtere display routines
    #*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*

    def info(self):
        unit_dict = {'central_frequency': 'GHz', 'bandwidth': 'GHz', 'beam_fwhm': 'arcmins', 'sampling_rate': 'Hz', 'pixel_size': 'mm', 'num_detectors': '', 'noise_type': '', 'white_noise_rms': 'uk.sqrt(s)', 'detector_yield': '', 'noise_margin': '', 'pol_modulation': '', 'channel_param_file': '', 'input_map_file': '', 'f_knee': 'mHz', 'noise_alpha': ''}
        prompt(colored("#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*\n", color="green"))
        self.prompt(f"CHANNEL: {self.params[channel_name]}:")
        for item in self.params:
            if item == 'noise':
                self.prompt("\tnoise/detector")
                for noise_item in list(self.params['noise']):
                    self.prompt(f"\t\t{noise_item}: {self.params['noise'][noise_item]} {unit_dict[noise_item]}")
            elif item == 'detectors':
                self.prompt("DETECTORS: {list(self.params['detectors'].keys())}")
            else:
                self.prompt(f"\t{item}: {self.params[item]} {unit_dict[item]}")
        prompt(colored("#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*\n", color="green"))

    #*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
    # DETECTOR OBJECT 
    #*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
    
    def get_detector(self, detector_name):
        detector_obj = Detector(self, detector_name)
        return detector_obj
