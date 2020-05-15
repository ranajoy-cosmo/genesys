import os
from termcolor import colored
from genesys import Genesys_Class
from genesys.instruments.detector import Detector
from .default_units import unit_dict

class Channel(Genesys_Class):
    """
    Class containing the parameters of a particular frequency channel
    """
    def __init__(self, inst_obj=None, channel_name=None):
        """
        Configure the channel params.
        The channel object can only be created from the instrument object
        """
        if inst_obj != None:
            self.params = {}
            self.params['channel_name'] = channel_name
            self.params.update(inst_obj.params['channels'][channel_name])
            # HWP and scan strategy
            pol_mod = self.params['pol_modulation']
            if pol_mod == 'passive':
                self.params['HWP'] = {'HWP_name': None, 'spin_rate': 0}
            else:
                self.params['HWP'] = {}
                self.params['HWP']['HWP_name'] = pol_mod
                self.params['HWP'].update(inst_obj.params['HWP'][pol_mod])
            self.params['scan_strategy'] = {}
            self.params['scan_strategy'].update(inst_obj.params['scan_strategy'])
            # All detector parameters from channel file
            channel_param_file_path = os.path.join(inst_obj.instrument_dir, inst_obj.params['channels'][channel_name]['channel_param_file'])
            self.params['detectors'] = self.load_param_file(file_path=channel_param_file_path)
        else:
            self.params = {}

    #*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
    # Paramtere display routines
    #*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*

    def info(self):
        self.prompt(colored(f"#{self.params['channel_name']:^29}#", 'green')) 
        for item in self.params:
            if item == 'noise':
                self.prompt("\tnoise/detector")
                for noise_item in list(self.params['noise']):
                    self.prompt(f"\t\t{noise_item:<19}{self.params['noise'][noise_item]} {unit_dict[noise_item]}")
            elif item == 'detectors':
                self.prompt("DETECTORS: {list(self.params['detectors'].keys())}")
            else:
                self.prompt(f"\t{item:<27}{self.params[item]} {unit_dict[item]}")
        prompt(colored("#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*\n", color="green"))

    #*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
    # DETECTOR OBJECT 
    #*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
    
    def get_detector(self, detector_name):
        detector_obj = Detector(self, detector_name)
        return detector_obj
