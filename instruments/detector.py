import os
import copy
from genesys.instruments.instrument import Instrument
from genesys.pointing import Pointing
from genesys.maps import Sky_Map
 
class Detector(Instrument):
    # Instrument ALREADY INHERITS FROM Genesys_Class
    def __init__(self, instrument_name=None, channel_name=None, detector_name=None, other=None):
        if detector_name != None:
            self.configure_detector(instrument_name, channel_name, detector_name)
        elif other != None:
            self.copy_attributes(other)
        else:
            pass

    def configure_detector(self, instrument_name, channel_name, detector_name):
        # GETTING THE INSTRUMENT PARAMS
        instrument_param_file_path = self.path_to_instrument_param_file(instrument_name)
        inst_params = self.load_param_file(file_path=instrument_param_file_path)
        channel_params = inst_params['channels'][channel_name]
        # GETTING THE DETECTOR PARAMS FOR THE ENTIRE CHANNEL
        detector_file_name = channel_params['detector_param_file']
        detector_param_file_path = self.path_to_detector_param_file(instrument_name, detector_file_name)
        detector_params = self.load_param_file(file_path=detector_param_file_path)
        # GETTING THE DETECTOR SPECIFIC PARAMS
        self.params = detector_params[detector_name]
        self.params['scan_strategy'] = copy.deepcopy(inst_params['scan_strategy'])
        self.params['channel_name'] = channel_name
        # ADDING CHANNEL-WIDE PARAMETERS IF NOT PRESENT FOR DETECTOR
        if 'noise_type' not in self.params['noise'] or self.params['noise']['noise_type'] is None:
            self.params['noise']['noise_type'] = channel_params['noise_type']
        if 'white_noise_rms' not in self.params['noise'] or self.params['noise']['white_noise_rms'] is None:
            self.params['noise']['white_noise_rms'] = channel_params['detector_NET']
        if self.params['noise']['noise_type'] == '1_over_f':
            if 'f_knee' not in self.params['noise'] or self.params['noise']['f_knee'] is None:
                self.params['noise']['f_knee'] = channel_params['f_knee']
            if 'noise_alpha' not in self.params['noise'] or self.params['noise']['noise_alpha'] is None:
                self.params['noise']['noise_alpha'] = channel_params['noise_alpha']
        if 'sampling_rate' not in self.params or self.params['sampling_rate'] is None:
            self.params['sampling_rate'] = channel_params['sampling_rate']
        if 'pol_modulation' not in self.params or self.params['pol_modulation'] is None:
            self.params['pol_modulation'] = channel_params['pol_modulation']
        if self.params['pol_modulation'] != 'scan':
            self.params['HWP'] = inst_params['half_wave_plates'][channel_params['pol_modulation']]
            self.params['HWP']['HWP_label'] = channel_params['pol_modulation']
        if 'offset' not in self.params or self.params['offset'] == None:
            self.params['offset'] = [0.0,0.0]

    def initialise_pointing(self):
        self.pointing = Pointing(self.params)

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

    def display_half_wave_plate(self):
        hwp = self.params['HWP']
        print("HALF WAVE PLATE:")
        for item in list(hwp.keys()):
            print("\t\t{}: {}".format(item, hwp[item]))

    def info(self):
        print("Channel {} and Detector {}".format(self.params['channel_name'], self.params['detector_name'])) 
    #*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
    # Path naming conventions
    #*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*

    def path_to_detector_param_file(self, instrument_name, detector_file_name):
        return os.path.join(self.path_to_instrument_dir(instrument_name), detector_file_name)
