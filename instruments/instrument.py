"""
class for handling instrument parameter loading and handling
The instrument is specified by a instrument_name folder.
The folder contains the following files
    instrument.ini -> Top level instrument specifications, scan strategies, polarisation modulator
    focal_plane.ini -> Focal plane and channel descriptions.
    detectors_<channel>.ini -> Individual detector configuration for each channel
If None is provided for focal_plane and scan_strategy, the default files will be loaded. Otherwise the custom files will be loaded
"""

import os
from .. import Genesys_Class

class Instrument(Genesys_Class):
    def __init__(self, other=None, instrument_name=None):
        """
        Constructor for Instrument class
        If other is given, the instance variables of other are copied into self
        If instrument_name is given, the intrument parameters are loaded from the instrument_dir
        """
        if other != None:
            self.copy_other(other)
        elif instrument_name != None:
            self.configure_instrument(instrument_name=instrument_name)
        else:
            pass

    def configure_instrument(self, instrument_name):
        """
        Reads in the instrument parameters from the files in the intstrument_folder
        """
        instrument_param_file_path = self.path_to_instrument_param_file(instrument_name)
        self.params = self.load_param_file(file_path=instrument_param_file_path)
        if 'half_wave_plates' in self.params.keys():
            self.hwp_list = list(self.params['half_wave_plates'].keys())
        self.channel_list = list(self.params['channels'].keys())

    #*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
    # Paramter display methods
    #*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*

    def display_scan_strategy(self):
        unit_dict = {'alpha': 'degrees', 'beta': 'degrees', 't_precession': 'seconds', 't_spin': 'seconds', 'duration': 'years'}
        scan_strategy = self.params['scan_strategy']
        for item in list(scan_strategy.keys()):
            print("{}: {} {}".format(item, scan_strategy[item], unit_dict[item]))

    def display_half_wave_plate(self):
        hwps = self.params['half_wave_plates']
        for hwp in list(hwps.keys()):
            print("{}:".format(hwp))
            for item in list(hwps[hwp].keys()):
                print("\t{}: {}".format(item, hwps[hwp][item]))

    def display_channel(self, channel_list=None):
        if channel_list == None:
            channel_list = self.channel_list
        for channel in channel_list:
            print("{}:".format(channel))
            for item in list(self.params['channels'][channel]):
                print("\t{}: {}".format(item, self.params['channels'][channel][item]))

    #*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
    # Path naming conventions
    #*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*

    def path_to_instrument_dir(self, instrument_dir_name):
        return os.path.join(self.global_paths['instruments_dir'], instrument_dir_name)

    def path_to_instrument_param_file(self, instrument_dir_name, instrument_param_file_name="instrument_params.yaml"):
        """
        Default instrument param file name: instrument.yaml
        """
        return os.path.join(self.path_to_instrument_dir(instrument_dir_name), instrument_param_file_name)

    #  #*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
    #  # Validating the parameters
    #  #*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
#
    #  def set_main_instrument_param_fields(self, main_instrument_param_fileds=None):
        #  """
        #  Sets the main parameter fields for the instrument
        #  """
        #  if main_instrument_param_fileds == None:
            #  self.main_instrument_param_fields = ['instrument_name', 'version', 'scan_strategy', 'half_wave_plates', 'channels']
        #  else:
            #  self.main_instrument_param_fields = main_instrument_param_fields
#
    #  def set_channel_param_fields(self, channel_param_fields=None):
        #  """
        #  Sets the parameter fields for the channels
        #  """
        #  if channel_param_fileds == None:
            #  self.channel_param_fields = ['central_frequency', 'bandwidth', 'beam_fwhm', 'pixel_size', 'num_detectors', 'detector_NET', 'detector_yield', 'noise_margin', 'pol_modulation', 'channel_param_file']
        #  else:
            #  self.channel_param_fields = channel_param_fields
#
    #  def set_scan_strategy_param_fields(self, scan_strategy_param_fields=None):
        #  """
        #  Sets the parameter fields for the scan strategy
        #  """
        #  if scan_strategy_param_fields == None:
            #  self.scan_strategy_param_fields = ['alpha', 'beta', 't_precession', 't_spin', 'duration', 'data_rate']
        #  else:
            #  self.scan_strategy_param_fields = scan_strategy_param_fields
#
    #  def set_hwp_param_fields(self, hwp_param_fields=None)
        #  """
        #  Sets the parameter fields for the half wave plate
        #  """
        #  if hwp_param_fields == None:
            #  self.hwp_param_fields = ['rotation', 'rpm']
        #  else:
            #  self.hwp_param_fields = hwp_param_fields
#
    #  def set_main_instrument_param_to_default(self):
        #  for
