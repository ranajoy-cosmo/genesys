import os
from genesys import Genesys_Class

class Instrument(Genesys_Class):
    """
    class for handling instrument parameter loading and handling
    The instrument is described in a instrument_name folder.
    The folder contains the following files
        instrument.yaml -> Top level instrument specifications, scan strategies, polarisation modulator
        <Channel_name>.yaml -> Individual detector parameters.
    """
    def __init__(self, instrument_name=None, other=None):
        """
        Constructor for Instrument class
        If other is given, the instance variables of other are copied into self
        If instrument_name is given, the intrument parameters are loaded from the corresponding directory
        """
        if instrument_name != None:
            self.configure_instrument(instrument_name=instrument_name)
        elif other != None:
            self.copy_attributes(other)
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
        print("SCAN STRATEGY:")
        for item in list(scan_strategy.keys()):
            print("\t{}: {} {}".format(item, scan_strategy[item], unit_dict[item]))

    def display_half_wave_plate(self):
        hwps = self.params['half_wave_plates']
        print("HALF WAVE PLATES:")
        for hwp in list(hwps.keys()):
            print("\t{}:".format(hwp))
            for item in list(hwps[hwp].keys()):
                print("\t\t{}: {}".format(item, hwps[hwp][item]))

    def display_channels(self, channel_list=None):
        if channel_list == None:
            channel_list = self.channel_list
        print("CHANNELS:")
        for channel in channel_list:
            print("\t{}:".format(channel))
            for item in list(self.params['channels'][channel]):
                print("\t\t{}: {}".format(item, self.params['channels'][channel][item]))

    def info(self, channel_list=None):
        self.display_scan_strategy()
        if 'half_wave_plates' in self.params.keys():
            self.display_half_wave_plate()
            self.display_channels(channel_list)

    #*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
    # Path naming conventions
    #*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*

    def path_to_instrument_dir(self, instrument_name):
        return os.path.join(self.global_paths['instruments_dir'], instrument_name)

    def path_to_instrument_param_file(self, instrument_name, instrument_param_file_name="instrument_params.yaml"):
        """
        Default instrument param file name: instrument.yaml
        """
        return os.path.join(self.path_to_instrument_dir(instrument_name), instrument_param_file_name)

    #*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
    # Validating the parameters
    #*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
