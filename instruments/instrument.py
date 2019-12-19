import os
from genesys import Genesys_Class
from genesys.instruments.channel import Channel

class Instrument(Genesys_Class):
    """
    class for handling instrument parameter loading and handling
    The instrument is described in a instrument_name folder.
    The folder contains the following files
        instrument.yaml -> Top level instrument specifications, scan strategies, polarisation modulator
        <Channel_name>.yaml -> Individual detector parameters.
    """
    def __init__(self, instrument_dir=None, other=None):
        """
        Constructor for Instrument class
        If other is given, the instance variables of other are copied into self
        If instrument_name is given, the intrument parameters are loaded from the corresponding directory
        """
        if instrument_dir != None:
            self.instrument_dir = instrument_dir
            self.configure_instrument()
        elif other != None:
            self.copy_attributes(other)
        else:
            pass

    def configure_instrument(self):
        """
        Reads in the instrument parameters from the files in the intstrument_folder
        """
        instrument_param_file_path = self.path_to_instrument_param_file()
        self.params = self.load_param_file(file_path=instrument_param_file_path)
        if 'HWP' in self.params.keys():
            self.hwp_list = list(self.params['HWP'].keys())
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
        hwps = self.params['HWP']
        print("HALF WAVE PLATES:")
        for hwp in list(hwps.keys()):
            print("\t{}:".format(hwp))
            for item in list(hwps[hwp].keys()):
                print("\t\t{}: {}".format(item, hwps[hwp][item]))

    def display_channels(self, channel_list=None):
        unit_dict = {'central_frequency': 'GHz', 'bandwidth': 'GHz', 'beam_fwhm': 'arcmins', 'sampling_rate': 'Hz', 'pixel_size': 'mm', 'num_detectors': '', 'noise_type': '', 'white_noise_rms': 'uk.sqrt(s)', 'detector_yield': '', 'noise_margin': '', 'pol_modulation': '', 'channel_param_file': ''}
        if channel_list == None:
            channel_list = self.channel_list
        print("CHANNELS:")
        for channel in channel_list:
            print("\t{}:".format(channel))
            for item in list(self.params['channels'][channel]):
                if item != 'noise':
                    print("\t\t{}: {} {}".format(item, self.params['channels'][channel][item], unit_dict[item]))
                else:
                    print("\t\tnoise/detector")
                    for noise_item in list(self.params['channels'][channel]['noise']):
                        print("\t\t\t{}: {} {}".format(noise_item, self.params['channels'][channel]['noise'][noise_item], unit_dict[noise_item]))

    def info(self, channel_list=None):
        self.display_scan_strategy()
        if 'HWP' in self.params.keys():
            self.display_half_wave_plate()
        self.display_channels(channel_list)

    #*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
    # Path naming conventions
    #*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*

    def path_to_instrument_dir(self):
        return os.path.join(self.global_paths['instruments_dir'], self.instrument_dir)

    def path_to_instrument_param_file(self):
        """
        Default instrument param file name: instrument.yaml
        """
        return os.path.join(self.path_to_instrument_dir(), "instrument_params.yaml")

    #*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
    # Loading channel
    #*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*

    def get_channel_object(self, channel_name):
        channel_obj = Channel(self, channel_name)
        return channel_obj

