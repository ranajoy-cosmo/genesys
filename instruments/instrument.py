import os
from genesys import Genesys_Class
from genesys.instruments.channel import Channel

class Instrument(Genesys_Class):
    """
    Class for handling instrument parameter loading and handling
    The instrument params are written in  <instrument_dir> directory, and contains the following files
        instrument.yaml -> Top level instrument specifications, scan strategy, polarisation modulator
        <channel_name>.yaml -> Individual detector parameters
    """
    def __init__(self, instrument_dir=None, other=None):
        """
        Order of preference: instrument_dir > other
        """
        if instrument_dir != None:
            self.instrument_dir = os.path.join(self.global_paths['base_dir'], 'instruments', instrument_dir)
            self.configure_instrument()
        elif other != None:
            self.copy_attributes(other)
        else:
            self.params = {}

    def configure_instrument(self):
        """
        READS IN THE INSTRUMENT PARAMETERS FROM THE FILES IN THE self.intstrument_dir
        """
        instrument_param_file_path = os.path.join(self.instrument_dir, 'instrument.yaml')
        self.params = self.load_param_file(file_path=instrument_param_file_path)
        # TO DO: WRITE A MODULE TO VALIDATE THE LOADED PARAMS

    #*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
    # Paramter display methods
    #*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*

    def display_scan_strategy(self):
        unit_dict = {'alpha': 'degrees', 'beta': 'degrees', 't_precession': 'seconds', 't_spin': 'seconds', 'duration': 'years'}
        scan_strategy = self.params['scan_strategy']
        self.prompt("\nSCAN STRATEGY:")
        for item in scan_strategy:
            self.prompt(f"\t{item}: {scan_strategy[item]} {unit_dict[item]}")

    def display_half_wave_plates(self):
        unit_dict = {'spin_rate': 'rpm', 'rotation': ''}
        hwps = self.params['HWP']
        self.prompt("\nHALF WAVE PLATES:")
        for hwp in hwps:
            self.prompt(f"\t{hwp}:")
            for item in list(hwps[hwp]):
                self.prompt(f"\t\t{item}: {hwps[hwp][item]} {unit_dict[item]}")

    def display_channels(self, channel_list=None):
        unit_dict = {'central_frequency': 'GHz', 'bandwidth': 'GHz', 'beam_fwhm': 'arcmins', 'sampling_rate': 'Hz', 'pixel_size': 'mm', 'num_detectors': '', 'noise_type': '', 'white_noise_rms': 'uk.sqrt(s)', 'detector_yield': '', 'noise_margin': '', 'pol_modulation': '', 'channel_param_file': '', 'input_map_file': ''}
        channels = self.params['channels']
        if channel_list == None:
            channel_list = channels.keys()
        self.prompt("\nCHANNELS:")
        for channel in channel_list:
            self.prompt(f"\t{channel}:")
            for item in list(channels[channel]):
                if item != 'noise':
                    self.prompt(f"\t\t{item}: {channels[channel][item]} {unit_dict[item]}")
                else:
                    self.prompt("\t\tnoise/detector")
                    for noise_item in list(channels[channel]['noise']):
                        self.prompt(f"\t\t\t{noise_item}: {channels[channel]['noise'][noise_item]} {unit_dict[noise_item]}")

    def info(self, channel_list=None):
        self.prompt(f"\nInstrument name: {self.params['instrument_name']}")
        self.prompt(f"Version: {self.params['version']}")
        self.display_scan_strategy()
        if 'HWP' in self.params:
            self.display_half_wave_plates()
        self.display_channels(channel_list)

    #*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
    # Loading channel
    #*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*

    def get_channel(self, channel_name):
        channel_obj = Channel(self, channel_name)
        return channel_obj

