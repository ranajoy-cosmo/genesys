import os
import h5py
from genesys import Genesys_Class
from genesys.instruments.channel import Channel
from .default_units import unit_dict

class Instrument(Genesys_Class):
    """
    Class for handling instrument parameter loading and handling
    The instrument params are written in  <instrument_dir> directory, and contains the following files
        instrument.yaml -> Top level instrument specifications, scan strategy, polarisation modulator
        <channel_name>.yaml -> Individual detector parameters
    """
    def __init__(self, instrument_dir=None):
        """
        instrument_dir relative to genesys/instruments
        """
        if instrument_dir != None:
            self.instrument_dir = os.path.join(self.global_paths['base_dir'], 'instruments', instrument_dir)
            instrument_param_file_path = os.path.join(self.instrument_dir, 'instrument.yaml')
            self.params = self.load_param_file(file_path=instrument_param_file_path)
        else:
            self.params = {}

    def get_satellite_velocity(self):
        vel_file_path = os.path.join(self.global_paths['data_dir'], self.params['velocity']['velocity_file'])
        with h5py.File(vel_file_path, 'r') as f_vel:
            time = f_vel['time'][:]
            xvel = f_vel['xvel'][:]
            yvel = f_vel['yvel'][:]
            zvel = f_vel['zvel'][:]
        return time, xvel, yvel, zvel

    #*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
    # Paramter display methods
    #*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*

    def display_scan_strategy(self):
        scan_strategy = self.params['scan_strategy']
        self.prompt("\nScan Strategy:")
        for item in scan_strategy:
            self.prompt(f"\t{item:<35}{scan_strategy[item]} {unit_dict[item]}")

    def display_half_wave_plates(self):
        if self.params['HWP'] == None:
            self.prompt("HWP absent")
        else:
            hwps = self.params['HWP']
            self.prompt("\nHALF WAVE PLATES:")
            for hwp in hwps:
                self.prompt(f"\t{hwp}:")
                for item in hwps[hwp]:
                    self.prompt(f"\t\t{item:<27}{hwps[hwp][item]} {unit_dict[item]}")

    def display_channels(self):
        channels = self.params['channels']
        self.prompt("\nCHANNELS:")
        for channel in channels:
            self.prompt(f"\t{channel}:")
            for item in channels[channel]:
                if item == 'noise':
                    self.prompt("\t\tnoise/detector:")
                    for noise_item in list(channels[channel]['noise']):
                        self.prompt(f"\t\t\t{noise_item:<19}{channels[channel]['noise'][noise_item]} {unit_dict[noise_item]}")
                else:
                    self.prompt(f"\t\t{item:<27}{channels[channel][item]} {unit_dict[item]}")

    def info(self, channel_list=None):
        self.prompt(f"\n{'Instrument name':<35}{self.params['instrument_name']}")
        self.prompt(f"{'Version':<35}{self.params['version']}")
        self.display_scan_strategy()
        self.display_half_wave_plates()
        self.display_channels()

    #*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
    # Loading channel
    #*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*

    def get_channel(self, channel_name):
        channel_obj = Channel(self, channel_name)
        return channel_obj
