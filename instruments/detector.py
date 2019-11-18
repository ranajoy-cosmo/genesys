import os
from genesys import Genesys_Class

class Detector(Genesys_Class):
    def __init__(self, other=None, instrument_name=None, channel_name=None, detector_name=None):
        if other is not None:
            self.copy_other(other)
        elif detector_name is not None:
            self.configure_detector(instrument_name, channel_name, detector_name)
        else:
            pass

    def configure_detector(instrument_name, channel_name, detector_name):
        instrument_param_file_path = self.path_to_instrument_param_file(instrument_name)
        inst_params = self.load_param_file(file_path=instrument_param_file_path)
        channel_file_name = inst_params['channels'][channel_name]['channel_param_file']
        channel_param_file_path = self.path_to_channel_param_file(instrument_name, channel_file_name)
        channel_params = self.load_param_file(file_path=channel_param_file_path)
        self.params = channel_params[detector_name]

    #*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
    # Path naming conventions
    #*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*

    def path_to_instrument_dir(self, instrument_name):
        return os.path.join(self.global_paths['instruments_dir'], instrument_name)

    def path_to_instrument_param_file(self, instrument_name):
        return os.path.join(self.path_to_instrument_dir(instrument_name), "instrument_params.yaml")

    def path_to_channel_param_file(self, instrument_name, channel_file_name):
        return os.path.join(self.path_to_instrument_dir(instrument_name), "instrument_params.yaml", channel_file_name)
