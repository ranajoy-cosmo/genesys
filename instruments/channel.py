import importlib
import os
import copy

class Channel():
    """
    Class containing the specifications of a particular frequency channel
    Data members:
        name: string name
        frequency: central frequency in GHz
        bandpass:
        num_det: # of detectors
        beam_fwhm:
        detector_NET:
        det_yield:
        noise_margin:
        detector_dict: Dictionary of all the individual detectors
    """
    def __init__(self, other=None, channel_file=None, channel_name=None):
        if other != None:
            self.copy_channel(other)
        elif channel_file != None:
            self.read_channel_from_file(channel_file, channel_name)
        else:
            self.fill_empty_channel(self)

        def copy_channel(self, other):

