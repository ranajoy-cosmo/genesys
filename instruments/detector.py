import os
import importlib

class Detector():
    def __init__(self, other=None, instrument_name=None, channel_name=None, detector_name=None):
        if other != None:
            self.copy_detector(other)
        elif file_name != None:
            self.load_detector_from_file(instrument_name, channel_name, detector_name)
        else:
            pass

    def copy_detector(self, other):
        pass

    def load_detector_from_file(instrument_name, channel_name, detector_name):

