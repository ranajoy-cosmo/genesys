import numpy as np
import copy
from genesys import Genesys_Class

"""
Class for handling timestream objects
"""

class Timestream(Genesys_Class):
    def __init__(self, other=None, sim_params=None, detector_params=None, segment=None):
        if other is not None:
            self.copy_attributes(other=other)
        elif detector_params is not None:
            self.configure_timestream(sim_params, detector_params, segment)
        else:
            self.ts = {}

    def configure_timestream(self, sim_params, detector_params, segment):
        pass
    #*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*

    def data_items(self):
        return list(self.data.keys())

    def get_note(self):
        print self.note

    def set_note(self, note):
        self.note = note

    #*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*

    def add_empty_data_items(self, data_items):
        current_data_items = self.data_items()
        for data_item in data_items:
            if not data_item in current_data_items: 
                self.data[data_item] = None
            else
                pass

    def overwrite_data_items_with_empty(self, data_items):
        for data_item in data_items:
            self.data[data_item] = None

    def remove_data_items(self, data_items):
        current_data_items = self.data_items()
        for data_item in data_items:
            if data_item in current_data_items: 
                del self.data[data_item]
            else:
                pass

    def add_data_item(self, data_item, data_values=None):
        self.data[data_item] = data_values


