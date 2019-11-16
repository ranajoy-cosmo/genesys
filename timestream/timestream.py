import numpy as np
import copy

"""
Class for handling timestream objects
"""

class Timestream():
    def __init__(self, other=None, filename=None, data_items=[], note=None):
        if other != None:
            self.copy_timestream(other)
        elif filename != None:
            pass
        else:
            self.data = dict.fromkeys(data_items)
            self.note = note

    def copy_timestream(self, other):
        self.data = copy.deepcopy(other.data)
        self.note = copy.deepcopy(other.note)

    def read_ts_from_file(self, file_name, data_items=None):
        pass

    def return_copy(self):
        return Timestream(self)

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


