import numpy as np
import copy
from genesys import Genesys_Class

"""
Class for handling timestream objects
"""

class TStream(Genesys_Class):
    def __init__(self, other=None, sim_config=None, detector_params=None, segment=None):
        if other is not None:
            self.copy_attributes(other=other)
        elif detector_params is not None:
            self.initialise_timestream(sim_config, detector_params, segment)
        else:
            self.ts = {}

    def initialise_timestream(self, sim_config, detector_params, segment):
        self.ts = {}.fromkeys(sim_config['ts_data_products'])
