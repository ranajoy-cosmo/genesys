import genesys.data_io.segment_distribution as segd
import numpy as np
from genesys.map_maker.config_files.default_config import config

#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
# Data selection
#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*

# t_segment is the length of each data segment
band_detector_list = {'140GHz' : ['0001a', '0001b', '0002a', '0002b']}
segment_length = 0.71337890625*24*60*60.0 / 2.0                     #seconds
segment_list = np.arange(64*8*2)

config.band_detector_segment_dict = segd.get_band_detector_time_dictionary_uniform(band_detector_list, segment_length, segment_list)

config.focal_plane_config = "cmb_bharat.baseline_1.focal_plane_baseline_1"
config.scan_strategy_config = "cmb_bharat.cmb_bharat_baseline"
