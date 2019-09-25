from genesys.utilities import Generic_Class
import genesys.data_io.segment_distribution as segd
import numpy as np

config = Generic_Class()

#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
# Naming and general configurations for this simulation run
#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*

# This will be name of the parent directory under which all output data products will be written.
config.sim_tag = "sim_test"
# This will be the name of the directory under which the simulated timestream data will be written.
config.scan_tag = "scan_test"

# This will be prepended to each data segment written down. For example when resampling, the tag "resim_" can be prepended to signal.
config.special_tag = ""

# Options for input_pol_type = ["IQU", "QU", "I", "_QU", "noise_only"]
# Default -> "IQU"
# A _ means that the input map has a component which will not be read"
config.sim_pol_type = "IQU"

# coordinate_system can be "ecliptic" OR "galactic". This assumes that the input map is always in the galactic reference frame.
config.coordinate_system = "galactic"

# Focal plane configuration file. This is going to be a single file containing the specifications of each detector that is going to scan the sky. Among the specification, the noise and beam properties are given.
config.focal_plane_config = "litebird.baseline_op1.lb_focal_plane_baseline_op1"

config.notes = ""

#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
# Scan Strategy 
#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*

# Mention the module name in scan_strategy_module in a 'dot' separated fashion assuming genesys.scan_strategy to be the parent.
config.scan_strategy_config = "litebird.litebird_baseline"

#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
# Noise 
#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*

#Options for noise_type = ["none", "white", "1_over_f"]
config.noise_type = "white"

#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
# Data selection
#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*

# t_segment is the length of each data segment
band_detector_list = {'HFT_3' : ['0001a', '0001b', '0002a', '0002b']}
segment_length = 0.507292*24*60*60.0                      #seconds
segment_list = np.arange(4)

config.band_detector_segment_dict = segd.get_band_detector_time_dictionary_uniform(band_detector_list, segment_length, segment_list)

#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
# Read & Write
#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*

# Available options for timestream_data_products = ["signal", "pointing_vec", "pol_ang", "noise", "mask"]
config.ts_data_products = ["time", "signal", "pointing_vec", "pol_ang"]
