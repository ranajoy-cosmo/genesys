from genesys.utilities import Generic_Class
import genesys.data_io.segment_distribution as segd
import numpy as np

config = Generic_Class()

#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
# Sim Action
#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
#config.simulate_ts = True/False. True means it will simulate the entire timestream signal and pointing from scratch
#if simulate_ts is True, this is where the timestream data will be written if it is asked for. If False, this is where the timestream data will be read from.
config.simulate_ts = True

# This will be name of the parent directory under which all output data products will be written.
config.sim_tag = "sim_test"
# This will be the name of the directory under which the simulated timestream data will be written.
config.scan_tag = "scan_test"
# This will be the name of the directory under which the reconstructed maps will be written.
config.recon_tag = "recon_test"

# This will be prepended to each data segment written down. For example when resampling, the tag "resim_" can be prepended to signal.
config.special_tag = ""

# Options for recon_pol_type = ["IQU", "QU", "I", "_QU", "noise_only"]
# Default -> "IQU"
# A _ means that the input map has a component which will not be read"
config.recon_pol_type = "IQU"

# noise_only is set to True when the map of the noise separately is required
config.noise_only = False

config.pair_difference = False
config.subtract_template = False

config.notes = "Testing"
#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
# Output Resolution
#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
#config.nside_out : nside of the output map
config.nside_out = 1024

#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
# Data selection
#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*

# t_segment is the length of each data segment
band_detector_list = {'HFT_3' : ['0001a', '0001b', '0002a', '0002b']}
segment_length = 0.71337890625*24*60*60.0                      #seconds
segment_list = np.arange(64*8)

config.band_detector_segment_dict = segd.get_band_detector_time_dictionary_uniform(band_detector_list, segment_length, segment_list)

#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
# Sim Action (Only if simulate_ts is True)
#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*

# Options for input_pol_type = ["IQU", "QU", "I", "_QU", "noise_only"]
# Default -> "IQU"
# A _ means that the input map has a component which will not be read"
config.sim_pol_type = "IQU"

# coordinate_system can be "ecliptic" OR "galactic". This assumes that the input map is always in the galactic reference frame.
config.coordinate_system = "galactic"

# Focal plane configuration file. This is going to be a single file containing the specifications of each detector that is going to scan the sky. Among the specification, the noise and beam properties are given.
config.focal_plane_config = "litebird.baseline_op1.lb_focal_plane_baseline_op1"

#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
# Scan Strategy (Only if simulate_ts is True but currently needs to be loaded always for the detector object to initialise correctly)
#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*

# Mention the module name in scan_strategy_module in a 'dot' separated fashion assuming genesys.scan_strategy to be the parent.
config.scan_strategy_config = "litebird.litebird_baseline"

#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
# Noise (Only if simulate_ts is True)
#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*

#Options for noise_type = ["none", "white", "1_over_f"]
config.noise_type = "white"

#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
# Read & Write (Only if simulate_ts is True)
#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*

# Available options for timestream_data_products = ["signal", "pointing_vec", "pol_ang", "noise", "mask"]
#  config.ts_data_products = ["time", "signal", "pointing_vec", "pol_ang"]
config.ts_data_products = []
