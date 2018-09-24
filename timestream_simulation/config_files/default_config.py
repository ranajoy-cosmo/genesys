from genesys.utilities import Generic_Class
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

config.simulate_data = True

# Options for input_pol_type = ["IQU", "QU", "I", "_QU", "noise_only"]
# Default -> "IQU"
# A _ means that the input map has a component which will not be read"
config.sim_pol_type = "IQU"

# coordinate_system can be "ecliptic" OR "galactic". This assumes that the input map is always in the galactic reference frame.
config.coordinate_system = "galactic"

# Focal plane configuration file. This is going to be a single file containing the specifications of each detector that is going to scan the sky. Among the specification, the noise and beam properties are given.
config.focal_plane_config = "litebird.baseline_op1.lb_focal_plane_baseline_op1_detailed"

# tod_type can be "signal" OR "gradient"
config.tod_type = "signal"
# if tod_type is chosen as "gradient", gradient_type can be the following combinations
# Zero order gradient. Useful for BPMM templates
# "zero_order"
# First order derivatives
# "gradient_co" OR "gradient_cross" OR ["gradient_co", "gradient_cross"]
# Second order derivatives
# "gradient_coXgradient_co" OR "gradient_crossXgradient_cross" OR ["gradient_coXgradient_co", "gradient_crossXgradient_cross"]
# "gradient_coXgradient_cross"
# Combinations of the gradients
config.gradient_type = []

config.notes = ""

#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
# Scan Strategy 
#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*

# If the scan strategy is defined in a file, mention the module name in scan_strategy_module in a 'dot' separated fashion assuming genesys.scan_strategy to be the parent.
# If you want to custom set the scan strategy here, set scan_strategy_module as None.
# Important: If scan_strategy_module is set to None, the custom strategy will be accepted.
config.scan_strategy_module = "litebird.litebird_baseline"

#* Example #*#*#*#*#*
# config.t_year = 365.25*24*60*60                   # seconds
# config.scan_strategy_name = "Example"
# config.t_prec = 4*24*60*60.0                      # seconds
# config.t_spin = 120.0                             # seconds
# config.sampling_rate = 170                        # Hz
#
# config.alpha = 30.0                               # degrees
# config.beta = 61.0                                # degrees
# config.scan_straegy_note = ""
#*#*#*#*#*#*#*#*#*#*

config.oversampling_rate = 1

config.nside_in = 1024

#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
# Noise 
#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*

#Options for noise_type = ["none", "white", "1_over_f"]
config.noise_type = "white"

#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
# Data selection
#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*

# t_segment is the length of each data segment
config.t_segment = 12*60*60.0                      #seconds

segment_list = np.arange(8)

config.detector_segment_dict = {'0001a': segment_list, '0001b': segment_list, '0002a': segment_list, '0002b': segment_list}

#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
# Read & Write
#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*

# Available options for timestream_data_products = ["signal", "pointing_vec", "pol_ang", "noise", "mask"]
config.ts_data_products = ["signal", "pointing_vec", "pol_ang"]

# Set output_dir to None if you want the global_path defined data_dir as your output_dir
config.output_dir = None

config.write_beam = False

#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
# Beam Parameters
#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*

#beam_type can be "pencil" OR "full_simulated" OR "from_file" 
config.beam_type = "pencil"
config.beam_cutoff = 4                                                             #fwhm
config.check_normalisation = True
