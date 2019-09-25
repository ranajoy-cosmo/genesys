from genesys.utilities import Generic_Class

scan_strategy = Generic_Class()

scan_strategy.scan_strategy_name = "CMB-Bharat baseline"

scan_strategy.t_year = 365.25*24*60*60          # seconds

scan_strategy.alpha = 50.0                      # degrees
scan_strategy.beta = 45.0                       # degrees
scan_strategy.t_prec = 4*24*60*60.0                  # seconds
scan_strategy.t_spin = 60.0                    # seconds

scan_strategy.polarisation_modulation = "instrumental_scanning"

#  scan_strategy.polarisation_modulation = "continuous_HWP"
#  scan_strategy.HWP_phase_ini = None              # degrees
#  scan_strategy.HWP_angular_speed = None          # RPM
#
#  scan_strategy.polarisation_modulation = "stepped_HWP"
#  scan_strategy.HWP_phase_ini = None              # degrees
#  scan_strategy.HWP_step = None                   # degrees
#  scan_strategy.HWP_step_duration = None          # seconds
