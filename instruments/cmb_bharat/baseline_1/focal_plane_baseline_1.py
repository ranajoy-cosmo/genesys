#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
# This module contains the instrumental specifications for a particular experiment.
#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*

from genesys.utilities import Generic_Class
import importlib

focal_plane = Generic_Class()

#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
# This section defines the general aspects of the instrument and the frequency bands
#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*

# Name for this instrument configuration
focal_plane.name = "CMB-Bharat_baseline_1"

focal_plane.band_list = ['LFT_1', 'LFT_2', 'LFT_3', 'LFT_4', 'LFT_5', 'LFT_6', 'LFT_7', 'LFT_8', 'LFT_9', 'HFT_1', 'HFT_2', '140GHz', 'HFT_4', 'HFT_5', 'HFT_6', 'HFT_7', 'HFT_8', 'HFT_9']

focal_plane.bands = {}

########## Band LFT_1
band_name = "LFT_1"
focal_plane.bands[band_name] = Generic_Class()
focal_plane.bands[band_name].band_frequency = 40.0
focal_plane.bands[band_name].band_width = 12.0
focal_plane.bands[band_name].beam_fwhm = 69.2
focal_plane.bands[band_name].num_det = 42.0
focal_plane.bands[band_name].det_NET = 87.0
focal_plane.bands[band_name].detector_yield = 0.8
focal_plane.bands[band_name].noise_margin = 1.15

########## Band LFT_2
band_name = "LFT_2"
focal_plane.bands[band_name] = Generic_Class()
focal_plane.bands[band_name].band_frequency = 50.0
focal_plane.bands[band_name].band_width = 15.0
focal_plane.bands[band_name].beam_fwhm = 56.9
focal_plane.bands[band_name].num_det = 56.0
focal_plane.bands[band_name].det_NET = 54.6
focal_plane.bands[band_name].detector_yield = 0.8
focal_plane.bands[band_name].noise_margin = 1.15

########## Band LFT_3
band_name = "LFT_3"
focal_plane.bands[band_name] = Generic_Class()
focal_plane.bands[band_name].band_frequency = 60.0
focal_plane.bands[band_name].band_width = 14.0
focal_plane.bands[band_name].beam_fwhm = 49.0
focal_plane.bands[band_name].num_det = 42.0
focal_plane.bands[band_name].det_NET = 48.7
focal_plane.bands[band_name].detector_yield = 0.8
focal_plane.bands[band_name].noise_margin = 1.15

########## Band LFT_4
band_name = "LFT_4"
focal_plane.bands[band_name] = Generic_Class()
focal_plane.bands[band_name].band_frequency = 68.0
focal_plane.bands[band_name].band_width = 16.0
focal_plane.bands[band_name].beam_fwhm = 40.8
focal_plane.bands[band_name].num_det = 56.0
focal_plane.bands[band_name].det_NET = 43.4
focal_plane.bands[band_name].detector_yield = 0.8
focal_plane.bands[band_name].noise_margin = 1.15

########## Band LFT_5
band_name = "LFT_5"
focal_plane.bands[band_name] = Generic_Class()
focal_plane.bands[band_name].band_frequency = 78.0
focal_plane.bands[band_name].band_width = 18.0
focal_plane.bands[band_name].beam_fwhm = 36.1
focal_plane.bands[band_name].num_det = 42.0
focal_plane.bands[band_name].det_NET = 40.0
focal_plane.bands[band_name].detector_yield = 0.8
focal_plane.bands[band_name].noise_margin = 1.15

########## Band LFT_6
band_name = "LFT_6"
focal_plane.bands[band_name] = Generic_Class()
focal_plane.bands[band_name].band_frequency = 89.0
focal_plane.bands[band_name].band_width = 20.0
focal_plane.bands[band_name].beam_fwhm = 32.3
focal_plane.bands[band_name].num_det = 56.0
focal_plane.bands[band_name].det_NET = 36.8
focal_plane.bands[band_name].detector_yield = 0.8
focal_plane.bands[band_name].noise_margin = 1.15

########## Band LFT_7
band_name = "LFT_7"
focal_plane.bands[band_name] = Generic_Class()
focal_plane.bands[band_name].band_frequency = 100.0
focal_plane.bands[band_name].band_width = 12.0
focal_plane.bands[band_name].beam_fwhm = 27.7
focal_plane.bands[band_name].num_det = 114.0
focal_plane.bands[band_name].det_NET = 40.1
focal_plane.bands[band_name].detector_yield = 0.8
focal_plane.bands[band_name].noise_margin = 1.15

########## Band LFT_8
band_name = "LFT_8"
focal_plane.bands[band_name] = Generic_Class()
focal_plane.bands[band_name].band_frequency = 119.0
focal_plane.bands[band_name].band_width = 36.0
focal_plane.bands[band_name].beam_fwhm = 23.7
focal_plane.bands[band_name].num_det = 114.0
focal_plane.bands[band_name].det_NET = 31.3
focal_plane.bands[band_name].detector_yield = 0.8
focal_plane.bands[band_name].noise_margin = 1.15

########## Band LFT_9
band_name = "LFT_9"
focal_plane.bands[band_name] = Generic_Class()
focal_plane.bands[band_name].band_frequency = 140.0
focal_plane.bands[band_name].band_width = 42.0
focal_plane.bands[band_name].beam_fwhm = 20.7
focal_plane.bands[band_name].num_det = 114.0
focal_plane.bands[band_name].det_NET = 29.6
focal_plane.bands[band_name].detector_yield = 0.8
focal_plane.bands[band_name].noise_margin = 1.15

########## Band HFT_1
band_name = "HFT_1"
focal_plane.bands[band_name] = Generic_Class()
focal_plane.bands[band_name].band_frequency = 100.0
focal_plane.bands[band_name].band_width = 12.0
focal_plane.bands[band_name].beam_fwhm = 37.0
focal_plane.bands[band_name].num_det = 222.0
focal_plane.bands[band_name].det_NET = 53.7
focal_plane.bands[band_name].detector_yield = 0.8
focal_plane.bands[band_name].noise_margin = 1.15

########## Band HFT_2
band_name = "HFT_2"
focal_plane.bands[band_name] = Generic_Class()
focal_plane.bands[band_name].band_frequency = 119.0
focal_plane.bands[band_name].band_width = 31.6
focal_plane.bands[band_name].beam_fwhm = 31.6
focal_plane.bands[band_name].num_det = 148.0
focal_plane.bands[band_name].det_NET = 38.4
focal_plane.bands[band_name].detector_yield = 0.8
focal_plane.bands[band_name].noise_margin = 1.15

########## Band 140GHz
band_name = "140GHz"
focal_plane.bands[band_name] = Generic_Class()
focal_plane.bands[band_name].band_id = band_name
focal_plane.bands[band_name].band_frequency = 140.0
focal_plane.bands[band_name].band_width = 42.0
focal_plane.bands[band_name].beam_fwhm = 27.6
focal_plane.bands[band_name].num_det = 222.0
focal_plane.bands[band_name].det_NET = 33.9
focal_plane.bands[band_name].detector_yield = 0.8
focal_plane.bands[band_name].noise_margin = 1.15
focal_plane.bands[band_name].detectors = importlib.import_module("genesys.instrument.cmb_bharat.baseline_1.140GHz").detectors

########## Band HFT_4
band_name = "HFT_4"
focal_plane.bands[band_name] = Generic_Class()
focal_plane.bands[band_name].band_frequency = 166.0
focal_plane.bands[band_name].band_width = 50.0
focal_plane.bands[band_name].beam_fwhm = 24.2
focal_plane.bands[band_name].num_det = 148.0
focal_plane.bands[band_name].det_NET = 31.8
focal_plane.bands[band_name].detector_yield = 0.8
focal_plane.bands[band_name].noise_margin = 1.15

########## Band HFT_5
band_name = "HFT_5"
focal_plane.bands[band_name] = Generic_Class()
focal_plane.bands[band_name].band_frequency = 195.0
focal_plane.bands[band_name].band_width = 58.0
focal_plane.bands[band_name].beam_fwhm = 21.7
focal_plane.bands[band_name].num_det = 222.0
focal_plane.bands[band_name].det_NET = 32.4
focal_plane.bands[band_name].detector_yield = 0.8
focal_plane.bands[band_name].noise_margin = 1.15

########## Band HFT_6
band_name = "HFT_6"
focal_plane.bands[band_name] = Generic_Class()
focal_plane.bands[band_name].band_frequency = 235.0
focal_plane.bands[band_name].band_width = 70.0
focal_plane.bands[band_name].beam_fwhm = 19.6
focal_plane.bands[band_name].num_det = 148.0
focal_plane.bands[band_name].det_NET = 36.4
focal_plane.bands[band_name].detector_yield = 0.8
focal_plane.bands[band_name].noise_margin = 1.15

########## Band HFT_7
band_name = "HFT_7"
focal_plane.bands[band_name] = Generic_Class()
focal_plane.bands[band_name].band_frequency = 280.0
focal_plane.bands[band_name].band_width = 84.0
focal_plane.bands[band_name].beam_fwhm = 13.2
focal_plane.bands[band_name].num_det = 338.0
focal_plane.bands[band_name].det_NET = 62.2
focal_plane.bands[band_name].detector_yield = 0.8
focal_plane.bands[band_name].noise_margin = 1.15

########## Band HFT_8
band_name = "HFT_8"
focal_plane.bands[band_name] = Generic_Class()
focal_plane.bands[band_name].band_frequency = 337.0
focal_plane.bands[band_name].band_width = 102.0
focal_plane.bands[band_name].beam_fwhm = 11.2
focal_plane.bands[band_name].num_det = 338.0
focal_plane.bands[band_name].det_NET = 78.3
focal_plane.bands[band_name].detector_yield = 0.8
focal_plane.bands[band_name].noise_margin = 1.15

########## Band HFT_9
band_name = "HFT_9"
focal_plane.bands[band_name] = Generic_Class()
focal_plane.bands[band_name].band_frequency = 402.0
focal_plane.bands[band_name].band_width = 92.0
focal_plane.bands[band_name].beam_fwhm = 9.7
focal_plane.bands[band_name].num_det = 338.0
focal_plane.bands[band_name].det_NET = 134.5
focal_plane.bands[band_name].detector_yield = 0.8
focal_plane.bands[band_name].noise_margin = 1.15

# Unit: GHz
focal_plane.band_frequency = {'LFT_1' : 40.0, 'LFT_2' : 50.0, 'LFT_3' : 60.0, 'LFT_4' : 68.0, 'LFT_5' : 78.0, 'LFT_6' : 89.0, 'LFT_7' : 100.0, 'LFT_8' : 119.0, 'LFT_9' : 140.0, 'HFT_1' : 100.0, 'HFT_2' : 119.0, 'HFT_3' : 140.0, 'HFT_4' : 166.0, 'HFT_5' : 195.0, 'HFT_6' : 235.0, 'HFT_7' : 280.0, 'HFT_8' : 337.0, 'HFT_9' : 402.0}

# Unit: GHz
focal_plane.band_width = {'LFT_1' : 12.0, 'LFT_2' : 15.0, 'LFT_3' : 14.0, 'LFT_4' : 16.0, 'LFT_5' : 18.0, 'LFT_6' : 20.0, 'LFT_7' : 12.0, 'LFT_8' : 36.0, 'LFT_9' : 42.0, 'HFT_1' : 12.0, 'HFT_2' : 36.0, 'HFT_3' : 42.0, 'HFT_4' : 50.0, 'HFT_5' : 58.0, 'HFT_6' : 70.0, 'HFT_7' : 84.0, 'HFT_8' : 102.0, 'HFT_9' : 92.0}

# Unit: arcmin
focal_plane.beam_fwhm = {'LFT_1' : 69.2, 'LFT_2' : 56.9, 'LFT_3' : 49.0, 'LFT_4' : 40.8, 'LFT_5' : 36.1, 'LFT_6' : 32.3, 'LFT_7' : 27.7, 'LFT_8' : 23.7, 'LFT_9' : 20.7, 'HFT_1' : 37.0, 'HFT_2' : 31.6, 'HFT_3' : 27.6, 'HFT_4' : 24.2, 'HFT_5' : 21.7, 'HFT_6' : 19.6, 'HFT_7' : 13.2, 'HFT_8' : 11.2, 'HFT_9' : 9.7}

focal_plane.num_det = {'LFT_1' : 42.0, 'LFT_2' : 56.0, 'LFT_3' : 42.0, 'LFT_4' : 56.0, 'LFT_5' : 42.0, 'LFT_6' : 56.0, 'LFT_7' : 114.0, 'LFT_8' : 114.0, 'LFT_9' : 114.0, 'HFT_1' : 222.0, 'HFT_2' : 148.0, 'HFT_3' : 222.0, 'HFT_4' : 148.0, 'HFT_5' : 222.0, 'HFT_6' : 148.0, 'HFT_7' : 338.0, 'HFT_8' : 338.0, 'HFT_9' : 338.0}

# Unit: uK.sqrt(s)
focal_plane.det_NET = {'LFT_1' : 87.0, 'LFT_2' : 54.6, 'LFT_3' : 48.7, 'LFT_4' : 43.4, 'LFT_5' : 40.0, 'LFT_6' : 36.8, 'LFT_7' : 40.1, 'LFT_8' : 31.3, 'LFT_9' : 29.6, 'HFT_1' : 53.7, 'HFT_2' : 38.4, 'HFT_3' : 33.9, 'HFT_4' : 31.8, 'HFT_5' : 32.4, 'HFT_6' : 36.4, 'HFT_7' : 62.2, 'HFT_8' : 78.3, 'HFT_9' : 134.5}

focal_plane.detector_yield = 0.8
focal_plane.noise_margin = 1.15

#{'LFT_1' : , 'LFT_2' : , 'LFT_3' : , 'LFT_4' : , 'LFT_5' : , 'LFT_6' : , 'LFT_7' : , 'LFT_8' : , 'LFT_9' : , 'HFT_1' : , 'HFT_2' : , 'HFT_3' : , 'HFT_4' : , 'HFT_5' : , 'HFT_6' : , 'HFT_7' : , 'HFT_8' : , 'HFT_9' : }
