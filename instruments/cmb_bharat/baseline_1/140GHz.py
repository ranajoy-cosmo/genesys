from genesys.utilities import Generic_Class
from genesys.data_io.data_io import get_detector_id

#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
# This section defines the individual detectors of the HFT_3 band
#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*

band_name = "140GHz"
detectors = {}

##### Detector 0001a
detector_name = '0001a'
detectors[detector_name] = Generic_Class()
detectors[detector_name].detector_id = get_detector_id(band_name, detector_name)
detectors[detector_name].beam_fwhm_major = 27.6         # arcmins
detectors[detector_name].beam_fwhm_minor = 27.6         # arcmins
detectors[detector_name].beam_angle = 0.0               # degrees
detectors[detector_name].pol_phase_ini = 0.0            # degrees
detectors[detector_name].focal_plane_pos = [0.0,0.0]    # arcmins
detectors[detector_name].offset = [0.0,0.0]             # arcsecs
detectors[detector_name].sampling_rate = 200             # Hz
detectors[detector_name].input_sky_map = "/mn/stornext/u3/ranajoyb/genesys/maps/map_files/test_map.fits"
detectors[detector_name].white_noise_sigma = 0.0       # uK / sqrt(sec)
detectors[detector_name].f_knee = None                  # mHz
detectors[detector_name].one_over_f_alpha = None
detectors[detector_name].white_noise_seed = None
detectors[detector_name].low_f_noise_seed = None
detectors[detector_name].input_beam_file = None

##### Detector 0001b
detector_name = '0001b'
detectors[detector_name] = Generic_Class()
detectors[detector_name].detector_id = get_detector_id(band_name, detector_name)
detectors[detector_name].beam_fwhm_major = 27.6
detectors[detector_name].beam_fwhm_minor = 27.6
detectors[detector_name].beam_angle = 0.0
detectors[detector_name].pol_phase_ini = 90.0
detectors[detector_name].focal_plane_pos = [0.0,0.0]
detectors[detector_name].offset = [0.0,0.0]
detectors[detector_name].sampling_rate = 200             # Hz
detectors[detector_name].input_sky_map = "/mn/stornext/u3/ranajoyb/genesys/maps/map_files/test_map.fits"
detectors[detector_name].white_noise_sigma = 0.0
detectors[detector_name].f_knee = None
detectors[detector_name].one_over_f_alpha = None
detectors[detector_name].white_noise_seed = None
detectors[detector_name].low_f_noise_seed = None
detectors[detector_name].input_beam_file = None

##### Detector 0001c
detector_name = '0002a'
detectors[detector_name] = Generic_Class()
detectors[detector_name].detector_id = get_detector_id(band_name, detector_name)
detectors[detector_name].beam_fwhm_major = 27.6
detectors[detector_name].beam_fwhm_minor = 27.6
detectors[detector_name].beam_angle = 0.0
detectors[detector_name].pol_phase_ini = 45.0
detectors[detector_name].focal_plane_pos = [0.0,0.0]
detectors[detector_name].offset = [0.0,0.0]
detectors[detector_name].sampling_rate = 200             # Hz
detectors[detector_name].input_sky_map = "/mn/stornext/u3/ranajoyb/genesys/maps/map_files/test_map.fits"
detectors[detector_name].white_noise_sigma = 0.0
detectors[detector_name].f_knee = None
detectors[detector_name].one_over_f_alpha = None
detectors[detector_name].white_noise_seed = None
detectors[detector_name].low_f_noise_seed = None
detectors[detector_name].input_beam_file = None

##### Detector 0001d
detector_name = '0002b'
detectors[detector_name] = Generic_Class()
detectors[detector_name].detector_id = get_detector_id(band_name, detector_name)
detectors[detector_name].beam_fwhm_major = 27.6
detectors[detector_name].beam_fwhm_minor = 27.6
detectors[detector_name].beam_angle = 0.0
detectors[detector_name].pol_phase_ini = 135.0
detectors[detector_name].focal_plane_pos = [0.0,0.0]
detectors[detector_name].offset = [0.0,0.0]
detectors[detector_name].sampling_rate = 200             # Hz
detectors[detector_name].input_sky_map = "/mn/stornext/u3/ranajoyb/genesys/maps/map_files/test_map.fits"
detectors[detector_name].white_noise_sigma = 0.0
detectors[detector_name].f_knee = None
detectors[detector_name].one_over_f_alpha = None
detectors[detector_name].white_noise_seed = None
detectors[detector_name].low_f_noise_seed = None
detectors[detector_name].input_beam_file = None
