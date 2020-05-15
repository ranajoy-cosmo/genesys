import numpy as np
import numpy.fft as fft
import math
from genesys import Genesys_Class
from genesys.numerical import Unit_Converter

uc = Unit_Converter()

class Noise(Genesys_Class):
    """
    Class for simulating noise TOD
    Refer to numerical recipes (NR hereafter) in C++, second edition
    """
    def __init__(self, noise_params):
        self.copy_params(noise_params)

    def simulate_noise(self, segment_length):
        n_samples = segment_length * self.params['sampling_rate']
        if self.params['noise_type'] == 'white':
            noise = self.get_white_noise(n_samples)
        if self.params['noise_type'] == '1_over_f':
            noise = self.get_one_over_f_noise(n_samples)
        return noise

    def get_white_noise(self, n_samples):
        noise_tod = np.random.normal(loc=0.0, scale=self.params['white_noise_rms'], size=n_samples)
        return noise_tod

    def get_one_over_f_noise(self, n_samples):
        """
        n_samples is assumed even. This assures padding gives a power of 2
        """
        n_samples_2 = self.get_higher_power_of_2(n_samples)
        pad_length = (n_samples_2 - n_samples) // 2
        fft_freq = fft.rfftfreq(n_samples_2, 1.0/self.params['sampling_rate'])
        noise_profile = self.get_one_over_f_profile(fft_freq)
        white_noise_tod = self.get_white_noise(n_samples_2)
        white_noise_fft = fft.rfft(white_noise_tod)
        one_over_f_noise_fft = white_noise_fft * noise_profile
        one_over_f_noise_tod = fft.irfft(one_over_f_noise_fft)[pad_length:-pad_length]
        return one_over_f_noise_tod

    def get_one_over_f_profile(self, freq):
        freq_mHz = uc.convert_unit(freq, 'frequency', 'Hz', 'mHz')
        noise_profile = np.ones(freq_mHz.size)
        noise_profile[1:] += (self.params['f_knee'] / freq_mHz[1:]) ** self.params['alpha']
        noise_profile[0] = 0.0
        return noise_profile

    def get_higher_power_of_2(self, n_samples):
        h_exp = math.ceil(np.log2(n_samples))
        return 2**h_exp

    #  def get_real_one_sided_psd(self, timestream):
        #  '''
           #  Takes the one sided power spectral density per unit time of the input timestream.
           #  NR : Section 12.0, Eq 12.0.14 for definition of one-sided PSD
           #  Units: timestream_units**2/Hz
           #  NR : Section 13.4 for PSD
           #  The PSD is normalized such that sum(PSD)/total_time = variance
           #  For white noise, mean(PSD)*(largest_f-smallest_f) = variance
        #  '''
        #  n_samples = timestream.size
        #  #NR : Section 12.1, Eq 12.1.5 but only one-sided for FREQ
        #  psd = np.abs(np.fft.rfft(timestream))**2
        #  norm = n_samples**2 / (n_samples / self.config.sampling_rate)
        #  psd /= norm
        #  psd[1:-1] *= 2
        #  return psd
#
#
    #  def get_one_sided_frequencies(self):
        #  n_samples = self.config.t_segment * self.config.sampling_rate
        #  freq = np.fft.rfftfreq(n_samples, 1.0/self.config.sampling_rate)
        #  return freq
#
#
    #  def simulate_timestream_noise_from_parameters(self, do_mc=False):
        #  #white_noise_sigma units : signal_unit*sqrt(s)
        #  #So, to get the SD of each sample point we must use white_noise_sigma*sqrt(sampling_rate) as the scale for the Gaussina RNG
        #  n_samples = self.config.t_segment * self.config.sampling_rate
#
        #  if self.config.noise_type=="white":
            #  noise = np.random.normal(loc=0.0, scale=self.config.white_noise_sigma*np.sqrt(self.config.sampling_rate), size=n_samples)
#
        #  #Needs updating. Action item
        #  if self.config.noise_type=="one_over_f":
            #  if do_mc:
                #  np.random.seed()
            #  else:
                #  np.random.seed(self.config.one_over_f_seed)
            #  freq = get_one_sided_frequencies(n_samples, sampling_rate)
            #  noise = np.random.normal(loc=0.0, scale=self.config.white_noise_sigma*np.sqrt(self.config.sampling_rate), size=n_samples)
            #  noise_fft = np.fft.rfft(noise)*one_over_f_profile(freq, f_knee, alpha)
            #  noise = np.fft.irfft(noise_fft)
#
        #  return noise
#
#
    #  def noise_psd_mc(self, mc_samples):
        #  #noise_variance units : timestream_unit/sqrt(Hz) ....... Ex : uK/sqrt(Hz)
        #  n_samples = self.config.t_segment * self.config.sampling_rate
        #  freq = self.get_one_sided_frequencies()
        #  noise_psd = np.zeros(freq.size)
#
        #  for i in range(mc_samples):
            #  #noise = np.random.normal(loc=0.0, scale=white_noise_sigma*np.sqrt(sampling_rate), size=n_samples)
            #  noise = simulate_timestream_noise_from_parameters(white_noise_sigma, total_time, sampling_rate, noise_type, f_knee, alpha)
            #  noise_psd += get_real_one_sided_psd(noise, sampling_rate)
#
        #  noise_psd /= mc_samples
#
        #  return noise_psd
#
#
    #  def one_over_f_profile(self, freq, f_knee, alpha):
        #  profile = np.ones(freq.size)
        #  valid = freq>0
        #  profile[valid] *= (1 + (f_knee/freq[valid])**(alpha))
#
        #  return profile
#
#
    #  def get_theoretical_noise_spectra(self, white_noise_sigma, total_time, sampling_rate, noise_type, f_knee=None, alpha=None):
        #  n_samples = total_time*sampling_rate
#
        #  if noise_type == "white_noise":
            #  noise_psd = np.full(n_samples/2 + 1, 2*white_noise_sigma**2)
        #  if noise_type == "full":
            #  freq = get_one_sided_frequencies(n_samples, sampling_rate)
            #  noise_psd = np.full(freq.size, 2*white_noise_sigma**2)*one_over_f_profile(freq, f_knee, 2*alpha)
#
        #  return noise_psd
#
#
