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
        one_over_f_noise_tod -= np.mean(one_over_f_noise_tod)
        return one_over_f_noise_tod

    def get_one_over_f_profile(self, freq):
        freq_mHz = uc.convert_unit(freq, 'frequency', 'Hz', 'mHz')
        noise_profile = np.ones(freq_mHz.size)
        noise_profile[1:] += (self.params['f_knee'] / freq_mHz[1:]) ** self.params['noise_alpha']
        noise_profile[0] = 0.0
        return noise_profile

    def get_higher_power_of_2(self, n_samples):
        h_exp = math.ceil(np.log2(n_samples))
        return 2**h_exp
