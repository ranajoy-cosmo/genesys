#! /usr/bin/env python

import numpy as np
import math
from genesys.utilities import Genesys_Class

"""
Refer to Numerical Recipes (NR hereafter) in C++, second edition
"""

class Noise(Genesys_Class):
    def __init__(self, config):
        self.config = Generic_Class()
        self.config.__dict__.update(config.__dict__)

    def initialise_for_segment(self, segment):
        # Double check this with the one in pointing.py
        delta_t = 1.0 / self.config.sampling_rate
        t_start = segment[1]
        t_stop = segment[2]
        t_duration = t_stop - t_start
        self.n_samples = int(math.ceil(t_duration / delta_t))

    def simulate_ts_noise(self):
        if self.config.noise_type == "white":
            noise = np.random.normal(loc=0.0, scale=self.config.white_noise_sigma, size=self.n_samples)

        return noise

    def get_real_one_sided_psd(self, timestream):
        '''
           Takes the one sided power spectral density per unit time of the input timestream.
           NR : Section 12.0, Eq 12.0.14 for definition of one-sided PSD
           Units: timestream_units**2/Hz
           NR : Section 13.4 for PSD
           The PSD is normalized such that sum(PSD)/total_time = variance
           For white noise, mean(PSD)*(largest_f-smallest_f) = variance
        '''
        n_samples = timestream.size
        #NR : Section 12.1, Eq 12.1.5 but only one-sided for FREQ
        psd = np.abs(np.fft.rfft(timestream))**2
        norm = n_samples**2 / (n_samples / self.config.sampling_rate)
        psd /= norm
        psd[1:-1] *= 2
        return psd 


    def get_one_sided_frequencies(self):
        n_samples = self.config.t_segment * self.config.sampling_rate
        freq = np.fft.rfftfreq(n_samples, 1.0/self.config.sampling_rate)
        return freq


    def simulate_timestream_noise_from_parameters(self, do_mc=False):
        #white_noise_sigma units : signal_unit*sqrt(s)
        #So, to get the SD of each sample point we must use white_noise_sigma*sqrt(sampling_rate) as the scale for the Gaussina RNG
        n_samples = self.config.t_segment * self.config.sampling_rate

        if self.config.noise_type=="white":
            noise = np.random.normal(loc=0.0, scale=self.config.white_noise_sigma*np.sqrt(self.config.sampling_rate), size=n_samples)

        #Needs updating. Action item
        if self.config.noise_type=="one_over_f":
            if do_mc:
                np.random.seed()
            else:
                np.random.seed(self.config.one_over_f_seed)
            freq = get_one_sided_frequencies(n_samples, sampling_rate)
            noise = np.random.normal(loc=0.0, scale=self.config.white_noise_sigma*np.sqrt(self.config.sampling_rate), size=n_samples)
            noise_fft = np.fft.rfft(noise)*one_over_f_profile(freq, f_knee, alpha)
            noise = np.fft.irfft(noise_fft)

        return noise


    def noise_psd_mc(self, mc_samples):
        #noise_variance units : timestream_unit/sqrt(Hz) ....... Ex : uK/sqrt(Hz)
        n_samples = self.config.t_segment * self.config.sampling_rate
        freq = self.get_one_sided_frequencies()
        noise_psd = np.zeros(freq.size)

        for i in range(mc_samples):
            #noise = np.random.normal(loc=0.0, scale=white_noise_sigma*np.sqrt(sampling_rate), size=n_samples)
            noise = simulate_timestream_noise_from_parameters(white_noise_sigma, total_time, sampling_rate, noise_type, f_knee, alpha)
            noise_psd += get_real_one_sided_psd(noise, sampling_rate)

        noise_psd /= mc_samples

        return noise_psd


    def one_over_f_profile(self, freq, f_knee, alpha):
        profile = np.ones(freq.size)
        valid = freq>0
        profile[valid] *= (1 + (f_knee/freq[valid])**(alpha))

        return profile


    def get_theoretical_noise_spectra(self, white_noise_sigma, total_time, sampling_rate, noise_type, f_knee=None, alpha=None):
        n_samples = total_time*sampling_rate

        if noise_type == "white_noise":
            noise_psd = np.full(n_samples/2 + 1, 2*white_noise_sigma**2)
        if noise_type == "full":
            freq = get_one_sided_frequencies(n_samples, sampling_rate)
            noise_psd = np.full(freq.size, 2*white_noise_sigma**2)*one_over_f_profile(freq, f_knee, 2*alpha)

        return noise_psd


