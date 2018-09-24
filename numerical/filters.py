#!/usr/bin/env python 

import numpy as np
import scipy as sp


def find_larger_closest_power_of_two(n):
    try:
        length = 2**(n-1).bit_length() # closest power of two >= n, only supported in python >= 2.7
    except AttributeError:
        length = 2**(int(np.ceil(np.log2(n))))
    return length

def pad_array(signal):
    length = signal.size
    padded_length = find_larger_closest_power_of_two(length) 
    pad_length = (padded_length - length)/2
    return np.lib.pad(signal, (pad_length, pad_length), 'reflect'), pad_length

def get_fft(signal, sample_rate=200.0):
    N = signal.size
    signal_padded = pad_array(signal)
    FREQ = np.fft.rfftfreq(N, 1.0/sample_rate)
    FFT = np.fft.rfft(signal)
    return FREQ, FFT

def butter_gain_low_pass(freq, filter_order, critical_freq):
    critical_freq = float(critical_freq)
    #freq = np.fft.fftfreq(ts_length, 1.0/sample_rate)
    gain = 1.0/np.sqrt(1.0+(freq/critical_freq)**(2*filter_order))
    return gain

def filter_butter(signal, sample_rate, freq_cutoff):
    signal_padded, pad_length = pad_array(signal)
    freq, signal_fft = get_fft(signal_padded, sample_rate)
    fft_filtered = signal_fft*butter_gain_low_pass(freq, 10, freq_cutoff)
    signal_filtered = np.fft.irfft(fft_filtered)[pad_length:pad_length+signal.size]
    return signal_filtered
