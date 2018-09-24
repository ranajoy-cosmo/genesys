# spectral analysis functions:
# originally written by Hannes.
# updated by Kate, October 2013.

import numpy as np
import scipy as sp
import scipy.signal
import matplotlib.pyplot as plt
import math
import copy

sample_rate = 200.0


###########################################################################################
# the core functions for psd, smoothing, etc.


def clean(data, cutoff=4, verbose=1):
    '''function which removes >= cutoff*sigma events and replaces
    with the mean value. data should be a list. cutoff defines
    that points above cutoff*standard deviation are replaced with
    the mean. Note this function aint going to be very good if the
    data has large 1/f noise.
    '''

    N = len(data)
    clean_data = copy.copy(data)
    mu = sp.mean(np.array(data))
    sigma = sp.std(np.array(data))
    num_cleanedpoints = 0
    for i in range(N):
        if abs(clean_data[i] - mu) >= cutoff * sigma:
            clean_data[i: i+1] = [mu]
            num_cleanedpoints = num_cleanedpoints + 1
        else:
            pass
    if verbose == 1:
        print('mean = %.3f\nSTD = %.3f\t%d' % (mu, sigma, num_cleanedpoints))
    return clean_data


def window(data, window='hanning', renormalize=1):
    ''' Windows the data with hanning window
        If renormalize is set to 1, the data is multiplied by
        a factor which makes the mean squared amplitude of the
        windowed data equal to the mean squared amplitude of the
        original data.  This is important for normalizing a PSD
        correctly.
    '''
    N = len(data)
    rescale = (8/3.)**0.5  # normalize the power for a hanning window
    if renormalize == 1:
        return data*np.hanning(N)*rescale
    else:
        return data*np.hanning(N)


def fft_coefficients_ranajoy(timestream, samplerate=200.0):
    ''' Calculates the fft amplitude coefficients for the input
        timestream.  timestream is a list or an array.
        data[0] = frequency
        data [1] = amplitude coefficients
        '''
    N = len(timestream)
    FFT = np.abs(np.fft.fft(timestream)[:N/2]) / N
    FFT[1:-1] *= 2
    FREQ = np.fft.fftfreq(N, 1/samplerate)[:N/2]
    return np.vstack((FREQ, FFT))


def fft_coefficients(timestream, samplerate=200.0):
    ''' Calculates the fft amplitude coefficients for the input
        timestream.  timestream is a list or an array.  The output
        data is a nested list.

        data[0] = frequency
        data [1] = amplitude coefficients
        '''
    # create data list and find length of samples
    data = []
    n = len(timestream)
    # create frequency array
    delta = 1 / samplerate
    deltaf = samplerate / n
    fc = samplerate / 2.0
    f = np.array(range(0, n / 2))
    frequency = fc * f / (n / 2)
    frequency = list(frequency)
    data.append(frequency)
    # take FFT and normalize so that the output is in units of amplitude
    # negative frequencies are thrown away and a factor of 2 is added to
    # compensate.
    FFT = sp.fft(timestream, n) / n
    FFT = np.abs(resize(FFT, (n / 2,)))
    FFT[0] = FFT[0]
    FFT[1: n / 2 - 1] = 2 * FFT[1: n / 2 - 1]
    data.append(sp.real(FFT))
    return data

def psd_func_real(timestream, sample_rate=sample_rates.bolo, truncate=False):
    '''
       Takes the power spectral density per unit time of the input timestream.
       Units: timestream_units**2/Hz
       The PSD is normalized such that sum(PSD)/total_time = variance
       For white noise, mean(PSD)*(largest_f-smallest_f) = variance
       If truncate=True, the array is truncated to the closest power of two.
       WARNING: this is a fast implentation that works for REAL tiemstream only
       Output:
        [frequencies, PSD]
    '''
    ts = timestream
    N = len(timestream)
    if truncate:
        if (N & N-1):
            a = int(math.log(N, 2))
            N = 2**a
            ts = timestream[:N]
    #rfftfreq is not supported at Minnesota, using the slightly slower fftfreq when rfftfreq is not available
    try:
        FREQ = np.abs(np.fft.rfftfreq(N, 1.0/sample_rate))
    except AttributeError:
        FREQ = np.abs(np.fft.fftfreq(N, 1.0/sample_rate)[:N/2+1])
    PSD = np.abs(np.fft.rfft(ts))**2
    norm = N**2 / (N / sample_rate)
    PSD /= norm
    PSD[1:-1] *= 2
    return np.vstack((FREQ, PSD))


def psd_func(timestream, samplerate=sample_rates.bolo):
    '''Takes the power spectral density of the input timestream.
       The format of timestream is a list.  It CANNOT be a nested list.
       The output data is:
       data[0] is the frequency list
       data[1] is the psd in units of sqrt(power)/sqrt(Hz)
    '''

    # take the largest power of 2 data points in array for psd
    # this makes the function quicker if len(array)= power of 2
    data = []  # data is a list
    # make sure timestream has an even number of samples. fft is faster
    n = len(timestream)
    # make frequency column
    delta = 1 / samplerate
    deltaf = samplerate / n
    fc = samplerate / 2.0  # critical frequency
    f = np.array(range(0, n / 2))
    frequency = fc * f / (n / 2)  # frequency bins
    frequency = list(frequency)
    data.append(frequency)
    # do the fft; throw away negative frequencies
    PSD = sp.fft(timestream, n)
    PSD = np.resize(PSD, (n / 2,))
    # Do normalization from Numerical Recipes in C pg 551
    norm = n * n * deltaf
    PSD[0] = PSD[0] * sp.conj(PSD[0]) / norm
    PSD[1: n / 2-1] = 2 * PSD[1:n / 2 - 1] * sp.conj(PSD[1: n / 2 - 1]) / norm
    PSD[n / 2 - 1] = PSD[n / 2 - 1] * sp.conj(PSD[n / 2 - 1]) / norm
    # dividing by n**2 makes the sum of the coefficients equal to the
    # mean squared amplitude of the time stream.  multiplying by n*delta =
    # delta f converts to PSD (ie per Hz)
    PSD = sp.real(PSD)
    PSD = list((PSD) ** 0.5)  # convert to sqrt(power_rms)/sqrt(Hz) ie Vrms/rtHz
    data.append(PSD)
    return data


def cospectrum(data, samplerate=25.e6 / 2 ** 16):
    '''computes the cospectrum of the input data[0] and data[1]
       data format is an embedded list.  Returns a dictionary in
       the DfMUX python library data format
        '''
    dataout = {}
    n = len(data[0])
    # make frequency column
    delta = 1 / samplerate
    deltaf = samplerate / n
    fc = samplerate / 2.0  # critical frequency
    f = np.array(range(0, n / 2))
    frequency = fc * f / (n / 2)  # frequency bins
    frequency = list(frequency)
    dataout['xdata'] = frequency
    dataout['xtitle'] = 'frequency'
    # take fft of input time stream / remove negative frequencies
    h11 = sp.fft(data[0], n)
    h11 = np.resize(h11, (n / 2,))
    h22 = sp.fft(data[1], n)
    h22 = np.resize(h22, (n / 2,))
    # compute auto spectra and cross spectra unnormalized
    #norm = n*n*deltaf
    h21 = sp.conj(h11) * h22
    h11[0] = h11[0] * sp.conj(h11[0])
    h11[1: n / 2 - 1] = 2 * h11[1:n / 2 - 1] * sp.conj(h11[1:n/2-1])
    h11[n / 2 - 1] = h11[n / 2 - 1]*sp.conj(h11[n / 2 - 1])
    h11 = sp.real(h11)
    h22[0] = h22[0] * sp.conj(h22[0])
    h22[1: n / 2 - 1] = 2 * h22[1: n / 2 - 1] * sp.conj(h22[1: n / 2 - 1])
    h22[n / 2 - 1] = h22[n / 2 - 1] * sp.conj(h22[n / 2 - 1])
    h22 = sp.real(h22)
    coherency = np.abs(h21)/(h11 * h22) ** 0.5
    coherency = sp.real(coherency)
    dataout['y1data'] = list(h11)
    dataout['y2data'] = list(h22)
    dataout['y3data'] = list(coherency)
    dataout['ylabel'] = 'data out'
    dataout['y1label'] = 'real'
    return dataout


def test_cospectrum():
    x = []
    i = 0
    n = 2 ** 14
    noise1 = sp.random.normal(0, 2, n)
    noise2 = sp.random.normal(0, 4, n)
    timestream = []
    timestream.append(noise1)
    timestream.append(noise2)
    CS = cospectrum(timestream)
    return CS


def test_psd_normalization():
    ''' This function tests the normalization of function psd. Mock data is
        one second of normal, mean zero, std = 2 data sampled at
        1kHz.  Since this is white noise, the white noise level of the PSD times
        the root of the bandwidth should give the rms amplitude of the
        data (in this case rt(2)).

        The normalization for a hanning window is also tested.  Windowing
        the data removes power from the time stream.  The data must be
        recalibrated in order to recover the best estimate of the white
        noise level.  For a hanning window the time stream must be multipled by
        root(8/3) before the PSD is taken.
        '''

    # make fake data, window, window and rescale
    x = sp.random.normal(0, 2, 10000)
    wrx = window(x, 'hanning', 1)
    ms_x = sp.mean(x ** 2)
    ms_wrx = sp.mean(np.array(wrx) ** 2)
    ratio = ms_x / ms_wrx
    print ('MSA of timestream = %.4f\t\nMSA of windowed timestream = %.4f\nratio = %.4f' % (ms_x, ms_wrx, ratio))
    # take PSDs
    x_psd = psd(x, 381.47)
    wrx_psd = psd(wrx, 381.47)
    pylab.subplot(2, 1, 1)
    pylab.title('Test psd normalization')
    pylab.xlabel('Sample')
    pylab.ylabel('Cnts')
    pylab.plot(x, 'bo', wrx, 'ro')
    pylab.subplot(2, 1, 2)
    pylab.title('PSD')
    pylab.xlabel('Frequency [Hz]')
    pylab.ylabel('Cnts/rtHz')
    pylab.loglog(x_psd[0], x_psd[1], 'b-', wrx_psd[0], wrx_psd[1], 'r-')
    pylab.show()
    #return PSD


def rebin(a, *args):
    #taken from sp cookbook.  It mimics REBIN in IDL.
    '''rebin ndarray data into a smaller ndarray of the same rank whose dimensions
    are factors of the original dimensions. eg. An array with 6 columns and 4 rows
    can be reduced to have 6,3,2 or 1 columns and 4,2 or 1 rows.
    example usages:
    >>> a=rand(6,4); b=rebin(a,3,2)
    >>> a=rand(6); b=rebin(a,2)
    '''
    if type(a) == list:
        a = np.array(a)
    else:
        pass
    shape = a.shape
    lenShape = len(shape)
    factor = np.asarray(shape) / np.asarray(args)
    evList = ['a.reshape('] + \
             ['args[%d],factor[%d],' % (i, i) for i in range(lenShape)] + \
             [')'] + ['.mean(%d)' % (i + 1) for i in range(lenShape)]
    #print ''.join(evList)
    return eval(''.join(evList))


def rebin_xlog(ARRAY, ppi=128):
    '''rebins ARRAY over the xaxis and rebins on a logscale.
    '''
    rebin_data = []
    for column in range(len(ARRAY)):
        r_data = []
        # leave 1st chunk of data alone
        for i in range(0, ppi):
            r_data.append(ARRAY[column][i])
        # average next chunks of data logrithmically
        k = ppi
        while k*2 <= len(ARRAY[column]):
            x = rebin(ARRAY[column][k: k * 2], ppi)
            x = list(x)
            r_data = r_data + x
            k = k * 2
        rebin_data.append(r_data)
    return rebin_data


def integratePSD(FREQ, PSD, peak=0, freq_center=24., bandwidth=1, verbose=0):
    ''' module to calculate the rms signal in a given
    bandwidth of a power spectral density.  The input data
    is an embedded list with the first list the frequency data.
    The second list is the PSD in units of rms/rtHz.
    '''

    # find the binwidth of one sample and the number of bins in the bandwidth
    deltaf = FREQ[2] - FREQ[1]
    if peak == 1:
        bandwidth = 0
    numpoints = int(bandwidth / deltaf)
    if numpoints == 0:
        numpoints = 1

    # find the index of the center freuquency
    for f in FREQ:
        if freq_center-deltaf / 2 < f < freq_center + deltaf / 2:
            fc_index = FREQ.index(f)
            center = f
        if f > freq_center + 2 * deltaf:
            break
    if peak == 1:
        return PSD[fc_index]
    # Calculate the average signal in the band
    band = []
    if numpoints == 1:
        band.append(PSD[fc_index])
    else:
        band = PSD[fc_index-numpoints / 2: fc_index + numpoints / 2]
    average_signal = 0
    for i in band:
        average_signal = average_signal + i ** 2
    average_signal = (average_signal*numpoints*deltaf) ** 0.5
    if verbose == 1:
        print('center_frequency = %.4f\ndeltaf = %3f\nnumpoints = %d' % (center, deltaf, numpoints))
    return average_signal

#################################################################################################
# Functions build on the above core designed specifically for dfmux library
#################################################################################################


def dfmux_psd(data, clean_data=1, detrend_data=1, window_data=1, frequency_smooth=1, plotthedata=1):
    '''Designed for use with the dfmux library.  Takes the input
       bolometer data of format:
       data['Ch1'] = list
       data['Ch2'] = list , etc
       and returns a dictionary of power spectral densities.
       Options: (0 is off, 1 is on)

           clean_data: removes 4 sigma outliers
           window_data: windows the data with a hanning window
           detrend_data: 1 removes the DC and 2 the linear term from the data
           frequency_smooth = rebins the frequency bins to give equal spaced
                              points on a log plot.
           plotthedata: plots the timestream and psd to screen

       The output format matches daq.write_alg_outputtofile
        '''
    # the returned data is a dictionary conforming to format applicable
    # to write_alg_outputtofile
    psddata = {}
    psddata['title'] = 'PSD of bolometers'
    psddata['xlabel'] = 'frequency [Hz]'
    psddata['ylabel'] = 'noise [ADC/rtHz]'
    datakeys = data.keys()
    index = 1
    for key in datakeys:
        x = np.array(data.get(key))
        # the clean up process
        if clean_data == 1:
            x = clean(x, 4, 0)
        if detrend_data == 1:
            x = sp.signal.detrend(x, type == 'constant')
        elif detrend_data == 2:
            x = detrend_linear(x, type == 'linear')
        if window_data == 1:
            x = window(x,  'hanning', 1)
        x = psd(x,  381.47)
        # smoothing in frequency domain
        if frequency_smooth == 1:
            x = rebin_xlog(x, ppi=16)
        keyname = 'y%1ddata' % index
        psddata[keyname] = x[1]
        keyname = 'y%1dlabel' % index
        psddata[keyname] = '%s' % key
        index = index + 1
    psddata['xdata'] = x[0]
    # plot the data
    index = 1
    if plotthedata == 1:
        for channel in datakeys:
            pylab.subplot(111)
            pylab.subplot(2, 1, 1)
            pylab.title('Timestream of channel %s' % channel)
            pylab.xlabel('Sample')
            pylab.ylabel('ADC counts')
            pylab.plot(data[channel], 'bo')
            pylab.subplot(2, 1, 2)
            pylab.title('PSD')
            ydata = 'y' + str(index) + 'data'
            pylab.loglog(psddata['xdata'], psddata[ydata], 'b-')
            pylab.xlabel("Frequency [Hz]")
            pylab.ylabel("Noise [ADC/rtHz]")
            print('Type any number, then enter to see the next timestream/psd')
            pylab.show()
            input()
            index = index + 1
    else:
        pass
    return psddata


def dfmux_Etau_amp(psd, freq_center, bandwidth, peak=0):
    ''' Module to find the amplitude of the frequency component corresponding
    to a voltage bias signal.
    '''
    # all the needed dictionaries
    data_in = psd.copy()  # make a copy otherwise original is tampered with
    data_out = {}  # create the outputdata
    # Integrate spectrums around freq_center
    N = len(data_in)
    for i in range(1, N + 1):
        datatag = 'y' + str(i) + 'data'
        if datatag in data_in.keys():
            labeltag = data_in['y' + str(i) + 'label']
            data_out[labeltag] = integratePSD(data_in['xdata'], data_in[datatag], peak, freq_center, bandwidth)
    return data_out


def dfmux_fft_coefficients(data, plotthedata=1):
    '''Takes the dictionary from sampling the DfMUX board and returns
       the amplitude coefficients of the timestreams
    #assumes sampling rate!
    '''
        # the returned data is a dictionary conforming to format applicable
    # to write_alg_outputtofile
    fftdata = {}
    fftdata['title'] = 'FFT amplitude coefficients of bolometer timestreams'
    fftdata['xlabel'] = 'frequency [Hz]'
    fftdata['ylabel'] = 'Amplitude [ADC]'
    index = 1
    for key in data.keys():
        x = data.get(key)
        x = fft_coefficients(x, 190.73)
        keyname = 'y%1ddata' % index
        fftdata[keyname] = x[1]
        keyname = 'y%1dlabel' % index
        fftdata[keyname] = '%s' % key
        index = index + 1
    fftdata['xdata'] = x[0]
    # plot the data
    index = 1
    if plotthedata == 1:
        for channel in data.keys():
            pylab.subplot(111)
            pylab.subplot(2, 1, 1)
            pylab.title('Timestream of channel %s' % channel)
            pylab.xlabel('Sample')
            pylab.ylabel('ADC counts')
            pylab.plot(data[channel], 'bo')
            pylab.subplot(2, 1, 2)
            pylab.title('FFT coefficients')
            ydata = 'y' + str(index) + 'data'
            pylab.loglog(fftdata['xdata'], fftdata[ydata], 'b-')
            pylab.xlabel("Frequency [Hz]")
            pylab.ylabel("Amplitude [ADC]")
            print('Type any number, then enter to see the next timestream/psd')
            pylab.show()
            input()
            index = index + 1
    else:
        pass
    return fftdata


def psd_readin(inputfile, columns=[0], window=1, smooth=1):
    '''Computes/Plots/writes to file the power spectra of the timestreams
    given in the columns of filename.  The smooth option runs the data
    through rebin_xlog.
    '''

    # define the sample rate (rate for FIR stage 5)
    samplerate = 381.47
    # read in the data. data is a list. it may be nested
    data = read_in_data_from_file(inputfile, columns)
    # window the data if window = 1
    if window == 1:
        # control cases with 1 or many timestreams
        num_datacolumns = len(data)
        if num_datacolumns > 20:  # this is my assumption for 1 timestream
            num_points = num_datacolumns
            data = list(np.hanning(num_points) * np.array(data))
        elif num_datacolumns > 1:
            # !!this option currently does not work if used with the smooth command
            x = []
            num_points = len(data[0])
            for i in data:
                ts = np.hanning(num_points) * np.array(data)
                x.append(list(ts))
            data = x
    else:
        pass
    # take the power spectral density
    PSD = psd_original(data, samplerate)
    # smooth the frequency domain data if smooth = 1
    if smooth == 1:
        PSD = rebin_xlog(PSD, ppi=8)
    else:
        pass
    # plot the data
    for i in range(1, len(PSD)):
        pylab.subplot(4, 2, i)
        pylab.loglog(PSD[0], PSD[i], 'k.', PSD[0], PSD[i], 'b-')
        pylab.xlabel("Freq [Hz]")
        pylab.ylabel("Noise")
    show()
    # write the data to a file with filename 'inputfile.psd'
    write(PSD, inputfile + '.psd')


def calculate_psd(timestream, rate=None):
    if rate is None:
        rate = sample_rates.bolo
    psd = sp.signal.detrend(timestream)
    psd = window(psd, 'hanning', 1)
    psd = psd_func(psd, rate)
    return psd


def calculate_psd_clean(timestream, rate=None, truncate=False):
    if rate is None:
        rate = sample_rates.bolo
    N = len(timestream)
    if truncate:
        if (N & N-1):
            a = int(math.log(N, 2))
            N = 2**a
    ts = sp.signal.detrend(timestream[:N])
    ts = window(ts, 'hanning', 1)
    psd = psd_func_real(ts, rate)
    return psd
