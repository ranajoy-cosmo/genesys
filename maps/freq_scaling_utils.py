import numpy as np
import healpy as hp
import matplotlib.pyplot as plt

#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
#* Physical constants in SI units
#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
c = 2.99792458 * 1e8                # m / s
k_B = 1.38064852 * 1e-23            # m^2 kg / s^2 / K
h = 6.62607004 * 1e-34              # m^2 kg / s
T_CMB = 2.7255                      # K
#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*

def convert_freq_units(nu_in, in_units, out_units):
    freq_unit_dict = {'mHz' : 1.0e-3, 'Hz' : 1.0, 'MHz' : 1.0e6, 'GHz' : 1.0e9, 'THz' : 1.0e12}
    conversion_factor = freq_unit_dict[in_units] / freq_unit_dict[out_units]
    return nu_in * conversion_factor

def convert_length_units(w_length_in, in_units, out_units):
    length_unit_dict = {'nm' : 1.0e-9, 'um' : 1.0e-6, 'mm' : 1.0e-3, 'm' : 1.0, 'Km' : 1.0e3}
    conversion_factor = length_unit_dict[in_units] / length_unit_dict[out_units]
    return w_length_in * conversion_factor

def frequency_to_wavelength(nu, in_unit, out_unit):
    nu_Hz = convert_freq_units(nu, in_unit, 'Hz')
    w_length_m = c / nu_Hz
    w_length = convert_length_units(w_length_m, 'm', out_unit)
    return w_length

def wavelength_to_frequency(w_length, in_unit, out_unit):
    w_length_m = convert_length_units(w_length, in_unit, 'm')
    nu_Hz = c / w_length_m
    nu = convert_freq_units(nu_Hz, 'Hz', out_unit)
    return nu

def B_nu_planck(nu, T, in_units):
    nu_Hz = convert_freq_units(nu, in_units, 'Hz')
    x = h*nu_Hz / (k_B*T)
    B = (2*h*nu_Hz**3 / c**2) / (np.exp(x) - 1.0)
    return B

def B_nu_RJ(nu, T, in_units):
    nu_Hz = convert_freq_units(nu, in_units, 'Hz')
    B = 2*nu_Hz**2*k_B*T / c**2
    return B

def uK_Planck(nu, T, in_units):
    nu_Hz = convert_freq_units(nu, in_units, 'Hz')
    x = h*nu_Hz / (k_B*T)
    B = (2*k_B*nu_Hz**2 / c**2) * np.exp(x) * x**2 / (np.exp(x) - 1)**2
    return B

def uK_RJ(nu, T, in_units):
    nu_Hz = convert_freq_units(nu, in_units, 'Hz')
    B = 2*k_B*nu_Hz**2 / c**2
    return B

def uk_Planck_to_uK_CMB(I_nu, nu, T, in_Units):
    conversion_factor = 1.0 / uK_Planck(nu, T, in_units)
    I_nu_uK_CMB = I_nu * conversion_factor
    return I_nu_uK_CMB

def uk_Planck_to_uK_RJ(I_nu, nu, T, in_units):
    conversion_factor = uK_RJ(nu, T, in_units) / uK_Planck(nu, T, in_units)
    I_nu_uK_CMB = I_nu * conversion_factor
    return I_nu_uK_CMB

#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
#* Scaling functions
#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
def ukRJ_to_ukCMB(nu, in_units='GHz'):
    nu_Hz = convert_freq_units(nu, in_units, 'Hz')
    x = h*nu_Hz / (k_B*T_CMB)
    factor = (np.exp(x) - 1.0)**2 / (x**2 * np.exp(x))
    return factor

def scale_dust(nu, nu_0, beta_dust, T_dust, in_units='GHz'):
    nu_Hz = convert_freq_units(nu, in_units, 'Hz')
    nu_0_Hz = convert_freq_units(nu_0, in_units, 'Hz')
    x = h*nu_Hz/(k_B*T_dust)
    x_0 = h*nu_0_Hz/(k_B*T_dust)
    return (nu_Hz/nu_0_Hz)**(beta_dust+1) * (np.exp(x_0) -1)/(np.exp(x) -1)

def scale_sync(nu, nu_0, beta_sync, in_units='GHz'):
    nu_Hz = convert_freq_units(nu, in_units, 'Hz')
    nu_0_Hz = convert_freq_units(nu_0, in_units, 'Hz')
    return (nu_Hz / nu_0_Hz)**beta_sync
#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*

#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
#* Testing scripts
#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*

def plot_spectrum_lambda(T_bb_list, spec_type):
    nu = np.logspace(-1,3,101)
    w_length = frequency_to_wavelength(nu, 'GHz', 'mm')

    spectra_dict = {'B_nu_planck' : B_nu_planck, 'B_nu_RJ' : B_nu_RJ, 'uK_Planck' : uK_Planck, 'uK_RJ' : uK_RJ}
    
    for T_bb in T_bb_list:
        I_nu = spectra_dict[spec_type](nu, T_bb, 'GHz')
        plt.loglog(w_length, I_nu, label=str(T_bb))

    plt.legend()

def plot_spectrum_freq(T_bb_list, spec_type):
    nu = np.logspace(-1,3,101)

    spectra_dict = {'B_nu_planck' : B_nu_planck, 'B_nu_RJ' : B_nu_RJ, 'uK_Planck' : uK_Planck, 'uK_RJ' : uK_RJ}
    
    for T_bb in T_bb_list:
        I_nu = spectra_dict[spec_type](nu, T_bb, 'GHz')
        plt.loglog(nu, I_nu, label=str(T_bb))

    plt.legend()
