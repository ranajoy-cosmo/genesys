import numpy as np
import healpy as hp
import matplotlib.pyplot as plt
from genesys.numerical.unit_conversion import Unit_Converter

uc = Unit_Converter()

#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
#* Physical constants in SI units
#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
c = 2.99792458 * 1e8                # m / s
k_B = 1.38064852 * 1e-23            # m^2 kg / s^2 / K
h = 6.62607004 * 1e-34              # m^2 kg / s
T_CMB_0 = 2.7255                    # K
#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*

def frequency_to_wavelength(nu, unit_in, unit_out, verbose=False):
    nu_Hz = uc.convert_unit(unit_type="frequency", quantity=nu, unit_in=unit_in, unit_out='Hz')
    w_length_m = c / nu_Hz
    w_length = uc.convert_unit(unit_type="length", quantity=w_length_m, unit_in='m', unit_out=unit_out)
    if verbose:
        print(f"{nu} {unit_in} = {w_length} {unit_out}")
    return w_length

def wavelength_to_frequency(w_length, unit_in, unit_out, verbose=False):
    w_length_m = uc.convert_unit(unit_type="length", quantity=w_length, unit_in=unit_in, unit_out='m')
    nu_Hz = c / w_length_m
    nu = uc.convert_unit(unit_type="frequency", quantity=nu_Hz, unit_in='Hz', unit_out=unit_out)
    if verbose:
        print(f"{w_length} {unit_in} = {nu} {unit_out}")
    return nu

#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*

def photon_energy(nu=None, w_length=None, unit_in=None, unit_out=None, verbose=False):
    if nu != None:
        nu_Hz = uc.convert_unit(unit_type="frequency", quantity=nu, unit_in=unit_in, unit_out='Hz')
        E_J = h * nu_Hz
        E = uc.convert_unit(unit_type="energy", quantity=E_J, unit_in="J", unit_out=unit_out)
        if verbose:
            print(f"{nu} {unit_in} -> {E} {unit_out}")
    if w_length != None:
        w_length_m = uc.convert_unit(unit_type="length", quantity=w_length, unit_in=unit_in, unit_out='m')
        E_J = h * c / w_length_m
        E = uc.convert_unit(unit_type="energy", quantity=E_J, unit_in="J", unit_out=unit_out)
        if verbose:
            print("{w_length} {unit_in} -> {E} {unit_out}")
    return E

def get_energy_fraction(nu, T):
    """
    Gives us the ratio h*nu/k_B*T
    nu is in Hz, T in Kelvin
    """
    energy_fraction = h*nu / (k_B*T)
    return energy_fraction

def B_nu_planck(nu, T, unit_in):
    nu_Hz = uc.convert_unit(unit_type="frequency", quantity=nu, unit_in=unit_in, unit_out='Hz')
    x = get_energy_fraction(nu_Hz, T)
    B_nu = (2*h*nu_Hz**3 / c**2) / (np.exp(x) - 1.0)
    return B_nu

def B_nu_RJ(nu, T, unit_in):
    nu_Hz = uc.convert_unit(unit_type="frequency", quantity=nu, unit_in=unit_in, unit_out='Hz')
    B = 2*nu_Hz**2*k_B*T / c**2
    return B

def uK_Planck(nu, T, unit_in):
    nu_Hz = uc.convert_unit(unit_type="frequency", quantity=nu, unit_in=unit_in, unit_out='Hz')
    x = get_energy_fraction(nu_Hz, T)
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
