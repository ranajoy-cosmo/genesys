"""
#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
# Unit Conversion Routine
#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
This module contains a single routine that converts between different units.
If you feel a certain new unit type or a new unit needs to be included, please contact the developer or hack it in yourself.
The syntax for calling the routine is:
    convert_unit(unit_type, quantity, unit_in, unit_out)
    For example:
    convert_unit('angle', 0.1, unit_in='radians', unit_out='degrees')
All convrsion routines work with scalars as well as with numpy arrays. The return type is the same as the input type.
#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
"""

import numpy as np
pi = np.pi

class Unit_Converter:
    def __init__(self):
        # THE UNIT TYPES, AND THE UNITS ARE DEFINED INDIVIDUALLY IN THIS SECTION
        self.all_unit_dicts = {
                'angle': {'radian' : 1.0, 'degree' : pi/180.0, 'arcmin' : pi/180.0/60.0, 'arcsec' : pi/180.0/60.0/60.0},
                'solid_angle': {'full_sky' : 4*pi, 'steradian' : 1.0, 'degree_square' : (pi/180.0)**2, 'arcmin_square' : (pi/180.0/60.0)**2, 'arcsec_square' : (pi/180.0/60.0/60.0)**2},
                'angular_velocity': {'rpm' : 1.0/60.0, 'rps' : 1.0, 'radians/min' : 2*pi/60.0, 'radians/sec' : 2*pi},
                'length': {'parsec' : 30856775814913700, 'light year' : 9460730472580800, 'km': 1000.0, 'm': 1.0, 'cm' : 1e-2, 'mm' : 1e-3, 'mu_m' : 1e-6, 'nm' : 1e-9, 'angstrom' : 1e-10},
                'time': {'siderial year' : 365.25*24*60*60.0, 'day' : 24*60*60.0, 'hour' : 60*60.0, 'minute' : 60.0, 'second' : 1.0, 'milli-second' : 1e-3, 'micro-second' : 1e-6},
                'frequency': {'mHz': 1.0e-3, 'Hz': 1.0, 'kHz' : 1.0e3, 'MHz': 1.0e6, 'GHz': 1.0e9, 'THz': 1.0e12},
                'energy': {'meV': 1.0e-3, 'eV': 1.0, 'keV': 1.0e3, 'MeV': 1.0e6, 'GeV': 1.0e9, 'TeV': 1.0e12, 'PeV': 1.0e15, 'J': 1.60218e-19, 'C': 3.82929e-20}}

    def conversion_factor(self, unit_type, unit_in, unit_out, verbose=False):
        unit_dict = self.all_unit_dicts[unit_type]
        conv_factor = unit_dict[unit_in] / unit_dict[unit_out]
        if verbose:
            print("{} to {} conversion factor: {}".format(unit_in, unit_out, conv_factor))
        return conv_factor

    def convert_unit(self, quantity, unit_type, unit_in, unit_out, verbose=False):
        conv_factor = self.conversion_factor(unit_type, unit_in, unit_out)
        converted_quantity = quantity * conv_factor
        if verbose:
            print("{} {} = {} {}".format(quantity, unit_in, converted_quantity, unit_out))
        return converted_quantity
