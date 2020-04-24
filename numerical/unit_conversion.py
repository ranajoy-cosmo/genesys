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
        # The unit types, and the units are defined individually in this section
        self.unit_dict = {}
        self.unit_dict['angle'] = {'radian': 1.0, 'degree': pi/180.0, 'arcmin': pi/180.0/60.0, 'arcsec': pi/180.0/60.0/60.0}
        self.unit_dict['solid_angle'] = {'full_sky': 4*pi, 'steradian': 1.0, 'degree_square': (pi/180.0)**2, 'arcmin_square': (pi/180.0/60.0)**2, 'arcsec_square': (pi/180.0/60.0/60.0)**2}
        self.unit_dict['angular_velocity'] = {'rpm': 1.0, 'rps': 60.0, 'radians/min': 1.0/(2*pi), 'radians/sec': 60.0/(2*pi), 'degree/min': 1.0/360.0, 'degree/sec': 60.0/360.0}
        self.unit_dict['length'] = {'parsec': 30856775814913700, 'light_year': 9460730472580800, 'km': 1000.0, 'm': 1.0, 'cm': 1e-2, 'mm': 1e-3, 'mu_m': 1e-6, 'nm': 1e-9, 'angstrom': 1e-10}
        self.unit_dict['time'] = {'year': 365*24*60*60.0, 'day': 24*60*60.0, 'hour': 60*60.0, 'minute': 60.0, 'sec': 1.0, 'milli-sec': 1e-3, 'micro-sec': 1e-6}
        self.unit_dict['frequency'] = {'mHz': 1.0e-3, 'Hz': 1.0, 'kHz': 1.0e3, 'MHz': 1.0e6, 'GHz': 1.0e9, 'THz': 1.0e12}
        self.unit_dict['energy'] = {'meV': 1.0e-3, 'eV': 1.0, 'keV': 1.0e3, 'MeV': 1.0e6, 'GeV': 1.0e9, 'TeV': 1.0e12, 'PeV': 1.0e15, 'J': 1.60218e-19, 'C': 3.82929e-20}

    def verify_units(self, unit_type, unit_in, unit_out):
        valid_unit_types = self.unit_dict.keys()
        assert unit_type in valid_unit_types, f"{unit_type} not a valid unit_type. List of valid unit types: {list(valid_unit_types)}"
        valid_units = self.unit_dict[unit_type].keys()
        assert unit_in in valid_units, f"{unit_in} not a valid unit for {unit_type}. List of valid units: {list(valid_units)}"
        assert unit_out in valid_units, f"{unit_out} not a valid unit for {unit_type}. List of valid units: {list(valid_units)}"

    def conversion_factor(self, unit_type, unit_in, unit_out, verbose=False):
        self.verify_units(unit_type, unit_in, unit_out)
        unit_dict = self.unit_dict[unit_type]
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
