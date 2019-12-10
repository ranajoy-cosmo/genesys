"""
#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
# Unit Conversion Routine
#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
This module contains a single routine that converts between different units.
If a unit type does not exist it is appropriately prompted and an exception raised.
If a certain unit does not exist it is appropriately prompted and an exception raised.
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

#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
# The unit types, and the units are defined individually in this section
unit_type_list = ['angle', 'solid_angle', 'length', 'frequency', 'time', 'energy']

all_unit_dicts = dict.fromkeys(unit_type_list)

all_unit_dicts['angle'] = {'radian' : 1.0, 'degree' : pi/180.0, 'arcmin' : pi/180.0/60.0, 'arcsec' : pi/180.0/60.0/60.0}
all_unit_dicts['solid_angle'] = {'full_sky' : 4*pi, 'steradian' : 1.0, 'degree_square' : (pi/180.0)**2, 'arcmin_square' : (pi/180.0/60.0)**2, 'arcsec_square' : (pi/180.0/60.0/60.0)**2}
all_unit_dicts['length'] = {'parsec' : 30856775814913700, 'light year' : 9460730472580800, 'km': 1000.0, 'm': 1.0, 'cm' : 1e-2, 'mm' : 1e-3, 'mu_m' : 1e-6, 'nm' : 1e-9, 'angstrom' : 1e-10}
all_unit_dicts['time'] = {'siderial year' : 365.25*24*60*60.0, 'day' : 24*60*60.0, 'hour' : 60*60.0, 'minute' : 60.0, 'second' : 1.0}
all_unit_dicts['frequency'] = {'mHz': 1.0e-3, 'Hz': 1.0, 'kHz' : 1.0e3, 'MHz': 1.0e6, 'GHz': 1.0e9, 'THz': 1.0e12}
all_unit_dicts['energy'] = {'meV': 1.0e-3, 'eV': 1.0, 'keV': 1.0e3, 'MeV': 1.0e6, 'GeV': 1.0e9, 'TeV': 1.0e12, 'PeV': 1.0e15, 'J': 1.60218e-19, 'C': 3.82929e-20}
#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*

#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
# The general unit conversion routine
#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*

def validate_unit_type(unit_type):
    assert (unit_type in unit_type_list), "Unit type \"{}\" does not exist. Please contact the developer to include it or hack it in yourself.".format(unit_type)
    return

def validate_unit(unit_type, unit):
    unit_dict = all_unit_dicts[unit_type]
    unit_list = list(unit_dict.keys())
    assert (unit in unit_list), "Unit \"{}\" does not exist. Please enter from among the following units\n{}.\nPlease contact the developer to include the required unit or hack it in yourself.".format(unit, unit_list)

def show_valid_unit_types():
    print("Valid unit types : {}".format(unit_type_list))

def show_valid_units(unit_type):
    validate_unit_type(unit_type)
    unit_dict = all_unit_dicts[unit_type]
    unit_list = list(unit_dict.keys())
    print("Valid units : {}".format(unit_list))

def conversion_factor(unit_type, unit_in, unit_out, verbose=False):
    validate_unit_type(unit_type)
    unit_dict = all_unit_dicts[unit_type]
    validate_unit(unit_type, unit_in)
    validate_unit(unit_type, unit_out)
    conv_factor = unit_dict[unit_in] / unit_dict[unit_out]
    if verbose:
        print("{} to {} conversion factor: {}".format(unit_in, unit_out, conv_factor))
    return conv_factor

def convert_unit(quantity, unit_type, unit_in, unit_out, verbose=False):
    conv_factor = conversion_factor(unit_type, unit_in, unit_out)
    converted_quantity = quantity * conv_factor
    if verbose:
        print("{} {} = {} {}".format(quantity, unit_in, converted_quantity, unit_out))
    return converted_quantity

