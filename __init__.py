import os
import sys
import copy
from termcolor import colored
from functools import wraps
from ruamel.yaml import YAML

# The its_owl_fat file is used by Ranajoy at ITA. Other users, please set the lines above to your preferred configuration file.
global_config_file = "ita_owl_fat.yaml"  # Only this line needs to change between systems

#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
# This block defines the generic class to be inherited by all other classes in the genesys ecosystem
#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*

class Genesys_Class:
    """
    AN EMPTY CLASS THAT CAN BE USED AS A TEMPLATE FOR ANY OBJECT.
    BASIC FEATURES:
        copy_other: COPY THE ATTRIBUTES OF ANOTHER OBJECT
    """

    def __init__(self):
        """
        Does nothing. Expect this to be overridden by the child class
        """
        pass

    def copy_attributes(self, other):
        """
        Copy just the instance attributes. Methods and class attributes are not copied
        """
        self.__dict__ = copy.deepcopy(other.__dict__) 

    def return_copy(self):
        """
        Return an entire copy of self
        """
        return copy.deepcopy(self)

#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*

def add_method(cls):
    """
    This method allows for adding an external method to a class
    From: https://medium.com/@mgarod/dynamically-add-a-method-to-a-class-in-python-c49204b85bd6
    """
    def decorator(func):
        @wraps(func) 
        def wrapper(self, *args, **kwargs): 
            return func(*args, **kwargs)
        setattr(cls, func.__name__, wrapper)
        # Note we are not binding func, but wrapper which accepts self but does exactly the same as func
        return func # returning func means func can still be used normally
    return decorator

# The subsequent methods are written with a @add_method(Genesys_Class) before its definition
# This allows them to be used both independently as themselves as well as as methods of Genesys_Class

#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
# PROMPTER UTILITY
#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
#
#  @add_method(Genesys_Class)
#  def prompt(text, outstream=sys.stdout, color=None):
    #  """
    #  Flush out any text to the outstream provided immediately.
    #  """
    #  if color:
        #  text = colored(text, color)
    #  outstream.write(text)
    #  outstream.flush()

#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
# PARAMETER LOADING METHOD
#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*

@add_method(Genesys_Class)
def load_param_file(file_handle=None, file_path=None):
    """
    Loads the parameters into a python dictionary from the file_handle
    Precedence:
        file_handle > file_path
    """
    yaml = YAML(typ='safe')
    if file_handle != None:
        param_data = yaml.load(file_handle)
    else:
        with open(file_path, 'r') as file_handle:
            param_data = yaml.load(file_handle)

    return param_data

#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
# THIS BLOCK DEFINES THE GLOBAL SYSTEM CONFIGURATION AND THE GLOBAL PATHS
#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*

current_dir = os.path.dirname(__file__)

with open(os.path.join(current_dir, 'global_config', global_config_file), 'r') as f:
    global_param_dict = load_param_file(file_handle=f)

global_system = global_param_dict['global_system']
global_paths = global_param_dict['global_paths']

# Forming the paths
global_paths['base_dir'] = current_dir
global_paths['output_dir'] = os.path.join(global_paths['data_dir'], global_paths['rel_output_dir']) 
global_paths['maps_dir'] = os.path.join(global_paths['data_dir'], global_paths['rel_maps_dir']) 
global_paths['spectra_dir'] = os.path.join(global_paths['data_dir'], global_paths['rel_spectra_dir']) 
global_paths['camb_params_dir'] = os.path.join(current_dir, 'spectra', 'camb_params') 
global_paths['instruments_dir'] = os.path.join(current_dir, 'instruments')

# Adds these dictionaries as class members of Genesys_Class
Genesys_Class.global_system = global_system
Genesys_Class.global_paths = global_paths

#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
# Path naming conventions to be followed throught the genesys ecosystem
#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
