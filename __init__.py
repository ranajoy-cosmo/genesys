import os
import sys
import copy
from termcolor import colored
from functools import wraps
from ruamel.yaml import YAML

#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
# This block is the only place a user needs to modify 
#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
# This can also be set in the global paths block. Check description there
storage_dir = "/mn/stornext/d14/comap/ranajoyb/Genesys_Files"

#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
# This block defines the generic class to be inherited by all other classes in the genesys ecosystem
#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*

class Genesys_Class:
    """
    An empty class that can be used as a template for any object.
    Basic features:
        copy_other: Copy the attributes of another object
    """
    def __init__(self):
        """
        Expect this to be overridden by the child class
        """
        self.params = {}

    def copy_params(self, params, param_list=None):
        """
        Copy params to self
        If param_list is specified, only those params are copied
        """
        if not hasattr(self, 'params'):
            self.params = {}
        if param_list == None:
            self.params.update(params)
        else:
            self.params.update({k: params[k] for k in param_list})

    def copy_attributes(self, other):
        """
        Copy just the instance variables including params. Class methods and class variables are not copied
        """
        self.__dict__.update(other.__dict__) 

    def return_copy(self):
        """
        Return an entire copy of self
        """
        return copy.deepcopy(self)

#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*

def add_method(cls):
    """
    THIS METHOD ALLOWS FOR ADDING AN EXTERNAL METHOD TO A CLASS
    FROM: https://medium.com/@mgarod/dynamically-add-a-method-to-a-class-in-python-c49204b85bd6
    """
    def decorator(func):
        @wraps(func) 
        def wrapper(self, *args, **kwargs): 
            return func(*args, **kwargs)
        setattr(cls, func.__name__, wrapper)
        # NOTE WE ARE NOT BINDING FUNC, BUT WRAPPER WHICH ACCEPTS self BUT DOES EXACTLY THE SAME AS func
        return func # RETURNING func MEANS func CAN STILL BE USED NORMALLY
    return decorator

# THE SUBSEQUENT METHODS ARE WRITTEN WITH A @add_method(Genesys_Class) BEFORE ITS DEFINITION
# THIS ALLOWS THEM TO BE USED BOTH INDEPENDENTLY AS THEMSELVES AS WELL AS AS METHODS OF Genesys_Class

#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
# PROMPTER UTILITY
#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*

@add_method(Genesys_Class)
def prompt(text, nature='general', outstream=sys.stdout):
    """
    FLUSH OUT ANY TEXT TO THE outstream PROVIDED IMMEDIATELY.
    """
    if nature == 'general':
        outstream.write(text+'\n')
    elif nature == 'info':
        outstream.write(colored('\u25CF ', 'green') + text+'\n')
    elif nature == 'warning':
        outstream.write(colored('\u25CF ', 'red') + text+'\n')
    else:
        raise Exception(f"Nature of prompt '{nature}' not recognised")
    outstream.flush()

#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
# PARAMETER LOADING METHOD
#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*

@add_method(Genesys_Class)
def load_param_file(file_path=None):
    """
    LOADS THE PARAMETERS INTO A PYTHON DICTIONARY
    """
    yaml = YAML(typ='safe')
    with open(file_path, 'r') as f:
        yaml_dict = yaml.load(f)
    return yaml_dict

@add_method(Genesys_Class)
def dump_param_file(yaml_dict, file_path=None):
    """
    DUMPS THE PYTHON DICTIONARY INTO A yaml FILE
    """
    yaml = YAML(typ='safe')
    with open(file_path, 'w') as f:
        yaml.dump(yaml_dict, f)


#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
# THIS BLOCK DEFINES THE GLOBAL PATHS
#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*

global_paths = {}

current_dir = os.path.dirname(__file__)
# THE DIRECTORY WHERE LARGE FILES ARE WRITTEN DOWN. FOR EXAMPLE THE scratch SYSTEM AT NERSC.
# PLEASE DO NOT SET THIS SOMEWHERE WITH VERY LIMITED STORAGE AS IT'LL FILL UP FAST.
# IDEALLY IT SHOULD HAVE SPACE OF THE ORDER OF A FEW TERABYTES.
#  storage_dir = current_dir

# FORMING THE PATHS
global_paths['base_dir'] = current_dir
global_paths['storage_dir'] = storage_dir
global_paths['output_dir'] = os.path.join(global_paths['storage_dir'], 'output') 
global_paths['maps_dir'] = os.path.join(global_paths['storage_dir'], 'maps') 
global_paths['spectra_dir'] = os.path.join(global_paths['storage_dir'], 'spectra')
global_paths['data_dir'] = os.path.join(global_paths['storage_dir'], 'data')
#  global_paths['camb_params_dir'] = os.path.join(current_dir, 'spectra', 'camb_params')
#  global_paths['instruments_dir'] = os.path.join(current_dir, 'instruments')

# ADDS THE DICTIONARY AS CLASS MEMBER OF GENESYS_CLASS
Genesys_Class.global_paths = global_paths
