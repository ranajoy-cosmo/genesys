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
    An empty class that can be used as a template for any object in the Genesys ecosystem.
    """
    def __init__(self):
        """
        Will be overridden by the child class if it has its own __init__
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
    This method allows for adding an external method to a class
    From: https://medium.com/@mgarod/dynamically-add-a-method-to-a-class-in-python-c49204b85bd6
    """
    def decorator(func):
        @wraps(func) 
        def wrapper(self, *args, **kwargs): 
            return func(*args, **kwargs)
        setattr(cls, func.__name__, wrapper)
        # NOTE we are not binding func, but wrapper which accepts self but does exactly the same as func
        return func # Returning func means func can still be used normally
    return decorator

# The subsequent methods are written with a @add_method(Genesys_Class) before its definition
# This allows them to be used both independently as themselves as well as as methods of Genesys_Class

#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
# Prompter utility
#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*

@add_method(Genesys_Class)
def prompt(text, nature='general', outstream=sys.stdout):
    """
    Flush out any text to the outstream provided immediately.
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
# Parameter loading method
#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*

@add_method(Genesys_Class)
def load_param_file(file_path=None):
    """
    Loads the parameters into a python dictionary
    """
    yaml = YAML(typ='safe')
    with open(file_path, 'r') as f:
        yaml_dict = yaml.load(f)
    return yaml_dict

@add_method(Genesys_Class)
def dump_param_file(yaml_dict, file_path=None):
    """
    Dumps the python dictionary into a yaml file
    """
    yaml = YAML(typ='safe')
    with open(file_path, 'w') as f:
        yaml.dump(yaml_dict, f)

#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
# This block defines the global paths
#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*

global_paths = {}

current_dir = os.path.dirname(__file__)
# The directory where large files are written down. for example the scratch system at nersc.
# Please do not set this somewhere with very limited storage as it'll fill up fast.
# Ideally it should have space of the order of a few terabytes.
# storage_dir = current_dir

# Forming the paths
global_paths['base_dir'] = current_dir
global_paths['storage_dir'] = storage_dir
global_paths['output_dir'] = os.path.join(global_paths['storage_dir'], 'output') 
global_paths['maps_dir'] = os.path.join(global_paths['storage_dir'], 'maps') 
global_paths['spectra_dir'] = os.path.join(global_paths['storage_dir'], 'spectra')
global_paths['data_dir'] = os.path.join(global_paths['storage_dir'], 'data')
#  global_paths['camb_params_dir'] = os.path.join(current_dir, 'spectra', 'camb_params')
#  global_paths['instruments_dir'] = os.path.join(current_dir, 'instruments')

# Adds the dictionary as class member of genesys_class
Genesys_Class.global_paths = global_paths
