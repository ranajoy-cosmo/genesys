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
    def __init__(self, other=None):
        """
        EXPECT THIS TO BE OVERRIDDEN BY THE CHILD CLASS
        """
        if other != None:
            self.copy_attributes(other)

    def copy_params(self, *params_list):
        """
        DEEPCOPY PARAMS TO SELF
        """
        if not hasattr(self, 'params'):
            self.params = {}
        for params in params_list:
            self.params.update(params)

    def copy_attributes(self, other):
        """
        COPY JUST THE INSTANCE VARIABLES. METHODS AND CLASS VARIABLES ARE NOT COPIED
        """
        self.__dict__.update(other.__dict__) 

    def return_copy(self):
        """
        RETURN AN ENTIRE COPY OF SELF
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
        outstream.write(text)
    elif nature == 'info':
        outstream.write(colored('INFO: ', 'yellow') + text)
    elif nature == 'warning':
        outstream.write(colored('WARNING: ', 'red') + text)
    outstream.flush()

#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
# PARAMETER LOADING METHOD
#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*

@add_method(Genesys_Class)
def load_param_file(file_path=None, file_handle=None):
    """
    LOADS THE PARAMETERS INTO A PYTHON DICTIONARY
    PRECEDENCE:
        file_path > file_handle
    """
    yaml = YAML(typ='safe')
    if file_path != None:
        with open(file_path, 'r') as file_handle:
            param_data = yaml.load(file_handle)
    else:
        param_data = yaml.load(file_handle)

    return param_data

@add_method(Genesys_Class)
def dump_param_file(yaml_dict, file_path=None, file_handle=None):
    """
    DUMPS THE PYTHON DICTIONARY INTO A yaml FILE
    PRECEDENCE:
        file_path > file_handle
    """
    yaml = YAML(typ='safe')
    if file_path != None:
        with open(file_path, 'w') as file_handle:
            yaml.dump(yaml_dump, file_handle)
    else:
        yaml.dump(yaml_dump, file_handle)


#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
# THIS BLOCK DEFINES THE GLOBAL SYSTEM CONFIGURATION AND THE GLOBAL PATHS
#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*

current_dir = os.path.dirname(__file__)

with open(os.path.join(current_dir, 'global_config', global_config_file), 'r') as f:
    global_param_dict = load_param_file(file_handle=f)

global_system = global_param_dict['global_system']
global_paths = global_param_dict['global_paths']

# FORMING THE PATHS
global_paths['base_dir'] = current_dir
global_paths['output_dir'] = os.path.join(global_paths['storage_dir'], global_paths['real_output_dir']) 
global_paths['maps_dir'] = os.path.join(global_paths['storage_dir'], global_paths['real_maps_dir']) 
global_paths['spectra_dir'] = os.path.join(global_paths['storage_dir'], global_paths['real_spectra_dir']) 
global_paths['data_dir'] = os.path.join(global_paths['storage_dir'], global_paths['real_data_dir']) 
global_paths['camb_params_dir'] = os.path.join(current_dir, 'spectra', 'camb_params') 
global_paths['instruments_dir'] = os.path.join(current_dir, 'instruments')

# ADDS THESE DICTIONARIES AS CLASS MEMBERS OF GENESYS_CLASS
Genesys_Class.global_system = global_system
Genesys_Class.global_paths = global_paths
