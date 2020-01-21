import os
import sys
import copy
from termcolor import colored
from functools import wraps

#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
# THIS BLOCK IS THE ONLY PLACE A USER NEEDS TO MODIFY 
#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
# THIS CAN ALSO BE SET IN THE GLOBAL PATHS BLOCK. CHECK DESCRIPTION THERE
storage_dir = "/mn/stornext/d14/comap/ranajoyb/Genesys_Files"

#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
# THIS BLOCK DEFINES THE GENERIC CLASS TO BE INHERITED BY ALL OTHER CLASSES IN THE genesys ECOSYSTEM
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
        self.params = {}
        if other != None:
            self.copy_attributes(other)

    def copy_params(self, *param_list):
        """
        COPY PARAMS TO SELF
        IF param_list IS SPECIFIED, ONLY THOSE PARAMS ARE COPIED
        """
        for param in param_list:
            self.params.update(param)

    def copy_attributes(self, other):
        """
        COPY JUST THE INSTANCE VARIABLES INCLUDING params. METHODS AND CLASS VARIABLES ARE NOT COPIED
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
    else:
        raise Exception("Nature of prompt -- no")
    outstream.flush()

#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
# PARAMETER LOADING METHOD
#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*

@add_method(Genesys_Class)
def load_yaml_file(file_path=None):
    """
    LOADS THE PARAMETERS INTO A PYTHON DICTIONARY
    """
    yaml = YAML(typ='safe')
    with open(file_path, 'r') as f:
        yaml_dict = yaml.load(f)
    return yaml_dict

@add_method(Genesys_Class)
def dump_yaml_file(yaml_dict, file_path=None):
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
