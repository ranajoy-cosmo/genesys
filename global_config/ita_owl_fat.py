import os
from genesys.utilities import Generic_Class

"""
In this module I define the application wide settings for a particular machine.
This particular case is that of the OWL 25-28 cluster ar ITA, Oslo.
The settings are stored as instances of the two objects
1) global_system : Contains the settings and information about the host system
2) global_paths : Contains the paths and links to important directories on the host system
"""

#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
# Global system settings
#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*

global_system = Generic_Class()
# Keys
# global_system.host_name : The hostname of the system the app is running on
# global_system.num_cores_per_socket : The number of cores(processors) per socket of the machine
# global_system.num_sockets_per_node : The number of sockets per node of the system
# global_system.num_threads_per_core : The maximum number of threads per core
# global_system.memory_per_node : The memory in GigaBytes per node of the machine

global_system.host_name = "ita_owl_fat"
global_system.num_processors_per_socket = 18
global_system.num_sockets_per_node = 4
global_system.num_threads_per_core = 2
global_system.memory_per_node = 1536.0                              #GB

# The following parameters are calculated from the values entered above
# Keys
# global_system.num_cores_per_node : The total number of cores per node
# global_system.num_threads_per_node : The maximum number of threads possible per node

global_system.num_cores_per_node = global_system.num_processors_per_socket * global_system.num_sockets_per_node
global_system.num_threads_per_node = global_system.num_cores_per_node * global_system.num_threads_per_core

#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
# Global path settings
#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*

global_paths = Generic_Class()
# Keys
# global_paths.home_dir : The home directory of the user
# global_paths.data_dir : The directory where large files are written down. For example the scratch system at NERSC.
# global_paths.base_dir : The directory of the simulation code
# The following directories are the locations on the user's system where large files get stored. Softlinks to these will be generated for easy access of the content of these directories.
# global_paths.output_dir : The directory for I/O during execution and storing simulation data
# global_paths.maps_dir : The directory under the data directory where the input maps are stored
# global_path.spectra_dir : The directory under the data directory where the CMB spectra are stored 

global_paths.home_dir = "/mn/stornext/u3/ranajoyb"
global_paths.data_dir = "/mn/stornext/d14/comap/ranajoyb"

global_paths.base_dir = os.path.join(global_paths.home_dir, 'genesys')
# The following paths will be set on the user's machine in the 'data' directory. Soft links to these will be created in the 'genesys' working directory.
global_paths.output_dir = os.path.join(global_paths.data_dir, 'genesys_output')
global_paths.maps_dir = os.path.join(global_paths.data_dir, 'genesys_maps')
global_paths.spectra_dir = os.path.join(global_paths.data_dir, 'genesys_spectra')
