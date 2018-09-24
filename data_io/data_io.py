import numpy as np
import os
import math
from genesys.global_config import global_paths

#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
#* Path naming routines
#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
def get_path_to_sim_dir(sim_tag):
    return os.path.join(global_paths.output_dir, sim_tag)

def get_path_to_scan_dir(sim_tag, scan_tag):
    return os.path.join(global_paths.output_dir, sim_tag, scan_tag)

def get_path_to_det_dir(sim_tag, scan_tag, det_name):
    return os.path.join(global_paths.output_dir, sim_tag, scan_tag, det_name)

def get_path_to_segment_dir(sim_tag, scan_tag, det_name, segment):
    """
    The segment name is given as an integer. It is converted by this routine to the appropriate 4 digit string preceded by 0s.
    By convention, 1 is added to the segment name. So, the 0th segment will be named '0001'.
    """
    segment_name = str(segment+1).zfill(4)
    return os.path.join(global_paths.output_dir, sim_tag, scan_tag, det_name, segment_name)
#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*


#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
#* File existence
#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
def check_if_files_exist(segment_dir, file_list):
    exists = True
    for file_name in file_list:
        if not os.path.exists(os.path.join(segment_dir, file_name)):
            exists = False
            break

    return exists
#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
