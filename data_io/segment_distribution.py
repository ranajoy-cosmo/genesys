import numpy as np
import os
import math
from genesys.global_config import global_paths

#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
#* Detector and segment distribution
#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
def get_local_detector_segment_dict(rank, size, all_detector_segment_dict):
    """
    This function determines the dictionary of detectors and segments for the particular MPI rank.
    The function takes as input the rank, the total number of MPI processes and the global dictionary of all detectors and segments.
    """
    # Counting the total number of segments
    num_segments = 0
    all_detector_segment_length = {}
    for detector in list(all_detector_segment_dict.keys()):
        segment_length = len(all_detector_segment_dict[detector])
        all_detector_segment_length[detector] = segment_length
        num_segments += segment_length

    num_segments_per_process = math.ceil(num_segments / size)

    start_index = rank * num_segments_per_process
    stop_index = (rank + 1) * num_segments_per_process

    # Right now, this is a slow hack job. Can be made much better.
    local_detector_segment_dict = {}
    count = -1
    seg_length = 0
    for detector in sorted(list(all_detector_segment_dict.keys())):
        for segment in all_detector_segment_dict[detector]:
            count += 1
            if count < start_index:
                pass
            elif count >= start_index and count < stop_index:
                if detector not in list(local_detector_segment_dict.keys()):
                    local_detector_segment_dict[detector] = [] 
                local_detector_segment_dict[detector].append(segment)
                seg_length += 1
            else:
                break
        if count >= stop_index:
            break

    return local_detector_segment_dict, seg_length

def get_bolo_pair_segment_list(rank, size, segment_list):
    num_segments = len(segment_list)

    if num_segments % size != 0:
        num_segments_per_process = num_segments/size + 1
    else:
        num_segments_per_process = num_segments/size

    pair_segment_list = segment_list[rank*num_segments_per_process : (rank+1)*num_segments_per_process]

    return pair_segment_list
#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
