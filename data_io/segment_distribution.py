import numpy as np
import os
import math
import itertools
from genesys.global_config import global_paths

#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
#* Detector and segment distribution
#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*

def get_local_band_detector_segment_dict(rank, size, band_detector_segment_dict):
    """
    This function determines the dictionary of detectors and segments for the particular MPI rank.
    The function takes as input the rank, the total number of MPI processes and the global dictionary of all detectors and segments.
    """
    # Counting the total number of segments
    num_segments = count_segments(band_detector_segment_dict)

    # Counting number of segments per rank and start and stop index
    num_segments_per_process = math.ceil(num_segments / size)
    start_index = rank * num_segments_per_process
    stop_index = (rank + 1) * num_segments_per_process

    # Right now, this is a slow hack job. Can be made much better.
    local_band_detector_segment_dict = {}
    count = -1
    seg_length = 0
    for band in sorted(list(band_detector_segment_dict.keys())):
        for detector in sorted(list(band_detector_segment_dict[band].keys())):
            for segment in band_detector_segment_dict[band][detector]:
                count += 1
                if count < start_index:
                    pass
                elif count >= start_index and count < stop_index:
                    if band not in list(local_band_detector_segment_dict.keys()):
                        local_band_detector_segment_dict[band] = {}
                    if detector not in list(local_band_detector_segment_dict[band].keys()):
                        local_band_detector_segment_dict[band][detector] = [] 
                    local_band_detector_segment_dict[band][detector].append(segment)
                    seg_length += 1
                else:
                    break
            if count >= stop_index:
                break

    return local_band_detector_segment_dict

#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*

def get_band_detector_time_dictionary_uniform(band_detector_dict, segment_length, segment_list):
    start_times = segment_list * segment_length
    stop_times = (segment_list + 1) * segment_length
    segment_time_list = list(zip(segment_list, start_times, stop_times))

    band_detector_time_dict = {}

    for band_name in list(band_detector_dict.keys()):
        band_detector_time_dict[band_name] = {}
        for detector_name in band_detector_dict[band_name]:
            band_detector_time_dict[band_name][detector_name] = segment_time_list

    return band_detector_time_dict

#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*

def count_segments(band_detector_segment_dict):
    # Counting the total number of segments
    num_segments = 0
    for band in list(band_detector_segment_dict.keys()):
        for detector in list(band_detector_segment_dict[band].keys()):
            num_segments += len(band_detector_segment_dict[band][detector])
    return num_segments

def get_detector_list(band_detector_segment_dict):
    detector_name_list = []
    for band in list(band_detector_segment_dict.keys()):
        for detector in list(band_detector_segment_dict[band].keys()):
            detector_name_list += [band+"--"+detector]
    return detector_name_list

# Needs to be modified for new definition of segment
def get_formatted_detector_segment_list(band_detector_segment_dict):
    formatted_det_seg_list = ""
    for band in list(band_detector_segment_dict.keys()):
        for detector in list(band_detector_segment_dict[band].keys()):
            formatted_det_seg_list += "{}--{}: {}\n".format(band, detector, ', '.join([str(seg) for seg in band_detector_segment_dict[band][detector]]))
    return formatted_det_seg_list

def get_total_simulation_time(band_detector_segment_dict):
    total_time = 0
    for band_name in list(band_detector_segment_dict.keys()):
        for detector_name in list(band_detector_segment_dict[band_name].keys()):
            total_time += sum([seg[2] - seg[1] for seg in band_detector_segment_dict[band_name][detector_name]])
    return total_time
#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
