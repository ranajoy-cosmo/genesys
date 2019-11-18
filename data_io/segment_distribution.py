import numpy as np
import os
import math
import itertools
from genesys import global_paths

#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
#* Detector and segment distribution
#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*

def get_local_channel_detector_segment_dict(rank, size, channel_detector_dict, num_segments_per_det):
    """
    This function determines the dictionary of detectors and segments for the particular MPI rank.
    The function takes as input the rank, the total number of MPI processes and the global dictionary of all detectors and segments.
    """
    # Counting the total number of segments
    num_segments_total = count_detectors(channel_detector_dict) * num_segments_per_det

    # Counting number of segments per rank and start and stop index
    num_segments_per_process = math.ceil(num_segments_total / size)
    start_index = rank * num_segments_per_process
    stop_index = (rank + 1) * num_segments_per_process

    # Right now, this is a slow hack job. Can be made much better.
    local_channel_detector_segment_dict = {}
    count = -1
    for channel in sorted(list(channel_detector_dict.keys())):
        for detector in channel_detector_dict[channel]:
            for segment in range(num_segments_per_det):
                count += 1
                if count < start_index:
                    pass
                elif count >= start_index and count < stop_index:
                    if channel not in list(local_channel_detector_segment_dict.keys()):
                        local_channel_detector_segment_dict[channel] = {}
                    if detector not in list(local_channel_detector_segment_dict[channel].keys()):
                        local_channel_detector_segment_dict[channel][detector] = [] 
                    local_channel_detector_segment_dict[channel][detector].append(segment)
                else:
                    break
            if count >= stop_index:
                break

    return local_channel_detector_segment_dict

#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*

def get_channel_list(channel_detector_segment_dict):
    # GET A LIST OF THE CHANNELS THAT ARE PRESENT
    channel_list = list(channel_detector_segment_dict.keys())
    return channel_list

def get_formatted_detector_list(channel_detector_segment_dict):
    # GET A LIST OF THE CHANNELS THAT ARE PRESENT
    formatted_detector_list = {}.fromkeys(get_channel_list(channel_detector_segment_dict))
    for channel in list(formatted_detector_list.keys()):
        formatted_detector_list[channel] = list(channel_detector_segment_dict[channel].keys())
    return formatted_detector_list

def count_detectors(channel_detector_segment_dict):
    # COUNT THE TOTAL NUMBER OF DETECTORS IN THE SIMULATION
    num_detectors = 0
    for channel in list(channel_detector_segment_dict.keys()):
        num_detectors += len(channel_detector_segment_dict[channel])
    return num_detectors

def count_segments(channel_detector_segment_dict):
    # COUNT THE TOTAL NUMBER OF SEGMENTS IN THE SIMULATION
    num_segments = 0
    for channel in list(channel_detector_segment_dict.keys()):
        for detector in list(channel_detector_segment_dict[channel].keys()):
            num_segments += len(channel_detector_segment_dict[channel])
    return num_segments
#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
