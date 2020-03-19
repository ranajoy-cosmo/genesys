import numpy as np
import os
import math
import itertools
from genesys import global_paths
from genesys import Genesys_Class

class Segment_Distributor(Genesys_Class):
    """
    This class creates a dictionary of the channels, detectors and segments such that they are most evenly distributed over the available MPI processes.
    The structure of the dictionary has the following hierarchy: segment - channel - detector
    {
    """
    def __init__(self, rank, size, channel_detector_dict, num_segments_per_det):
        self.channel_detector_segment_dict = self.get_channel_detector_segment_dict(rank, size, channel_detector_dict, num_segments_per_det)

    def get_channel_detector_segment_dict(self, rank, size, channel_detector_dict, num_segments_per_det):
        """
        This function determines the dictionary of detectors and segments for the particular MPI rank.
        The function takes as input the rank, the total number of MPI processes and the global dictionary of all detectors and segments.
        """
        # COUNTING THE TOTAL NUMBER OF SEGMENTS.
        num_detectors = 0
        for channel_name in list(channel_detector_dict.keys()):
            num_detectors += len(channel_detector_dict[channel_name])
        num_segments_total = num_detectors * num_segments_per_det

        # COUNTING NUMBER OF SEGMENTS PER RANK AND START AND STOP INDEX
        num_segments_per_process = math.ceil(num_segments_total / size)
        start_index = rank * num_segments_per_process
        stop_index = (rank + 1) * num_segments_per_process

        # RIGHT NOW, THIS IS A SLOW HACK JOB. CAN BE MADE MUCH BETTER.
        channel_detector_segment_dict = {}
        count = -1
        for channel_name in sorted(list(channel_detector_dict.keys())):
            for detector_name in channel_detector_dict[channel_name]:
                for segment in range(num_segments_per_det):
                    count += 1
                    if count < start_index:
                        pass
                    elif count >= start_index and count < stop_index:
                        if channel_name not in list(channel_detector_segment_dict.keys()):
                            channel_detector_segment_dict[channel_name] = {}
                        if detector_name not in list(channel_detector_segment_dict[channel_name].keys()):
                            channel_detector_segment_dict[channel_name][detector_name] = [] 
                        channel_detector_segment_dict[channel_name][detector_name].append(segment+1)
                    else:
                        break
                if count >= stop_index:
                    break

        return channel_detector_segment_dict

    #*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*

    def channel_list(self):
        # GET A LIST OF THE CHANNELS THAT ARE PRESENT
        return list(self.channel_detector_segment_dict.keys())

    def num_channels(self):
        return len(self.channel_list())

    def detectors_in_channel(self, channel_name):
        # GET A LIST OF THE DETECTOR IN A PARTICULAR CHANNELS
        return list(self.channel_detector_segment_dict[channel_name].keys())

    def num_detectors_in_channel(self, channel_name):
        return len(self.detectors_in_channel(channel_name))

    def num_detectors(self):
        # COUNT THE TOTAL NUMBER OF DETECTORS IN THE SIMULATION
        num_detectors = 0
        for channel_name in self.channel_list():
            num_detectors += len(self.detectors_in_channel(channel_name))
        return num_detectors

    def formatted_detector_list(self):
        # GET A FORMATTED LIST OF THE CHANNELS AND DETECTORS
        fmt_detector_list = {}.fromkeys(self.channel_list())
        for channel_name in list(self.channel_list()):
            fmt_detector_list[channel_name] = self.detectors_in_channel(channel_name)
        return fmt_detector_list

    def num_segments(self):
        # COUNT THE TOTAL NUMBER OF SEGMENTS IN THE SIMULATION
        num_segments = 0
        for channel_name in self.channel_list():
            for detector_name in self.detectors_in_channel(channel_name):
                num_segments += len(self.channel_detector_segment_dict[channel_name][detector_name])
        return num_segments
    #*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
