import numpy as np
import copy
from genesys import Genesys_Class
from genesys.numerical.unit_conversion import Unit_Converter

uc = Unit_Converter()

class Segment_Distributor(Genesys_Class):
    def __init__(self, segment_list, size, rank, segment_length):
        self.segment_list_global = copy.copy(segment_list)
        num_total_segment = len(segment_list)
        if num_total_segment%size and rank == 0:
            self.prompt("Segments are not evenly distributed on processes.", nature='warning') 
        self.get_local_segment_list(segment_list, size, rank)
        self.observation_time_global = num_total_segment * segment_length
        self.observation_time_local = len(self.segment_list_local) * segment_length

    def get_local_segment_list(self, segment_list, size, rank):
        num_total_segment = len(segment_list)
        if num_total_segment % size:
            num_segments_per_process = num_total_segment // size + 1
        else:
            num_segments_per_process = num_total_segment // size
        self.segment_list_local = segment_list[rank*num_segments_per_process:(rank+1)*num_segments_per_process]

    def count_detectors(self, channel_detector_dict):
        det_len_list = [len(dets) for dets in channel_detector_dict.values()]
        return sum(det_len_list)

def print_time_in_units(time_in, unit_in, units_out):
    print(f"{time_in} {unit_in} =")
    for unit_out in units_out:
        time_out = uc.convert_unit(time_in, 'time', unit_in, unit_out)
        print(f"\t{time_out} {unit_out}")
