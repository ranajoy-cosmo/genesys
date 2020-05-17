import numpy as np
import copy
from genesys import Genesys_Class
from genesys.numerical.unit_conversion import Unit_Converter

uc = Unit_Converter()

class Segment_Distributor(Genesys_Class):
    def __init__(self, data_block_list, num_segments_per_data_block, size, rank):
        self.data_block_list_local = self.get_local_data_block_list(data_block_list, size, rank)
        self.num_data_blocks_local = len(self.data_block_list_local)
        self.num_data_blocks_global = len(data_block_list)
        self.num_segments_local = self.num_data_blocks_local * num_segments_per_data_block
        self.num_segments_global = self.num_data_blocks_global * num_segments_per_data_block

    def get_local_data_block_list(self, data_block_list, size, rank):
        num_blocks = len(data_block_list)
        min_len = num_blocks // size
        num_proc_max = num_blocks - min_len*size 
        if rank < num_proc_max:
            start = rank*(min_len+1)
            stop = (rank+1)*(min_len+1)
        else:
            start = num_proc_max*(min_len+1) + (rank - num_proc_max)*min_len
            stop = num_proc_max*(min_len+1) + (rank + 1 - num_proc_max)*min_len
        data_block_list_local = data_block_list[start:stop]
        return data_block_list_local

    def get_segment_list(self, data_block, num_segments_per_data_block):
        """
        List of the segments in the data block
        """
        segment_list = np.arange(1, num_segments_per_data_block+1) + (data_block-1)*num_segments_per_data_block
        return segment_list

    def set_observation_times(self, segment_length):
        self.observation_time_global = self.num_segments_global * segment_length
        self.observation_time_local = self.num_segments_local * segment_length

    def count_detectors(self, channel_detector_dict):
        det_len_list = [len(dets) for dets in channel_detector_dict.values()]
        return sum(det_len_list)

def print_time_in_units(time_in, unit_in, units_out):
    print(f"{time_in} {unit_in} =")
    for unit_out in units_out:
        time_out = uc.convert_unit(time_in, 'time', unit_in, unit_out)
        print(f"\t{time_out} {unit_out}")
