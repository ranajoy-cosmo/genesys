import numpy as np
import healpy as hp
import sys
import os

def get_ts_signal(bolo_list, segment, config):
    if config.ind_bolo:
        t_stream = get_signal_individual_bolo(bolo_list[0], segment, config)
    if config.pair_difference:
        t_stream = get_signal_individual_bolo(bolo_list[0], bolo_list[1], segment, config)
    if config.pair_difference_subtract_template:
