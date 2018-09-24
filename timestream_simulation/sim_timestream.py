#!/usr/bin/env python

import numpy as np
import sys
import os
import importlib
import time
import pickle
from termcolor import colored
import genesys.data_io.segment_distribution as segd
import genesys.data_io.data_io as dio
from genesys.timestream_simulation.detector import Detector
from genesys.global_config import global_paths
from genesys.utilities import prompt

"""
This is the master executable for the timestream generation pipeline.
Data distribution is performed from here.
This code is embarassingly parallel. that is, there is no explicit communication between any processes.
There are two options for running this
1) Serial : A single process will iterate over all detectors and data segments.
2) MPI : The different detectors and data segments will be distributed to different processes.
How to execute this file.
python sim_timestream.py <CONFIG_FILE_PATH> <RUN_TYPE>
where, CONFIG_FILE_PATH is the absolute path to the config file to be passes to the code and RUN_TYPE is either 'run_serial' or 'run_mpi'
"""

#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
#* The master section which distributes the data and collects it
#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*

def run_simulation():
    count = 0
    cumulative_time_taken = 0

    for detector_name in local_detector_name_list:
        detector = Detector(detector_name, config)
        for segment in local_detector_segment_dict[detector_name]:
            segment_start_time = time.time()
            count += 1
            prompt("Rank {} doing detector {} and segment {}\n".format(rank, detector_name, segment+1))
            #  if config.sim_type == "signal":
                #  bolo.simulate_timestream_signal(segment)
            #  else:
                #  bolo.simulate_timestream_template(segment)
            segment_stop_time = time.time()
            time_taken = segment_stop_time - segment_start_time
            cumulative_time_taken += time_taken
            prompt("Rank {} finished {} of {} segments.\nTime taken: {}s. Total time taken: {}, Projected time: {}\n".format(rank, count, local_num_segments, time_taken, cumulative_time_taken, cumulative_time_taken*local_num_segments/count))

def make_data_dirs():
    """
    This function will create all the necessary directories under which the simulated data will be written to.
    This function is, under normal circumstances, only called by rank 0, and a barrier is set so that all necessary directories are generated before simulation begins.
    Known issues: Rarely, parallel and independent runs of the code might conflict leading to a directory not being written. Report if this occurs.
    """
    # The main sim directory under which all simulated timestream data and corresponding maps will be written to.
    sim_dir = dio.get_path_to_sim_dir(config.sim_tag)
    print(sim_dir)
    if not os.path.exists(sim_dir):
        os.makedirs(sim_dir, exist_ok=True)

    # The scan directory under which all simulated timestream data will be written to.
    scan_dir = dio.get_path_to_scan_dir(config.sim_tag, config.scan_tag)
    if not os.path.exists(scan_dir):
        os.makedirs(scan_dir, exist_ok=True)
#
    #  # Iterate through the local list of detectors and segments and generate the directories when appropriate.
    #  # If sim_mode is set to 'overwrite', existing data will be deleted.
    #  for detector_name in list(config.detector_segment_dict.keys()):
        #  det_dir = dio.get_path_to_det_dir(config.sim_tag, config.scan_tag, detector_name)
        #  if not os.path.exists(det_dir):
            #  os.makedirs(det_dir, exist_ok=True)
        #  for segment in config.detector_segment_dict[detector_name]:
            #  segment_dir = dio.get_path_to_segment_dir(config.sim_tag, config.scan_tag, detector_name, segment)
            #  if not os.path.exists(segment_dir):
                #  os.makedirs(segment_dir, exist_ok=True)

def start_message():
    display_string = colored("\n#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*\n", color='blue')
    display_string += colored("#* ", color='blue') + colored("GENESYS\n", color='green')
    display_string += colored("#* ", color='blue') + colored("BEGINNING TIMESTREAM SIMULATION\n", color='green')
    display_string += colored("#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*\n", color='blue')
    display_string += "RUN TYPE: {}\n".format(run_type)
    display_string += "# OF PROCESSES: {}\n".format(size)
    display_string += "SIMULATION TAG: {}\n".format(config.sim_tag)
    display_string += "SCAN TAG: {}\n".format(config.scan_tag)
    display_string += "SIMULATION POLARISATION TYPE: {}\n".format(config.sim_pol_type)
    display_string += "COORDINATE SYSTEM: {}\n".format(config.coordinate_system)
    display_string += "SIMULATION TOD TYPE: {}\n".format(config.tod_type)
    if config.tod_type == "gradient":
        display_string += "GRADIENT TYPES: {}\n".format(config.gradient_type)
    detector_list = sorted(list(config.detector_segment_dict.keys()))
    display_string += "DETECTOR LIST: {}\n".format(', '.join(detector_list))
    display_string += "# OF DETECTORS: {}\n".format(len(detector_list))
    formatted_det_seg_list = ""
    num_seg = 0
    for detector in detector_list:
        num_seg += len(config.detector_segment_dict[detector])
        formatted_det_seg_list += "{}: {}\n".format(detector, ', '.join([str(seg) for seg in config.detector_segment_dict[detector]]))
    display_string += "DETECTOR-SEGMENT LIST:\n{}".format(formatted_det_seg_list)
    display_string += "TOTAL NUMBER OF SEGMENTS : {}\n".format(num_seg)
    display_string += "SCAN STRATEGY NAME: {}\n".format(config.scan_strategy_name)
    display_string += "PRECESSION TIME: {} seconds\n".format(config.t_prec)
    display_string += "SPIN TIME: {} seconds\n".format(config.t_spin)
    display_string += "SAMPLING RATE: {} Hz\n".format(config.sampling_rate)
    display_string += "ALPHA: {} degrees\n".format(config.alpha)
    display_string += "BETA: {} degrees\n".format(config.beta)
    display_string += "WRITE FIELD: {}\n".format(config.ts_data_products)
    display_string += "NOTES: {}\n".format(config.notes)
    display_string += colored("#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*\n\n", color='blue')
    prompt(display_string)

#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
#* Main function definition. This is where the code begins when executed
#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
if __name__=="__main__":
    start_time_main = time.time()

    # run_type can be either "mpi" or "serial"
    run_type = sys.argv[1]
    assert run_type in ["mpi", "serial"], "Please enter \"mpi\" OR \"serial\" as the second argument."
    # config_file is provided as a dot separated path to the config file with genesys as the parent directory
    config_file = sys.argv[2]

    # This block just sets up the rank amd size of the simulation run
    if run_type == "mpi":
        from mpi4py import MPI          # This module is imported here to be compatible with systems not having mpi support
        comm = MPI.COMM_WORLD
        size = comm.Get_size()          # Number of MPI processes
        rank = comm.Get_rank()          # Rank of this particular MPI process. Ranges from 0 to (size - 1)
    else:
        size = 1
        rank = 0
    
    # Importing the general simulation parameters
    config = importlib.import_module(config_file).config
    if config.scan_strategy_module != None:
        scan_strategy = importlib.import_module("genesys.scan_strategy." + config.scan_strategy_module).scan_strategy
        config.__dict__.update(scan_strategy.__dict__)

    # Getting the local list of detectors and segments
    local_detector_segment_dict, local_num_segments = segd.get_local_detector_segment_dict(rank, size, config.detector_segment_dict)
    local_detector_name_list = sorted(list(local_detector_segment_dict.keys()))
    detector_name_list = sorted(list(config.detector_segment_dict.keys()))

    # Performing checks
    check_flag = 1
    for detector_name in local_detector_name_list:
        det = Detector(detector_name, config)
        check_flag *= det.check_flag

    if run_type == "mpi":
        check_flag = comm.allreduce(check_flag, op=MPI.PROD)
    else:
        check_flag = check_flag_local

    assert check_flag==1, "Did not run simulation as all checks were not met. Terminating run. Rank {} of {}\n".format(rank, size)

    if rank == 0:
        start_message()
        prompt(colored("\n#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*\n", color='blue'))
        for detector_name in detector_name_list:
            det = Detector(detector_name, config)
            det.display_params()
        prompt(colored("#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*\n\n", color='blue'))
        make_data_dirs()

    if run_type == "mpi":
        comm.Barrier()

    formatted_det_seg_list = ""
    for detector_name in local_detector_name_list:
        formatted_det_seg_list += "{}: {}\n".format(detector_name, ', '.join([str(seg) for seg in config.detector_segment_dict[detector_name]]))
    display_string = colored("#* RANK: {} #*#*#*#*#*#*#*#*#*#*\n".format(rank), color="green")
    display_string += "DOING DETECTORS: {}\n".format(', '.join(local_detector_name_list))
    display_string += "LOCAL DETECTOR-SEGMENT LIST : {}".format(formatted_det_seg_list)
    display_string += "TOTAL NUMBER OF SEGMENTS : {}\n".format(local_num_segments)
    display_string += colored("#*#*#*#*#*#*#*#*#*#*\n\n", color="green")
    prompt(display_string)

    #  Running the main simulation routine. This is common to both mpi and serial runs.
    run_simulation()

    stop_time_main = time.time()
    prompt("***RANK {} COMPLETE***. Total time taken : {}s\n".format(rank, stop_time_main - start_time_main))
#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
