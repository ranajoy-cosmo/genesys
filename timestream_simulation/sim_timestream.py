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
#!/usr/bin/env python

import sys
import os
import importlib
import time
from termcolor import colored
from ..data_io import segment_distribution as segd
from ..data_io import data_io as dio
from .detector import Detector
from ..utilities import unit_conversion as uc
from ..utilities import prompt

#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
#* The master section which distributes the data and collects it
#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*

def run_simulation():
    count = 0
    cumulative_time_taken = 0

    for band_name in list(local_band_detector_segment_dict.keys()):
        for detector_name in list(local_band_detector_segment_dict[band_name].keys()):
            detector = Detector(band_name, detector_name, config)
            for segment in local_band_detector_segment_dict[band_name][detector_name]:
                segment_start_time = time.time()
                count += 1
                detector.segment_start_prompt(rank, segment)
                detector.simulate_timestream_data(segment)
                segment_stop_time = time.time()
                time_taken = segment_stop_time - segment_start_time
                cumulative_time_taken += time_taken
                prompt("Rank {} finished {} of {} segments.\nTime taken: {}s. Total time taken: {}, Projected time: {}\n".format(rank, count, local_num_segments, time_taken, cumulative_time_taken, cumulative_time_taken*local_num_segments/count))

#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
#* The start messages
#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*

def global_start_message():
    display_string = colored("\n#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*\n", color='blue')
    display_string += colored("#* ", color='blue') + colored("GENESYS\n", color='green')
    display_string += colored("#* ", color='blue') + colored("BEGINNING TIMESTREAM SIMULATION\n", color='green')
    display_string += colored("#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*\n", color='blue')
    display_string += "RUN TYPE: {}\n".format(run_type)
    display_string += "# OF PROCESSES: {}\n".format(size)
    display_string += "SIMULATION TAG: {}\n".format(config.sim_tag)
    display_string += "SCAN TAG: {}\n".format(config.scan_tag)
    display_string += "SPECIAL TAG: {}\n".format(config.special_tag)
    display_string += "SIMULATION POLARISATION TYPE: {}\n".format(config.sim_pol_type)
    display_string += "COORDINATE SYSTEM: {}\n".format(config.coordinate_system)
    display_string += "# OF DETECTORS: {}\n".format(len(detector_list))
    display_string += "DETECTOR LIST: {}\n".format(', '.join(detector_list))
    display_string += "TOTAL NUMBER OF SEGMENTS : {}\n".format(num_segments)
    display_string += "TOTAL SIMULATION TIME : {}s {}d {}y\n".format(total_simulation_time, uc.convert_unit('time', total_simulation_time, 'second', 'day'), uc.convert_unit('time', total_simulation_time, 'second', 'siderial year'))
    display_string += "SCAN STRATEGY NAME: {}\n".format(config.scan_strategy_name)
    display_string += "POLARISATION MODULATION TYPE: {}\n".format(config.polarisation_modulation)
    if config.polarisation_modulation == "continuous_HWP":
        display_string += "HWP initial phase: {} degrees\n".format(config.HWP_phase_ini)
        display_string += "HWP angular speed: {} RPM\n".format(config.HWP_rpm)
    if config.polarisation_modulation == "stepped_HWP":
        display_string += "HWP initial phase: {} degrees\n".format(config.HWP_phase_ini)
        display_string += "HWP step size: {} degrees\n".format(config.HWP_step)
        display_string += "HWP step duration: {} s\n".format(config.HWP_step_duration)
    display_string += "PRECESSION PERIOD: {} seconds\n".format(config.t_prec)
    display_string += "SPIN PERIOD: {} seconds\n".format(config.t_spin)
    display_string += "ALPHA: {} degrees\n".format(config.alpha)
    display_string += "BETA: {} degrees\n".format(config.beta)
    display_string += "TOD WRITE PRODUCTS: {}\n".format(config.ts_data_products)
    display_string += "NOTES: {}\n".format(config.notes)
    display_string += colored("#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*\n\n", color='blue')
    prompt(display_string)

def local_start_message():
    display_string = colored("#* RANK: {} #*#*#*#*#*#*#*#*#*#*\n".format(rank), color="green")
    display_string += "# OF DETECTORS: {}\n".format(len(local_detector_list))
    display_string += "DETECTOR LIST: {}\n".format(', '.join(local_detector_list))
    display_string += "LOCAL NUMBER OF SEGMENTS: {}\n".format(local_num_segments)
    display_string += "TOTAL LOCAL SIMULATION TIME : {}s {}d {}y\n".format(local_simulation_time, uc.convert_unit('time', local_simulation_time, 'second', 'day'), uc.convert_unit('time', local_simulation_time, 'second', 'siderial year'))
    display_string += colored("#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*\n\n", color="green")
    prompt(display_string)

#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
#* Main function definition. This is where the code begins when executed
#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*

if __name__=="__main__":
    start_time_main = time.time()

    # run_type can be either "mpi" or "serial"
    run_type = sys.argv[1]
    assert run_type in ["mpi", "serial"], "Please enter \"mpi\" OR \"serial\" as the second argument."
    # This block just sets up the rank amd size of the simulation run
    if run_type == "mpi":
        from mpi4py import MPI          # This module is imported here to be compatible with systems not having mpi support
        comm = MPI.COMM_WORLD
        size = comm.Get_size()          # Number of MPI processes
        rank = comm.Get_rank()          # Rank of this particular MPI process. Ranges from 0 to (size - 1)
    else:
        size = 1
        rank = 0

    # config_file is provided as a dot separated path to the config file with genesys as the parent directory
    config_file = sys.argv[2]
    # Importing the general simulation parameters
    config = importlib.import_module(config_file).config
    # Importing the scan strategy
    scan_strategy = importlib.import_module("genesys.scan_strategy." + config.scan_strategy_config).scan_strategy
    config.__dict__.update(scan_strategy.__dict__)
        
    # Getting the global list of detectors and segments
    num_segments = segd.count_segments(config.band_detector_segment_dict)
    detector_list = segd.get_detector_list(config.band_detector_segment_dict)
    total_simulation_time = segd.get_total_simulation_time(config.band_detector_segment_dict)
    # Getting the local list of detectors and segments
    local_band_detector_segment_dict = segd.get_local_band_detector_segment_dict(rank, size, config.band_detector_segment_dict)
    local_num_segments = segd.count_segments(local_band_detector_segment_dict)
    local_detector_list = segd.get_detector_list(local_band_detector_segment_dict)
    local_simulation_time = segd.get_total_simulation_time(local_band_detector_segment_dict)

    # Global display message
    if rank == 0:
        global_start_message()
    # Local display message
    local_start_message()
    # Waiting for all the processes to print out their respectiv start messagee
    if run_type == "mpi":
        comm.Barrier()

    # Making the parent directories for the sim
    if rank == 0:
        dio.make_top_data_directories(config, dir_list=['sim_dir', 'scan_dir'], verbose=True)

    # Displaying the parameters ofthe individual detectors
    if rank == 0:
        for band_name in list(local_band_detector_segment_dict.keys()):
            for detector_name in list(local_band_detector_segment_dict[band_name].keys()):
                detector = Detector(band_name, detector_name, config, sim_run=False)
                detector.display_params()
    # Waiting for all the processes to print out their respective detector parameters
    if run_type == "mpi":
        comm.Barrier()

    #  Running the main simulation routine. This is common to both mpi and serial runs.
    run_simulation()

    stop_time_main = time.time()
    prompt(colored("***RANK {} COMPLETE***. Total time taken : {}s\n", color="green").format(rank, stop_time_main - start_time_main))
    # Waiting for all the processes to finish
    if run_type == "mpi":
        comm.Barrier()
    if rank == 0:
        prompt(colored("***SIMULATION COMPLETE***\n", color="green"))
#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
