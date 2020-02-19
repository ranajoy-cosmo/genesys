"""
This is the master executable for the timestream generation pipeline.
Data distribution is performed from here.
This code is embarassingly parallel. that is, there is no explicit communication between any processes.
There are two options for running this
1) Serial : A single process will iterate over all detectors and data segments.
2) MPI : The different detectors and data segments will be distributed to different processes.
How to execute this file:
python sim_timestream.py <CONFIG_FILE_PATH> <RUN_TYPE>
where, CONFIG_FILE_PATH is the absolute path to the config file to be passes to the code and RUN_TYPE is either 'run_serial' or 'run_mpi'
"""
#!/usr/bin/env python

import argparse
import os
import time
import h5py
from termcolor import colored
from genesys import Genesys_Class
from genesys import global_paths, load_param_file, prompt
from genesys.instruments.instrument import Instrument
from genesys.timestream.timestream import TStream
from genesys.data_io.segment_distribution import Data_Segment
from genesys.data_io.data_io import Data_IO
from genesys.numerical.unit_conversion import Unit_Converter

#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
#* THE MASTER SECTION WHICH DISTRIBUTES THE DATA AND COLLECTS IT
#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*

def run_simulation():
    prompt("")
    prompt(colored("#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#", color='green'))
    prompt(colored("#*", color='green') + 18*" " + colored("BEGINNING OBSERVATION", color='green') + 18*" " + colored("*#", color='green'))
    prompt(colored("#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#", color='green'))
    count = 0
    cumulative_time_taken = 0

    for channel_name in data_seg_local.channel_list():
        channel = instrument.get_channel(channel_name)
        dio.set_path_for_channel(channel_name)
        dio.make_channel_directory()
        for detector_name in data_seg_local.detectors_in_channel(channel_name):
            detector = channel.get_detector(detector_name)
            detector.load_map(sim_config['sim_pol_type'])
            dio.set_path_for_detector(detector_name)
            dio.make_detector_directory()
            for segment in data_seg_local.channel_detector_segment_dict[channel_name][detector_name]:
                if in_args.verbosity >= 1:
                    segment_start_time = time.time()
                    prompt(f"Rank {rank}: doing segment {segment} of detector {detector_name} in channel {channel_name}", nature='info')
                    count += 1
                t_stream = TStream()
                detector.prepare_for_observation(t_stream, sim_config['segment_length'], segment)
                detector.observe_sky(t_stream, sim_config['coordinate_system'], sim_config['sim_pol_type'])
                dio.write_t_stream_to_file(t_stream, segment, ['signal', 'theta', 'phi', 'psi', 'hwp_psi'])
                if in_args.verbosity >= 1:
                    segment_stop_time = time.time()
                    time_taken = segment_stop_time - segment_start_time
                    cumulative_time_taken += time_taken
                    projected_time = cumulative_time_taken*data_seg_local.num_segments()/count
                    prompt(f"Rank {rank}: finished {count}/{data_seg_local.num_segments()} segments. Time taken: {time_taken}s. Total time taken: {cumulative_time_taken}, Projected time: {projected_time}", nature='info')

#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
#* THE START MESSAGES
#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*

def global_start_message():
    uc = Unit_Converter()
    prompt(colored("#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#", color='green'))
    prompt(colored("#*", color='green') + 25*" " + colored("GENESYS", color='green') + 25*" " + colored("*#", color='green'))
    prompt(colored("#*", color='green') + 13*" " + colored("BEGINNING TIMESTREAM SIMULATION", color='green') + 13*" " + colored("*#", color='green'))
    prompt(colored("#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#", color='green'))
    prompt(f"RUN TYPE: {in_args.run_type}")
    prompt(f"# OF PROCESSES: {size}")
    prompt(f"SIMULATION TAG: {sim_config['sim_tag']}")
    prompt(f"TOD TAG: {sim_config['tod_tag']}")
    prompt(f"SPECIAL TAG: {sim_config['special_tag']}")
    prompt(f"SIMULATION POLARISATION TYPE: {sim_config['sim_pol_type']}")
    prompt(f"COORDINATE SYSTEM: {sim_config['coordinate_system']}")
    prompt(f"CHANNEL LIST: {data_seg_global.channel_list()}")
    prompt(f"# OF CHANNELS: {data_seg_global.num_channels()}")
    prompt(f"DETECTOR LIST: {data_seg_global.formatted_detector_list()}")
    prompt(f"# OF DETECTORS: {data_seg_global.num_detectors()}")
    prompt(f"TOTAL NUMBER OF SEGMENTS : {data_seg_global.num_segments()}")
    prompt(f"SEGMENT LENGTH: {sim_config['segment_length']}s = {uc.convert_unit(sim_config['segment_length'], 'time', 'sec', 'day')}d = {uc.convert_unit(sim_config['segment_length'], 'time', 'sec', 'siderial year')}y")
    sim_time = data_seg_global.num_segments() * sim_config['segment_length']
    prompt(f"TOTAL SIMULATION TIME : {sim_time}s = {uc.convert_unit(sim_time, 'time', 'sec', 'day')}d = {uc.convert_unit(sim_time, 'time', 'sec', 'siderial year')}y")
    sim_time = data_seg_global.num_segments() * sim_config['segment_length'] / data_seg_global.num_detectors()
    prompt(f"TOTAL SIMULATION TIME PER DETECTOR: {sim_time}s = {uc.convert_unit(sim_time, 'time', 'sec', 'day')}d = {uc.convert_unit(sim_time, 'time', 'sec', 'siderial year')}y")
    prompt(colored("#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*\n", color='green'))

def local_start_message():
    uc = Unit_Converter()
    prompt(colored(f"#* RANK: {rank} #*#*#*#*#*#*#*#*#*#*", color="green"))
    prompt(f"CHANNEL LIST: {data_seg_local.channel_list()}")
    prompt(f"# OF CHANNELS: {data_seg_local.num_channels()}")
    prompt(f"DETECTOR LIST: {data_seg_local.formatted_detector_list()}")
    prompt(f"# OF DETECTORS: {data_seg_local.num_detectors()}")
    prompt(f"LOCAL NUMBER OF SEGMENTS : {data_seg_local.num_segments()}")
    sim_time = data_seg_local.num_segments() * sim_config['segment_length']
    prompt(f"LOCAL SIMULATION TIME : {sim_time}s = {uc.convert_unit(sim_time, 'time', 'sec', 'day')}d = {uc.convert_unit(sim_time, 'time', 'sec', 'siderial year')}y")
    prompt(colored("#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*\n", color="green"))

#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
#* MAIN FUNCTION DEFINITION. THIS IS WHERE THE CODE BEGINS WHEN EXECUTED
#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*

if __name__=="__main__":
    start_time_main = time.time()

    #*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
    # PARSING ARGUMENTS
    #*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
    parser = argparse.ArgumentParser()
    parser.add_argument('config_file', type=str, help='Name of config file under relative path config_files/')
    parser.add_argument('run_type', type=str, choices=['mpi','serial'], help='Type of run')
    parser.add_argument('--verbosity', '-v', type=int, choices=[0,1,2], default=0, help='Verbosity of the code')
    in_args = parser.parse_args()
    #*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*

    #*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
    # SETS UP THE RANK AMD SIZE OF THE SIMULATION RUN
    #*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
    # SETTING UP THE rank AND size
    if in_args.run_type == "mpi":
        from mpi4py import MPI          # This module is imported here to be compatible with systems not having mpi support
        comm = MPI.COMM_WORLD
        size = comm.Get_size()          # Number of MPI processes
        rank = comm.Get_rank()          # Rank of this particular MPI process. Ranges from 0 to (size - 1)
    else:
        size = 1
        rank = 0
    #*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*

    #*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
    # SETTING UP THE CONFIG, INSTRUMENT, DETECTOR_SEGMENTS AND DATA DIRECTORIES
    #*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
    # SIM CONFIG
    current_dir = os.path.dirname(__file__)
    sim_config_file = os.path.join(current_dir,'config_files', in_args.config_file) # sim_config_file_name is located in config_files directory
    sim_config = load_param_file(file_path=sim_config_file)
    # LOADING THE INSTRUMENT OBJECT
    instrument = Instrument(instrument_dir=sim_config['instrument_name'])
    # GETTING THE LOCAL AND GLOBAL DETECTOR-SEGMENTS
    data_seg_global = Data_Segment(0, 1, sim_config['channel_detector_dict'], sim_config['num_segments_per_det'])
    data_seg_local = Data_Segment(rank, size, sim_config['channel_detector_dict'], sim_config['num_segments_per_det'])
    # DATA DIRECTORIES
    dio = Data_IO(sim_config)
    if rank == 0:
        dio.make_top_data_directories(dir_list=['sim_dir', 'tod_dir'])
    #*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*

    #*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
    # THE OPENING DISPLAY MESSAGES
    #*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
    # GLOBAL DISPLAY MESSAGE
    if in_args.verbosity >= 1:
        if rank == 0:
            global_start_message()
        # LOCAL DISPLAY MESSAGE
        local_start_message()
    #*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*

    #*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
    # THIS BLOCK DOES THE INSTRUMENT DISPLAY MESSAGES
    #*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
    if in_args.verbosity == 2:
        if rank == 0:
            prompt(colored("#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#", color='green'))
            prompt(colored("#*", color='green') + 18*" " + colored("INSTRUMENT PARAMETERS", color='green') + 18*" " + colored("*#", color='green'))
            prompt(colored("#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#", color='green'))
            instrument.info()
            prompt(colored("#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#", color='green'))
    #*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*

    #*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
    # THIS BLOCK DOES THE DETECTOR DISPLAY MESSAGES
    #*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
    if in_args.verbosity == 2:
        if rank == 0:
            prompt(colored("#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#", color='green'))
            prompt(colored("#*", color='green') + 19*" " + colored("DETECTOR PARAMETERS", color='green') + 19*" " + colored("*#", color='green'))
            prompt(colored("#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#", color='green'))
            for channel_name in data_seg_global.channel_list():
                channel = instrument.get_channel(channel_name)
                for detector_name in data_seg_global.detectors_in_channel(channel_name):
                    detector = channel.get_detector(detector_name)
                    detector.info()
            prompt(colored("#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#", color='green'))
    #*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*

    #  RUNNING THE MAIN SIMULATION ROUTINE. THIS IS COMMON TO BOTH MPI AND SERIAL RUNS.
    run_simulation()

    stop_time_main = time.time()
    prompt(colored(f"***RANK {rank} COMPLETE***. Total time taken : {stop_time_main - start_time_main}s", color="green"))
    # Waiting for all the processes to finish
    if in_args.run_type == "mpi":
        comm.Barrier()
    if rank == 0:
        prompt(colored("***SIMULATION COMPLETE***\n", color="green"))
#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
