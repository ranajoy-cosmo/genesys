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
#  from genesys.timestream.timestream import TStream
from genesys.data_io.segment_distribution import Data_Segment
from genesys.data_io.data_io import Data_IO
from genesys.numerical.unit_conversion import Unit_Converter

#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
#* THE MASTER SECTION WHICH DISTRIBUTES THE DATA AND COLLECTS IT
#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*

def run_simulation():
    count = 0
    cumulative_time_taken = 0

    for channel_name in data_seg_local.channel_list():
        channel = instrument.get_channel_object(channel_name)
        dio.set_path_for_channel(channel_name)
        dio.make_channel_directory()
        for detector_name in data_seg_local.detectors_in_channel(channel_name):
            detector = channel.get_detector_object(detector_name)
            detector.initialise_pointing()
            detector.load_map(sim_config['sim_pol_type'])
            dio.set_path_for_detector(detector_name)
            dio.make_detector_directory()
            for segment in data_seg_local.channel_detector_segment_dict[channel_name][detector_name]:
                prompt("Rank {} doing segment {} of detector {} in channel {}".format(rank, segment, detector_name, channel_name))
                segment_start_time = time.time()
                count += 1
                detector.initialise_timestream_segment(sim_config, segment)
                detector.observe_sky()
                write_tstream(detector.ts.t_stream, dio.get_path_to_segment_file(segment))
                segment_stop_time = time.time()
                time_taken = segment_stop_time - segment_start_time
                cumulative_time_taken += time_taken
                prompt("Rank {} finished {} of {} segments.\nTime taken: {}s. Total time taken: {}, Projected time: {}\n".format(rank, count, data_seg_local.num_segments(), time_taken, cumulative_time_taken, cumulative_time_taken*data_seg_local.num_segments()/count))

def write_tstream(ts, file_name):
    prompt(list(ts.keys()))
    f = h5py.File(file_name, 'a')
    for item in list(ts.keys()):
        f.create_dataset(item, data=ts[item])
    f.close()
#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
#* THE START MESSAGES
#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*

def global_start_message():
    uc = Unit_Converter()
    prompt(colored("\n#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#", color='blue'))
    prompt(colored("#*", color='blue') + 25*" " + colored("GENESYS", color='green') + 25*" " + colored("*#", color='blue'))
    prompt(colored("#*", color='blue') + 13*" " + colored("BEGINNING TIMESTREAM SIMULATION", color='green') + 13*" " + colored("*#", color='blue'))
    prompt(colored("#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#", color='blue'))
    prompt("RUN TYPE: {}".format(in_args.run_type))
    prompt("# OF PROCESSES: {}".format(size))
    prompt("SIMULATION TAG: {}".format(sim_config['sim_tag']))
    prompt("TOD TAG: {}".format(sim_config['tod_tag']))
    prompt("SPECIAL TAG: {}".format(sim_config['special_tag']))
    prompt("SIMULATION POLARISATION TYPE: {}".format(sim_config['sim_pol_type']))
    prompt("COORDINATE SYSTEM: {}".format(sim_config['coordinate_system']))
    prompt("CHANNEL LIST: {}".format(data_seg_global.channel_list()))
    prompt("# OF CHANNELS: {}".format(data_seg_global.num_channels()))
    prompt("DETECTOR LIST: {}".format(data_seg_global.formatted_detector_list()))
    prompt("# OF DETECTORS: {}".format(data_seg_global.num_detectors()))
    prompt("TOTAL NUMBER OF SEGMENTS : {}".format(data_seg_global.num_segments()))
    prompt("SEGMENT LENGTH: {}s = {}d = {}y".format(sim_config['segment_length'], uc.convert_unit(sim_config['segment_length'], 'time', 'sec', 'day'), uc.convert_unit(sim_config['segment_length'], 'time', 'sec', 'siderial year')))
    sim_time = data_seg_global.num_segments() * sim_config['segment_length']
    prompt("TOTAL SIMULATION TIME : {}s = {}d = {}y".format(sim_time, uc.convert_unit(sim_time, 'time', 'sec', 'day'), uc.convert_unit(sim_time, 'time', 'sec', 'siderial year')))
    sim_time = data_seg_global.num_segments() * sim_config['segment_length'] / data_seg_global.num_detectors()
    prompt("TOTAL SIMULATION TIME PER DETECTOR: {}s = {}d = {}y".format(sim_time, uc.convert_unit(sim_time, 'time', 'sec', 'day'), uc.convert_unit(sim_time, 'time', 'sec', 'siderial year')))
    prompt(colored("#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*\n", color='blue'))

def local_start_message():
    uc = Unit_Converter()
    prompt(colored("#* RANK: {} #*#*#*#*#*#*#*#*#*#*".format(rank), color="green"))
    prompt("CHANNEL LIST: {}".format(data_seg_local.channel_list()))
    prompt("# OF CHANNELS: {}".format(data_seg_local.num_channels()))
    prompt("DETECTOR LIST: {}".format(data_seg_local.formatted_detector_list()))
    prompt("# OF DETECTORS: {}".format(data_seg_local.num_detectors()))
    prompt("LOCAL NUMBER OF SEGMENTS : {}".format(data_seg_local.num_segments()))
    sim_time = data_seg_local.num_segments() * sim_config['segment_length']
    prompt("LOCAL SIMULATION TIME : {}s = {}d = {}y".format(sim_time, uc.convert_unit(sim_time, 'time', 'sec', 'day'), uc.convert_unit(sim_time, 'time', 'sec', 'siderial year')))
    prompt(colored("#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*\n", color="green"))

#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
#* MAIN FUNCTION DEFINITION. THIS IS WHERE THE CODE BEGINS WHEN EXECUTED
#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*

if __name__=="__main__":
    start_time_main = time.time()

    parser = argparse.ArgumentParser()
    parser.add_argument('config_file', type=str, help='Name of config file under relative path config_files/')
    parser.add_argument('run_type', type=str, choices=['mpi','serial'], help='Type of run')
    parser.add_argument('--verbosity', '-v', type=int, choices=[0,1,2], default=0, help='Verbosity of the code')
    in_args = parser.parse_args()

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
            prompt(colored("\n#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#", color='blue'))
            prompt(colored("#*", color='blue') + 18*" " + colored("INSTRUMENT PARAMETERS", color='green') + 18*" " + colored("*#", color='blue'))
            prompt(colored("#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#", color='blue'))
            instrument.info()
            prompt(colored("#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#", color='blue'))
    #*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*

    #*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
    # THIS BLOCK DOES THE DETECTOR DISPLAY MESSAGES
    #*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
    if in_args.verbosity == 2:
        if rank == 0:
            prompt(colored("\n#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#", color='blue'))
            prompt(colored("#*", color='blue') + 19*" " + colored("DETECTOR PARAMETERS", color='green') + 19*" " + colored("*#", color='blue'))
            prompt(colored("#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#", color='blue'))
            for channel_name in data_seg_global.channel_list():
                channel = instrument.get_channel(channel_name)
                for detector_name in data_seg_global.detectors_in_channel(channel_name):
                    detector = channel.get_detector(detector_name)
                    detector.info()
            prompt(colored("#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#", color='blue'))
    #*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*

    #  #  Running the main simulation routine. This is common to both mpi and serial runs.
    #  run_simulation()
#
    #  stop_time_main = time.time()
    #  prompt(colored("***RANK {} COMPLETE***. Total time taken : {}s\n", color="green").format(rank, stop_time_main - start_time_main), flush=True)
    #  # Waiting for all the processes to finish
    #  if run_type == "mpi":
        #  comm.Barrier()
    #  if rank == 0:
        #  prompt(colored("***SIMULATION COMPLETE***\n", color="green"))
#  #*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
#  #*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
