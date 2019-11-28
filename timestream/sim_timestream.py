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

import sys
import os
import time
import h5py
from termcolor import colored
from functools import partial
from genesys import Genesys_Class
from genesys import global_paths, load_param_file
from genesys.instruments.instrument import Instrument
from genesys.instruments.detector import Detector
from genesys.timestream.timestream import TStream
from genesys.data_io.segment_distribution import Data_Segment
from genesys.data_io.data_io import Data_IO
from genesys.numerical import unit_conversion as uc

#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
#* THE MASTER SECTION WHICH DISTRIBUTES THE DATA AND COLLECTS IT
#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*

def run_simulation():
    count = 0
    cumulative_time_taken = 0

    for channel_name in data_seg_local.channel_list():
        dio.set_path_for_channel(channel_name)
        dio.make_channel_directory()
        for detector_name in data_seg_local.detectors_in_channel(channel_name):
            detector = Detector(instrument_name=sim_config['instrument_name'], channel_name=channel_name, detector_name=detector_name)
            detector.initialise_pointing()
            detector.load_map(sim_config['sim_pol_type'])
            dio.set_path_for_detector(detector_name)
            dio.make_detector_directory()
            for segment in data_seg_local.channel_detector_segment_dict[channel_name][detector_name]:
                print("Rank {} doing segment {} of detector {} in channel {}".format(rank, segment, detector_name, channel_name))
                segment_start_time = time.time()
                count += 1
                ts = TStream(sim_config=sim_config, detector=detector, segment=segment+1)
                ts.generate_scan_tstream()
                write_tstream(ts.ts, dio.get_path_to_segment_file(segment))
                segment_stop_time = time.time()
                time_taken = segment_stop_time - segment_start_time
                cumulative_time_taken += time_taken
                print("Rank {} finished {} of {} segments.\nTime taken: {}s. Total time taken: {}, Projected time: {}\n".format(rank, count, data_seg_local.num_segments(), time_taken, cumulative_time_taken, cumulative_time_taken*data_seg_local.num_segments()/count))

def write_tstream(ts, file_name):
    print(list(ts.keys()))
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
    print(colored("\n#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#", color='blue'))
    print(colored("#*", color='blue') + 25*" " + colored("GENESYS", color='green') + 25*" " + colored("*#", color='blue'))
    print(colored("#*", color='blue') + 13*" " + colored("BEGINNING TIMESTREAM SIMULATION", color='green') + 13*" " + colored("*#", color='blue'))
    print(colored("#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#", color='blue'))
    print("RUN TYPE: {}".format(run_type))
    print("# OF PROCESSES: {}".format(size))
    print("SIMULATION TAG: {}".format(sim_config['sim_tag']))
    print("TOD TAG: {}".format(sim_config['tod_tag']))
    print("SPECIAL TAG: {}".format(sim_config['special_tag']))
    print("SIMULATION POLARISATION TYPE: {}".format(sim_config['sim_pol_type']))
    print("COORDINATE SYSTEM: {}".format(sim_config['coordinate_system']))
    print("CHANNEL LIST: {}".format(data_seg_global.channel_list()))
    print("# OF CHANNELS: {}".format(data_seg_global.num_channels()))
    print("DETECTOR LIST: {}".format(data_seg_global.formatted_detector_list()))
    print("# OF DETECTORS: {}".format(data_seg_global.num_detectors()))
    print("TOTAL NUMBER OF SEGMENTS : {}".format(data_seg_global.num_segments()))
    print("SEGMENT LENGTH: {}s = {}d = {}y".format(sim_config['segment_length'], uc.convert_unit('time', sim_config['segment_length'], 'second', 'day'), uc.convert_unit('time', sim_config['segment_length'], 'second', 'siderial year')))
    sim_time = data_seg_global.num_segments() * sim_config['segment_length']
    print("TOTAL SIMULATION TIME : {}s = {}d = {}y".format(sim_time, uc.convert_unit('time', sim_time, 'second', 'day'), uc.convert_unit('time', sim_time, 'second', 'siderial year')))
    sim_time = data_seg_global.num_segments() * sim_config['segment_length'] / data_seg_global.num_detectors()
    print("TOTAL SIMULATION TIME PER DETECTOR: {}s = {}d = {}y".format(sim_time, uc.convert_unit('time', sim_time, 'second', 'day'), uc.convert_unit('time', sim_time, 'second', 'siderial year')))
    print(colored("#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*\n", color='blue'))

def local_start_message():
    print(colored("#* RANK: {} #*#*#*#*#*#*#*#*#*#*".format(rank), color="green"))
    print("CHANNEL LIST: {}".format(data_seg_local.channel_list()))
    print("# OF CHANNELS: {}".format(data_seg_local.num_channels()))
    print("DETECTOR LIST: {}".format(data_seg_local.formatted_detector_list()))
    print("# OF DETECTORS: {}".format(data_seg_local.num_detectors()))
    print("LOCAL NUMBER OF SEGMENTS : {}".format(data_seg_local.num_segments()))
    sim_time = data_seg_local.num_segments() * sim_config['segment_length']
    print("LOCAL SIMULATION TIME : {}s = {}d = {}y".format(sim_time, uc.convert_unit('time', sim_time, 'second', 'day'), uc.convert_unit('time', sim_time, 'second', 'siderial year')))
    print(colored("#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*\n", color="green"))

#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
#* MAIN FUNCTION DEFINITION. THIS IS WHERE THE CODE BEGINS WHEN EXECUTED
#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*

if __name__=="__main__":
    start_time_main = time.time()

    #*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
    # THIS BLOCK JUST SETS UP THE RANK AMD SIZE OF THE SIMULATION RUN
    #*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
    run_type = sys.argv[1]
    # RUN_TYPE CAN BE EITHER "mpi" OR "serial"
    assert run_type in ["mpi", "serial"], "Please enter \"mpi\" OR \"serial\" as the second argument."
    # NOW SETTING UP THE rank AND size
    if run_type == "mpi":
        from mpi4py import MPI          # This module is imported here to be compatible with systems not having mpi support
        comm = MPI.COMM_WORLD
        size = comm.Get_size()          # Number of MPI processes
        rank = comm.Get_rank()          # Rank of this particular MPI process. Ranges from 0 to (size - 1)
    else:
        size = 1
        rank = 0
    #*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*

    #*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
    # THIS BLOCK LOADS THE CONFIGURATION FILE AND INITIALISES THE INSTRUMENT OBJECT
    #*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
    # sim_config_file_name IS LOCATED IN config_files DIRECTORY
    sim_config_file_name = sys.argv[2]
    current_dir = os.path.dirname(__file__)
    sim_config_file = os.path.join(current_dir,'config_files', sim_config_file_name)
    sim_config = load_param_file(file_path=sim_config_file)
    verbosity = sim_config['verbosity']
    #*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
        
    #*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
    # THIS BLOCK GETS THE LOCAL AND GLOBAL DICTIONARY OF BANDS, DETECTORS AND SEGMENTS
    #*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
    data_seg_global = Data_Segment(0, 1, sim_config['channel_detector_dict'], sim_config['num_segments_per_det'])
    data_seg_local = Data_Segment(rank, size, sim_config['channel_detector_dict'], sim_config['num_segments_per_det'])
    #*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*

    #*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
    # THIS BLOCK LOADS THE INSTRUMENT AND THE LOCAL DETECTOR OBJECTS
    #*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
    # Loading the instrument object
    instrument = Instrument(instrument_name=sim_config['instrument_name'])
    #*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*

    #*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
    # THIS BLOCK DOES THE OPENING DISPLAY MESSAGES
    #*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
    # GLOBAL DISPLAY MESSAGE
    if verbosity >= 1:
        if rank == 0:
            global_start_message()
        # LOCAL DISPLAY MESSAGE
        local_start_message()
        sys.stdout.flush()
    #*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*

    #*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
    # THIS BLOCK MAKES THE PARENT DATA DIRECTORIES
    #*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
    dio = Data_IO(sim_config)
    if rank == 0:
        dio.make_top_data_directories(dir_list=['sim_dir', 'tod_dir'])
    #*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*

    #*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
    # THIS BLOCK DOES THE INSTRUMENT DISPLAY MESSAGES 
    #*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
    if verbosity >= 2:
        if rank == 0:
            print(colored("\n#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#", color='blue'))
            print(colored("#*", color='blue') + 18*" " + colored("INSTRUMENT PARAMETERS", color='green') + 18*" " + colored("*#", color='blue'))
            print(colored("#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#", color='blue'))
            instrument.info()
            print(colored("#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#", color='blue'))
        sys.stdout.flush()
    #*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*

    #*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
    # THIS BLOCK DOES THE DETECTOR DISPLAY MESSAGES 
    #*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
    if verbosity >= 2:
        if rank == 0:
            print(colored("\n#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#", color='blue'))
            print(colored("#*", color='blue') + 19*" " + colored("DETECTOR PARAMETERS", color='green') + 19*" " + colored("*#", color='blue'))
            print(colored("#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#", color='blue'))
            for channel_name in data_seg_global.channel_list():
                for detector_name in data_seg_global.detectors_in_channel(channel_name):
                    detector = Detector(instrument_name=sim_config['instrument_name'], channel_name=channel_name, detector_name=detector_name)
                    detector.info()
            print(colored("#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#", color='blue'))
        sys.stdout.flush()

    #  Running the main simulation routine. This is common to both mpi and serial runs.
    run_simulation()

    stop_time_main = time.time()
    print(colored("***RANK {} COMPLETE***. Total time taken : {}s\n", color="green").format(rank, stop_time_main - start_time_main), flush=True)
    # Waiting for all the processes to finish
    if run_type == "mpi":
        comm.Barrier()
    if rank == 0:
        print(colored("***SIMULATION COMPLETE***\n", color="green"))
#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
