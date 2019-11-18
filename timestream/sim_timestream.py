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
from termcolor import colored
from genesys import Genesys_Class
from genesys import global_paths, load_param_file, prompt
from genesys.instruments.instrument import Instrument
from genesys.data_io import segment_distribution as segd
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

    for channel_name in list(local_channel_detector_segment_dict.keys()):
        for detector_name in list(local_channel_detector_segment_dict[channel_name].keys()):
            detector = Detector(channel_name, detector_name, config)
            for segment in local_channel_detector_segment_dict[channel_name][detector_name]:
                segment_start_time = time.time()
                count += 1
                detector.segment_start_prompt(rank, segment)
                detector.simulate_timestream_data(segment)
                segment_stop_time = time.time()
                time_taken = segment_stop_time - segment_start_time
                cumulative_time_taken += time_taken
                prompt("Rank {} finished {} of {} segments.\nTime taken: {}s. Total time taken: {}, Projected time: {}\n".format(rank, count, local_num_segments, time_taken, cumulative_time_taken, cumulative_time_taken*local_num_segments/count))

#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
#* THE START MESSAGES
#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*

def global_start_message():
    display_string = colored("\n#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#\n", color='blue')
    display_string += colored("#*", color='blue') + 25*" " + colored("GENESYS", color='green') + 25*" " + colored("*#\n", color='blue')
    display_string += colored("#*", color='blue') + 13*" " + colored("BEGINNING TIMESTREAM SIMULATION", color='green') + 13*" " + colored("*#\n", color='blue')
    display_string += colored("#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#\n", color='blue')
    display_string += "RUN TYPE: {}\n".format(run_type)
    display_string += "# OF PROCESSES: {}\n".format(size)
    display_string += "SIMULATION TAG: {}\n".format(sim_config['sim_tag'])
    display_string += "TOD TAG: {}\n".format(sim_config['scan_tag'])
    display_string += "SPECIAL TAG: {}\n".format(sim_config['special_tag'])
    display_string += "SIMULATION POLARISATION TYPE: {}\n".format(sim_config['sim_pol_type'])
    display_string += "COORDINATE SYSTEM: {}\n".format(sim_config['coordinate_system'])
    display_string += "CHANNEL LIST: {}\n".format(channel_list_global)
    display_string += "# OF CHANNELS: {}\n".format(len(channel_list_global))
    display_string += "DETECTOR LIST: {}\n".format(formatted_detector_list_global)
    display_string += "# OF DETECTORS: {}\n".format(num_detectors_global)
    display_string += "TOTAL NUMBER OF SEGMENTS : {}\n".format(num_segments_global)
    display_string += "TOTAL SIMULATION TIME : {}s {}d {}y\n".format(simulation_time_global, uc.convert_unit('time', simulation_time_global, 'second', 'day'), uc.convert_unit('time', simulation_time_global, 'second', 'siderial year'))
    #  display_string += "POLARISATION MODULATION TYPE: {}\n".format(sim_config['polarisation_modulation'])
    #  if config.polarisation_modulation == "continuous_HWP":
        #  display_string += "HWP initial phase: {} degrees\n".format(sim_config['HWP_phase_ini'])
        #  display_string += "HWP angular speed: {} RPM\n".format(config.HWP_rpm)
    #  if config.polarisation_modulation == "stepped_HWP":
        #  display_string += "HWP initial phase: {} degrees\n".format(config.HWP_phase_ini)
        #  display_string += "HWP step size: {} degrees\n".format(config.HWP_step)
        #  display_string += "HWP step duration: {} s\n".format(config.HWP_step_duration)
    display_string += colored("#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*\n\n", color='blue')
    prompt(display_string)

def local_start_message():
    display_string = colored("#* RANK: {} #*#*#*#*#*#*#*#*#*#*\n".format(rank), color="green")
    display_string += "CHANNEL LIST: {}\n".format(channel_list_local)
    display_string += "# OF CHANNELS: {}\n".format(len(channel_list_local))
    display_string += "DETECTOR LIST: {}\n".format(formatted_detector_list_local)
    display_string += "# OF DETECTORS: {}\n".format(num_detectors_local)
    display_string += "LOCAL NUMBER OF SEGMENTS: {}\n".format(num_segments_local)
    display_string += "TOTAL LOCAL SIMULATION TIME : {}s {}d {}y\n".format(simulation_time_local, uc.convert_unit('time', simulation_time_local, 'second', 'day'), uc.convert_unit('time', simulation_time_local, 'second', 'siderial year'))
    display_string += colored("#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*\n\n", color="green")
    prompt(display_string)

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
    #*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
        
    #*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
    # THIS BLOCK GETS THE LOCAL AND GLOBAL DICTIONARY OF BANDS, DETECTORS AND SEGMENTS
    #*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
    # GETTING THE GLOBAL LIST OF DETECTORS AND SEGMENTS
    channel_detector_segment_dict_global = segd.get_local_channel_detector_segment_dict(0, 1, sim_config['channel_detector_dict'], sim_config['num_segments_per_det'])
    channel_list_global = segd.get_channel_list(channel_detector_segment_dict_global)
    formatted_detector_list_global = segd.get_channel_list(channel_detector_segment_dict_global)
    num_detectors_global = segd.count_detectors(channel_detector_segment_dict_global)
    num_segments_global = segd.count_segments(channel_detector_segment_dict_global)
    simulation_time_global = num_segments_global * sim_config['segment_length']
    # GETTING THE LOCAL LIST OF DETECTORS AND SEGMENTS
    channel_detector_segment_dict_local = segd.get_local_channel_detector_segment_dict(rank, size, sim_config['channel_detector_dict'], sim_config['num_segments_per_det'])
    channel_list_local = segd.get_channel_list(channel_detector_segment_dict_local)
    formatted_detector_list_local = segd.get_channel_list(channel_detector_segment_dict_local)
    num_detectors_local = segd.count_detectors(channel_detector_segment_dict_local)
    num_segments_local = segd.count_segments(channel_detector_segment_dict_global)
    num_segments_local = segd.count_segments(channel_detector_segment_dict_local)
    simulation_time_local = num_segments_local * sim_config['segment_length']
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
    if rank == 0:
        global_start_message()
    # LOCAL DISPLAY MESSAGE
    display_string = colored("\n#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#\n", color='blue')
    display_string += colored("#*", color='blue') + 16*" " + colored("PROCESS DATA DISTRIBUTION", color='green') + 16*" " + colored("*#\n", color='blue')
    display_string += colored("#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#\n", color='blue')
    local_start_message()
    display_string += colored("#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#\n", color='blue')
    # WAITING FOR ALL THE PROCESSES TO PRINT OUT THEIR RESPECTIV START MESSAGEE
    if run_type == "mpi":
        comm.Barrier()
    #*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*

    #*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
    # THIS BLOCK MAKES THE PARENT DATA DIRECTORIES
    #*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
    if rank == 0:
        global_dio = Data_IO(sim_params)
        global_dio.make_top_data_directories(dir_list=['sim_dir', 'tod_dir'])
    #*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*

    #*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
    # THIS BLOCK DOES THE INSTRUMENT DISPLAY MESSAGES 
    #*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
    if rank == 0:
        display_string = colored("\n#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#\n", color='blue')
        display_string += colored("#*", color='blue') + 18*" " + colored("INSTRUMENT PARAMETERS", color='green') + 18*" " + colored("*#\n", color='blue')
        display_string += colored("#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#\n", color='blue')
        instrument.display_scan_strategy()
        instrument.display_half_wave_plate()
        instrument.display_channels()
        display_string += colored("#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#\n", color='blue')
        #  for band_name in list(local_band_detector_segment_dict.keys()):
            #  for detector_name in list(local_band_detector_segment_dict[band_name].keys()):
                #  detector = Detector(band_name, detector_name, config, sim_run=False)
                #  detector.display_params()
    #  # Waiting for all the processes to print out their respective detector parameters
    #  if run_type == "mpi":
        #  comm.Barrier()
#
    #  #  Running the main simulation routine. This is common to both mpi and serial runs.
    #  run_simulation()
#
    #  stop_time_main = time.time()
    #  prompt(colored("***RANK {} COMPLETE***. Total time taken : {}s\n", color="green").format(rank, stop_time_main - start_time_main))
    #  # Waiting for all the processes to finish
    #  if run_type == "mpi":
        #  comm.Barrier()
    #  if rank == 0:
        #  prompt(colored("***SIMULATION COMPLETE***\n", color="green"))
#  #*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
#  #*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
