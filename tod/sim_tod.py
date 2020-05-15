"""
This is the master executable to simulate TOD.
Data distribution is performed here.
This code is embarassingly parallel, that is, there is no explicit communication between any processes.
There are two options for running this:
    1) serial: a single process will iterate over all detectors and data segments.
    2) mpi: the different detectors and data segments will be distributed to different processes.
How to execute this file:
python sim_tod.py <config_file> <run_type> -v <verbosity>
where, config_file is the name of the config file to be passes to the code, run_type is either 'serial' or 'mpi', and verbosity is self explanatory
"""

import argparse
import os
import time
import numpy as np
from termcolor import colored
from genesys import load_param_file, prompt
from genesys.gen_io import GenIO
from genesys.instruments import Instrument
from genesys.pointing import Pointing
from genesys.numerical import Unit_Converter
from genesys.gen_io.segment_distribution import Segment_Distributor
from genesys.maps.sky_map import Sky_Map

#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
#* The simulation loop through all detectors
#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
def run_simulation():
    """
    The actual simulation is performed here.
    By convention, observation starts at segment #1. 
    Only segments are distributed among mpi processes.
    Each process therefore loops over all channels and detectors
    The order of the loops are:
        segment -> channel -> detector
    """
    pointing_obj = Pointing(instrument.params['scan_strategy'])

    input_maps_dict = load_maps()
    for data_block in sd.data_block_list_local:
        prompt(f"Rank {rank:^6} starting data block {data_block}", nature='info')
        for channel_name in config['channel_detector_dict'].keys():
            prompt(f"Rank {rank:^6} \tobserving with channel {channel_name}", nature='info')
            channel = instrument.get_channel(channel_name)
            io_obj.open_tod_file(channel_name, data_block, 'w')
            pol_ang = np.radians([channel.params['detectors'][det_name]['pol_phase_ini'] for det_name in config['channel_detector_dict'][channel_name]])
            det_list = np.string_(', '.join(config['channel_detector_dict'][channel_name]))
            mbangs = np.zeros(len(config['channel_detector_dict'][channel_name]))
            channel_commons = {'fsamp': channel.params['sampling_rate'], 'hwp_rpm': channel.params['HWP']['spin_rate'] , 'segl': config['segment_length'], 'det': det_list, 'polang': pol_ang, 'coords': config['coordinate_system'], 'pol_mod': channel.params['pol_modulation'], 'num_segments': config['num_segments_per_data_block'], 'mbangs': mbangs}
            io_obj.write_channel_common(channel_commons)
            segment_list = sd.get_segment_list(data_block, config['num_segments_per_data_block'])
            for segment in segment_list: 
                t_start = get_t_start(segment, config['segment_length'])
                pointing_obj.generate_obv_quaternion(t_start, channel.params['sampling_rate'], config['segment_length'], config['coordinate_system'])
                tod = {}
                if channel.params['pol_modulation'] != "passive":
                    tod['psi_hwp'] = pointing_obj.get_hwp_phase(t_start, channel.params['sampling_rate'], config['segment_length'], channel.params['HWP']['spin_rate'])
                segment_commons = {'time_0': t_start}
                io_obj.write_segment_common(segment, segment_commons)
                for detector_name in config['channel_detector_dict'][channel_name]:
                    detector = channel.get_detector(detector_name)
                    detector.observe_sky(pointing_obj, input_maps_dict[channel_name], tod, config['segment_length'])
                    scalars = {'gain': 1.0, 'sigma0': detector.params['noise']['white_noise_rms'], 'fknee': detector.params['noise']['f_knee'], 'alpha': detector.params['noise']['noise_alpha']}
                    io_obj.write_detector_data(segment, detector_name, tod, scalars)
            io_obj.close_tod_file()

#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
#* Messages and helper scripts 
#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*

def get_t_start(segment, segment_length):
    t_start = (segment - 1) * segment_length
    return t_start

def load_maps():
    """
    This hack is to take advantage of the massive memory on the owl system
    Loading all the input maps and keeping them in a dictionary
    """
    input_maps_dict = {}
    for channel_name in config['channel_detector_dict'].keys():
        map_file_name = instrument.params['channels'][channel_name]['input_map_file']
        input_maps_dict[channel_name] = Sky_Map(map_file_name, field=(0,1,2))
    return input_maps_dict

uc = Unit_Converter()
conv_u = uc.convert_unit

def global_start_message():
    print_header()
    print_header(text='GENESYS')
    print_header(text='SIMULATING SATELLITE TOD')
    print_header()
    prompt(f"{'RUN TYPE':<35}{in_args.run_type}")
    prompt(f"{'# OF PROCESSES':<35}{size}")
    prompt(f"{'SIMULATION TAG':<35}{config['sim_tag']}")
    prompt(f"{'SPECIAL TAG':<35}{config['special_tag']}")
    prompt(f"{'SIMULATION POLARISATION TYPE':<35}{config['pol_type']}")
    prompt(f"{'COORDINATE SYSTEM':<35}{config['coordinate_system']}")
    prompt(f"{'DETECTOR LIST':<35}{config['channel_detector_dict']}")
    prompt(f"{'# of DETECTORS':<35}{num_dets}")
    seg_len = config['segment_length']
    prompt(f"{'SEGMENT LENGTH':<35}{seg_len}sec = {conv_u(seg_len, 'time', 'sec', 'hour'):.2f}h = {conv_u(seg_len, 'time', 'sec', 'day'):.4f}d")
    prompt(f"{'NUMBER OF DATA BLOCKS':<35}{sd.num_data_blocks_global}")
    prompt(f"{'NUMBER OF SEGMENTS PER DATA BLOCK':<35}{config['num_segments_per_data_block']}")
    obv_len = sd.observation_time_global
    prompt(f"{'TOTAL OBSERVATION TIME':<35}{obv_len}sec = {conv_u(obv_len, 'time', 'sec', 'day'):.2f}d = {conv_u(obv_len, 'time', 'sec', 'year'):.4f}y")
    print_header()

def local_start_message():
    print_header(text=f'Rank {rank}', text_align='left')
    prompt(f"{'LOCAL DATA BLOCKS':<35}{sd.data_block_list_local}")
    prompt(f"{'LOCAL NUMBER OF DATA BLOCKS':<35}{sd.num_data_blocks_local}")
    obv_len = sd.observation_time_local
    prompt(f"{'OBSERVATION TIME':<35}{obv_len}sec = {conv_u(obv_len, 'time', 'sec', 'day'):.2f}d = {conv_u(obv_len, 'time', 'sec', 'year'):.2f}y")
    print_header()

def instrument_start_message():
    print_header()
    print_header(text='INSTRUMENT PARAMETERS')
    print_header()
    instrument.info()
    print_header()

def detector_start_message():
    print_header()
    print_header(text='DETECTOR PARAMETERS')
    print_header()
    for channel_name in config['channel_detector_dict'].keys():
        channel = instrument.get_channel(channel_name)
        for detector_name in config['channel_detector_dict'][channel_name]:
            detector = channel.get_detector(detector_name)
            detector.info()
    print_header()

def print_header(text=None, text_colour=None, length=81, text_align='center'):
    # default length = 81
    if text == None:
        prompt(40 * colored('#*', 'green') + colored('#', 'green'))
    else:
        if text_align == 'left':
            prompt(colored('# ', 'green') + colored(f'{text:<78}', text_colour) + colored('#', 'green'))
        else:
            prompt(colored('# ', 'green') + colored(f'{text:^78}', text_colour) + colored('#', 'green'))

#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
#* Main function definition. This is where the code begins when executed
#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*

if __name__=="__main__":
    start_time_main = time.time()

    #*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
    # Parsing arguments
    #*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
    parser = argparse.ArgumentParser()
    parser.add_argument('config_file', type=str, help='Name of config file under relative path config_files/')
    parser.add_argument('run_type', type=str, choices=['mpi','serial'], help='Type of run')
    parser.add_argument('--verbosity', '-v', type=int, choices=[0,1,2], default=0, help='Verbosity of the code')
    in_args = parser.parse_args()
    #*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*

    #*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
    # Sets up the rank amd size of the simulation run
    #*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
    # Setting up the rank and size
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
    # Setting up the config, instrument, detector_segments and data directories
    #*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
    # Sim config
    current_dir = os.path.dirname(__file__)
    config_file = os.path.join(current_dir,'config_files', in_args.config_file)
    config = load_param_file(file_path=config_file)
    # segment_list can be given as an explicit list, range(), np.arange or some other convenient way
    if type(config['data_block_list']) != list:
            config['data_block_list'] = list(eval(config['data_block_list']))
    # Instrument object
    instrument = Instrument(instrument_dir=config['instrument_dir'])
    # Segment distributor
    sd = Segment_Distributor(config['data_block_list'], config['num_segments_per_data_block'], size, rank)
    sd.set_observation_times(config['segment_length'])
    num_dets = sd.count_detectors(config['channel_detector_dict'])
    # I/O handling object
    io_obj = GenIO(config)
    if rank == 0:
        io_obj.make_top_data_directories(dir_list=['sim_dir', 'tod_dir'])
    #*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*

    #*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
    # The opening display messages
    #*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
    if rank == 0:
        global_start_message()
        instrument_start_message()
        detector_start_message()
    if in_args.run_type == "mpi":
        comm.Barrier()
        time.sleep(rank*0.1)
        local_start_message()
        comm.Barrier()
    #*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*

    # Everything is set up now for the simulation run.
    prompt(f"*** Rank {rank:^6} BEGINNING OBSERVATION ***", nature='info')
    run_simulation()
    stop_time_main = time.time()
    prompt(f"*** RANK {rank:^6} COMPLETE ***          Time taken: {stop_time_main - start_time_main:.2f}s", nature='info')

    # Waiting for all the processes to finish
    if in_args.run_type == "mpi":
        comm.Barrier()
    if rank == 0:
        print_header()
        print_header(text='SIMULATION COMPLETE')
        print_header()
#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
