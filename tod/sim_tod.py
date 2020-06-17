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
    tod = {}

    block_count = 0
    det_count = 0
    proc_time_start = time.time()
    num_tot_det = sd.count_detectors(config['channel_detector_dict'])

    for data_block in sd.data_block_list_local:
        block_count += 1
        prompt(f"Rank {rank:^6} starting data block {data_block}. ({block_count}/{sd.num_data_blocks_local})", nature='info')
        segment_list = sd.get_segment_list(data_block, config['num_segments_per_data_block'])
        if config['add_orbital_dipole']:
            t_start_block = get_t_start(segment_list[0], config['segment_length'])
            sat_vel_dict = get_satellite_velocity(t_start_block, config['segment_length'], segment_list)
        else:
            sat_vel_dict = dict.fromkeys(segment_list, [0.0,0.0,0.0])
        for channel_name in config['channel_detector_dict'].keys():
            if in_args.verbosity == 2:
                chan_time_start = time.time()
                det_count += len(config['channel_detector_dict'][channel_name])
            prompt(f"Rank {rank:^6} \tobserving with channel {channel_name}", nature='info')
            channel = instrument.get_channel(channel_name)
            io_obj.open_tod_file(channel_name, data_block, 'w')
            channel_common, channel_common_attributes = get_channel_common(channel.params, channel_name)
            io_obj.write_channel_common(channel_common, channel_common_attributes)
            for segment in segment_list: 
                t_start = get_t_start(segment, config['segment_length'])
                segment_common = {'time': t_start, 'vsun': sat_vel_dict[segment]}
                segment_common_attributes = {'vsun': {'info': '[x,y,z]', 'coords': config['coordinate_system']}, 'time': {'type': 'second'}}
                if channel.params['pol_modulation'] != "passive":
                    tod['psi_hwp'] = pointing_obj.get_hwp_phase(t_start=t_start, fsamp=channel.params['sampling_rate'], seg_len=config['segment_length'], hwp_spin_rate=channel.params['HWP']['spin_rate'])
                    segment_common['psi_hwp'] = tod['psi_hwp'][0]
                pointing_obj.generate_obv_quaternion(t_start=t_start, fsamp=channel.params['sampling_rate'], seg_len=config['segment_length'], coords=config['coordinate_system'])
                io_obj.write_segment_common(segment, segment_common, segment_common_attributes)
                for detector_name in config['channel_detector_dict'][channel_name]:
                    detector = channel.get_detector(detector_name)
                    detector.observe_sky(pointing_obj, input_maps_dict[channel_name], tod, config['segment_length'], sat_vel_dict[segment], config['noise_type'], config['tod_write_field'])
                    scalars = {'gain': 1.0, 'sigma0': detector.params['noise']['white_noise_rms'], 'fknee': detector.params['noise']['f_knee'], 'alpha': detector.params['noise']['noise_alpha']}
                    io_obj.write_tod(segment, detector_name, tod, config['tod_write_field'])
                    io_obj.write_tod_scalars(segment, detector_name, scalars)
            io_obj.close_tod_file()
            if in_args.verbosity == 2:
                chan_time_stop = time.time()
                chan_time_taken = chan_time_stop - chan_time_start
                elapsed_time = chan_time_stop - proc_time_start
                estimated_time = elapsed_time * (sd.num_data_blocks_local * num_tot_det / det_count) 
                prompt(f"Rank {rank:^6} \ttook {chan_time_taken:.2f}s for {channel_name}. Elapsed time {elapsed_time:.2f}s. Estimated total time {estimated_time:.2f}s", nature='info')

#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
#* Helper scripts 
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
        map_file_name = os.path.join(config['map_directory'], config['channel_map_dict'][channel_name])
        input_maps_dict[channel_name] = Sky_Map(map_file_name, field=(0,1,2))
    return input_maps_dict

def get_satellite_velocity(t_start_block, segment_length, segment_list):
    time, xvel, yvel, zvel = instrument.get_satellite_velocity()
    vel_dict = {}
    if config['coordinate_system'] == "galactic":
        pt_obj = Pointing()
        rot_quat = pt_obj.quaternion_coordinate_transformation('E', 'G')
    for segment in segment_list:
        t_start = t_start_block + segment_length*(segment - segment_list[0])
        t_stop = t_start + segment_length
        time_slice = time[(time >= t_start) & (time <= t_stop)]
        if time_slice != []:
            time_selected = time_slice[0]
        else:
            time_selected = time[time < t_start][-1]
        vel_index = np.where(time == time_selected)[0][0]
        vel = np.array([xvel[vel_index], yvel[vel_index], zvel[vel_index]])
        if config['coordinate_system'] == "galactic":
            vel = pt_obj.rotate_vector(rot_quat, vel)
        vel_dict[segment] = vel
    return vel_dict

def get_channel_common(params, channel_name):
    polang = {det_name: params['detectors'][det_name]['pol_phase_ini'] for det_name in config['channel_detector_dict'][channel_name]}
    mbangs = dict.fromkeys(config['channel_detector_dict'][channel_name], 0.0)
    det = np.string_(', '.join(config['channel_detector_dict'][channel_name]))
    channel_common = {'fsamp': params['sampling_rate'], 'hwp_rpm': params['HWP']['spin_rate'] , 'segl': config['segment_length'], 'det': det, 'polang': list(polang.values()), 'coords': config['coordinate_system'], 'pol_mod': params['pol_modulation'], 'mbangs': list(mbangs.values())}
    channel_common_attributes = {'polang': {'legend': list(polang.keys())}, 'mbangs': {'legned': list(mbangs.keys())}} 
    return channel_common, channel_common_attributes

uc = Unit_Converter()
conv_u = uc.convert_unit

#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
#* Messages
#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*

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
    prompt(f"{'NOISE TYPE':<35}{config['noise_type']}")
    prompt(f"{'ORBITAL DIPOLE':<35}{config['add_orbital_dipole']}")
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
    channel_list = list(config['channel_detector_dict'].keys())
    instrument.info(channel_list)
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
    config['tod_write_field'] = config['tod_special_fields'] + ["signal", "theta", "phi", "psi"]
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
        if in_args.verbosity >=1:
            instrument_start_message()
        if in_args.verbosity >=2:
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
