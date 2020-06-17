import argparse
import os
import sys
import time
import numpy as np
import healpy as hp
from termcolor import colored
from genesys import load_param_file, prompt
from genesys.gen_io import GenIO
from genesys.pointing import Pointing
from genesys.gen_io.segment_distribution import Segment_Distributor
import covariance_matrix_utils as cov_ut
from genesys.maps.sky_map import Sky_Map

def run_map_maker():

    pointing_obj = Pointing()

    for channel_name in config['channel_detector_dict'].keys():
        prompt(f"Rank {rank:^6} starting channel {channel_name}", nature='info')
        npix = hp.nside2npix(config['nside'][channel_name])
        inv_cov_matrix_local, b_matrix_local, hitmap_local = cov_ut.get_empty_matrices(config['nside'][channel_name], config['pol_type'])
        for data_block in sd.data_block_list_local:
            prompt(f"Rank {rank:^6} \tstarting data block {data_block}", nature='info')
            io_obj.open_tod_file(channel_name, data_block, 'r')
            chan_com, chan_com_attr = io_obj.read_channel_common(['fsamp', 'hwp_rpm', 'segl', 'pol_mod'])
            segment_list = sd.get_segment_list(data_block, config['num_segments_per_data_block'])
            for segment in segment_list:
                tod = {}
                seg_com, seg_com_attr = io_obj.read_segment_common(segment, ['time', 'psi_hwp'])
                if chan_com['pol_mod'] == 'passive':
                    tod['psi_hwp'] = 0
                else:
                    tod['psi_hwp'] = seg_com['psi_hwp'] + pointing_obj.get_hwp_phase(0.0, chan_com['fsamp'], chan_com['segl'], chan_com['hwp_rpm'])
                for detector_name in config['channel_detector_dict'][channel_name]:
                    tod.update(io_obj.read_tod(segment, detector_name, ['signal', 'theta', 'phi', 'psi']))
                    hitpix = hp.ang2pix(config['nside'][channel_name], tod['theta'], tod['phi'])
                    cov_ut.add_to_mm_matrices(hitpix, tod['psi'], tod['psi_hwp'], tod['signal'], inv_cov_matrix_local, b_matrix_local, hitmap_local, npix, config['pol_type'])
            io_obj.close_tod_file()
        prompt(f"Rank {rank:^6} acquired channel {channel_name} data.", nature='info')
        if in_args.run_type == "mpi":
            # Distributing and gathering the segments of the matrices in the proper process
            # The sky pixels are chunked and are handled by individual processes
            inv_cov_matrix_local_block = distribute_matrix(inv_cov_matrix_local, "cov_matrix", config['nside'][channel_name])
            del inv_cov_matrix_local
            b_matrix_local_block = distribute_matrix(b_matrix_local, "b_matrix", config['nside'][channel_name])
            del b_matrix_local
            hitmap_local_block = distribute_matrix(hitmap_local, "hitmap", config['nside'][channel_name])
            del hitmap_local

            prompt(f"Rank {rank:^6} making covariance matrix for channel {channel_name}.", nature='info')
            cov_matrix_local_block = cov_ut.get_covariance_matrix(inv_cov_matrix_local_block, hitmap_local_block, config['pol_type'])

            prompt(f"Rank {rank:^6} estimating sky map for channel {channel_name}.", nature='info')
            sky_map_local_block = cov_ut.get_sky_map(cov_matrix_local_block, b_matrix_local_block, hitmap_local_block, config['pol_type'])
            cov_matrix_local_block_flat = cov_ut.flatten_block_matrix(cov_matrix_local_block, config['pol_type'])

            prompt(f"Rank {rank:^6} gathering maps for channel {channel_name}.", nature='info')
            sky_map = gather_map_segments(sky_map_local_block, 'sky_map', config['nside'][channel_name])
            cov_matrix  = gather_map_segments(cov_matrix_local_block_flat, 'cov_matrix', config['nside'][channel_name])
            hitmap  = gather_map_segments(hitmap_local_block, 'hitmap', config['nside'][channel_name])
            if rank == 0:
                io_obj.write_map(sky_map, 'sky_map', channel_name)
                io_obj.write_map(cov_matrix, 'cov_map', channel_name)
                io_obj.write_map(hitmap, 'hit_map', channel_name)
        else:
            cov_matrix = cov_ut.get_covariance_matrix(inv_cov_matrix_local, hitmap_local, config['pol_type'])
            sky_map = cov_ut.get_sky_map(cov_matrix_full, b_matrix_local, hitmap_local, config['pol_type'])
            cov_ut.flatten_block_matrix(cov_matrix_full, config['pol_type'])
            io_obj.write_map(sky_map, 'sky_map', channel_name)
            io_obj.write_map(cov_ut.flatten_block_matrix(cov_matrix, config['pol_type']), 'cov_map', channel_name)
            io_obj.write_map(hitmap_local, 'hit_map', channel_name)


#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
        
#Where all the distribution and gathering takes place
def distribute_matrix(local_full_matrix, matrix_type, nside):
    npix = hp.nside2npix(nside)
    dim, ind_elements = cov_ut.get_dim(config['pol_type'])
    segment_length = int(npix/size)

    if matrix_type == "hitmap":
        local_segmented_matrix = np.zeros(segment_length) 
    if matrix_type == "cov_matrix":
        local_segmented_matrix = np.zeros((segment_length, ind_elements)) 
    if matrix_type == "b_matrix":
        local_segmented_matrix = np.zeros((segment_length, dim)) 

    for i in range(size):
        start = i*segment_length
        stop = (i+1)*segment_length
        comm.Reduce(local_full_matrix[start:stop], local_segmented_matrix, MPI.SUM, root=i)
        
    return local_segmented_matrix

def gather_map_segments(map_segment, map_type, nside):
    npix = hp.nside2npix(nside)
    dim, ind_elements = cov_ut.get_dim(config['pol_type'])
    segment_length = npix/size 

    total_map = None

    if map_type == "sky_map":
        if rank == 0:
            total_map = np.empty(npix*dim, dtype=map_segment.dtype)
        if dim == 1:
            comm.Gather(map_segment, total_map, root=0)
        else:
            map_segment = map_segment.T.flatten()
            comm.Gather(map_segment, total_map, root=0)
            if rank == 0:
                total_map = total_map.reshape((npix,dim)).T

    if map_type == "hitmap":
        if rank == 0:
            total_map = np.empty(npix, dtype=map_segment.dtype)
        comm.Gather(map_segment, total_map, root=0)

    if map_type in  ["cov_matrix", "inv_cov_matrix"]:
        if rank == 0:
            total_map = np.empty(npix*ind_elements, dtype=map_segment.dtype)
        map_segment = map_segment.flatten()
        comm.Gather(map_segment, total_map, root=0)
        if rank == 0:
            total_map = total_map.reshape((npix,ind_elements)).T

    return total_map

#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*

def global_start_message():
    print_header()
    print_header(text='GENESYS')
    print_header(text='SIMULATING SATELLITE TOD')
    print_header()
    prompt(f"{'RUN TYPE':<35}{in_args.run_type}")
    prompt(f"{'# OF PROCESSES':<35}{size}")
    prompt(f"{'SIMULATION TAG':<35}{config['sim_tag']}")
    prompt(f"{'SPECIAL TAG':<35}{config['special_tag']}")
    prompt(f"{'NSIDE':<35}{config['nside']}")
    prompt(f"{'DETECTOR LIST':<35}{config['channel_detector_dict']}")
    prompt(f"{'NUMBER OF DATA BLOCKS':<35}{sd.num_data_blocks_global}")
    prompt(f"{'NUMBER OF SEGMENTS PER DATA BLOCK':<35}{config['num_segments_per_data_block']}")
    print_header()

def local_start_message():
    print_header(text=f'Rank {rank}', text_align='left')
    prompt(f"{'LOCAL DATA BLOCKS':<35}{sd.data_block_list_local}")
    prompt(f"{'LOCAL NUMBER OF DATA BLOCKS':<35}{sd.num_data_blocks_local}")
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
    config_file = os.path.join(current_dir,'config_files', in_args.config_file)
    config = load_param_file(file_path=config_file)
    # segment_list can be given as an explicit list, range(), np.arange or some other convenient way
    if type(config['data_block_list']) != list:
            config['data_block_list'] = list(eval(config['data_block_list']))
    # Segment distributor
    sd = Segment_Distributor(config['data_block_list'], config['num_segments_per_data_block'], size, rank)
    # I/O handling object
    io_obj = GenIO(config)
    if rank == 0:
        io_obj.make_top_data_directories(dir_list=['recon_dir'])
    #*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*

    #*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
    # The opening display messages
    #*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
    if rank == 0:
        global_start_message()
    if in_args.run_type == "mpi":
        comm.Barrier()
        time.sleep(rank*0.1)
        local_start_message()
        comm.Barrier()
    #*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*

    # Everything is set up now for the map-makingrun.
    prompt(f"*** Rank {rank:^6} BEGINNING MAP-MAKING ***", nature='info')
    run_map_maker()
    stop_time_main = time.time()
    prompt(f"*** RANK {rank:^6} COMPLETE ***          Time taken: {stop_time_main - start_time_main:.2f}s", nature='info')

    # Waiting for all the processes to finish
    if in_args.run_type == "mpi":
        comm.Barrier()
    if rank == 0:
        print_header()
        print_header(text='MAP-MAKING COMPLETE')
        print_header()
