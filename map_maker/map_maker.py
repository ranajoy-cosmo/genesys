#!/usr/bin/env python

import numpy as np
import healpy as hp
import sys
import time
import importlib
from termcolor import colored
import genesys.data_io.segment_distribution as segd
import genesys.data_io.data_io as dio
from genesys.utilities import Generic_Class, prompt
import genesys.utilities.unit_conversion as uc
from genesys.timestream_simulation.detector import Detector
import genesys.map_maker.covariance_matrix_utils as cov_ut

def run_map_maker():
    #Defining the dimensions of the arrays and matrices
    npix = hp.nside2npix(config.nside_out)
    dim, ind_elements = cov_ut.get_dim(config.recon_pol_type)

    # The matrices local to the process
    inv_cov_matrix_local = np.zeros((npix, ind_elements), dtype=np.float)
    b_matrix_local = np.zeros((npix, dim), dtype=np.float)
    hitmap_local = np.zeros(npix, dtype=np.float)

    count = 0
    total_time_taken = 0
    # Iterating over the detectors
    for band_name in list(local_band_detector_segment_dict.keys()):
        for detector_name in list(local_band_detector_segment_dict[band_name].keys()):
            if config.pair_difference:
                detector_a = Detector(band_name, dio.get_detector_a_name(detector_name), config)
                detector_b = Detector(band_name, dio.get_detector_b_name(detector_name), config)
            else:
                detector = Detector(band_name, detector_name, config)
            for segment in local_band_detector_segment_dict[band_name][detector_name]:
                segment_start_time = time.time()
                count += 1
                if config.pair_difference:
                    detector_a.segment_start_prompt(rank, segment)
                    detector_b.segment_start_prompt(rank, segment)
                    signal, hitpix, pol_ang = acquire_difference_signal(detector_a, detector_b, segment)
                else:
                    detector.segment_start_prompt(rank, segment)
                    signal, hitpix, pol_ang = acquire_signal(detector, segment)
                if config.subtract_template:
                    for template_name in config.template_list:
                        if config.pair_difference:
                            template_tod = acquire_template_signal(detector_a, segment, template_name)
                            template_amplitude = get_template_amplitude(detector_a, segment, template_name)
                        else:
                            template_tod = acquire_template_signal(detector, segment, template_name)
                            template_amplitude = get_template_amplitude(detector, segment, template_name)
                        signal -= template_amplitude * template_tod
                    del template_tod
                cov_ut.add_to_mm_matrices(hitpix, pol_ang, signal, inv_cov_matrix_local, b_matrix_local, hitmap_local, npix, config.recon_pol_type)
                segment_stop_time = time.time()
                time_taken = segment_stop_time - segment_start_time
                total_time_taken += segment_stop_time - segment_start_time
                prompt("Rank : {}, Time taken : {}. Total time take : {}, Projected time : {}, Finished {} of {}\n".format(rank, time_taken, total_time_taken, total_time_taken*local_num_segments/count, count, local_num_segments))

    del signal
    del pol_ang
    del hitpix

    if run_type == "mpi":
        # Distributing and gathering the segments of the matrices in the proper process
        # The sky pixels are chunked and are handled by individual processes
        inv_cov_matrix_local_segment = distribute_matrix(inv_cov_matrix_local, "cov_matrix")
        del inv_cov_matrix_local
        b_matrix_local_segment = distribute_matrix(b_matrix_local, "b_matrix")
        del b_matrix_local
        hitmap_local_segment = distribute_matrix(hitmap_local, "hitmap")
        del hitmap_local

        # Inverting the local segment of the inverse covariance matrix
        cov_matrix_local_segment = cov_ut.get_covariance_matrix(inv_cov_matrix_local_segment, hitmap_local_segment, config.recon_pol_type)

        # estimating the local sky segment
        sky_map_local_segment = cov_ut.get_sky_map(cov_matrix_local_segment, b_matrix_local_segment, hitmap_local_segment, config.recon_pol_type)
        cov_matrix_local_segment = cov_ut.flatten_block_matrix(cov_matrix_local_segment, config.recon_pol_type)

        sky_map_total = gather_map_segments(sky_map_local_segment, 'sky_map')
        if rank == 0:
            dio.write_sky_map(config, sky_map_total, "sky_map")
            del sky_map_total
        hitmap_total = gather_map_segments(hitmap_local_segment, 'hitmap')
        if rank == 0:
            dio.write_sky_map(config, hitmap_total, "hitmap")
            del hitmap_total
        cov_matrix_total = gather_map_segments(cov_matrix_local_segment, 'cov_matrix')
        if rank == 0:
            dio.write_sky_map(config, cov_matrix_total, "covariance_matrix")
            del cov_matrix_total
        inv_cov_matrix_total = gather_map_segments(inv_cov_matrix_local_segment, 'inv_cov_matrix')
        if rank == 0:
            dio.write_sky_map(config, inv_cov_matrix_total, "inverse_covariance_matrix")

    else:
        cov_matrix = cov_ut.get_covariance_matrix(inv_cov_matrix_local, hitmap_local, config.recon_pol_type)
        sky_map = cov_ut.get_sky_map(cov_matrix, b_matrix_local, hitmap_local, config.recon_pol_type)
        dio.write_sky_map(config, sky_map, "sky_map")
        dio.write_sky_map(config, hitmap_local, "hitmap")
        dio.write_sky_map(config, cov_ut.flatten_block_matrix(cov_matrix, config.recon_pol_type).T, "covariance_matrix")
        dio.write_sky_map(config, inv_cov_matrix_local.T, "inverse_covariance_matrix")

#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*

def acquire_signal(detector, segment):
    if config.simulate_ts:
        time_stream_dict = detector.simulate_timestream_data(segment, return_field=["signal", "pointing_vec", "pol_ang"])
    else:
        time_stream_dict = detector.read_timestream_data(segment, return_field=["signal", "pointing_vec", "pol_ang"], noise_only=config.noise_only)
    hitpix = cov_ut.pointing_vec_to_hitpix(time_stream_dict["pointing_vec"], config.nside_out)
    return time_stream_dict["signal"], hitpix, time_stream_dict["pol_ang"]

def acquire_difference_signal(det_a, det_b, segment):
    if config.simulate_ts:
        t_stream_a = det_a.simulate_timestream_signal(segment, return_field=["signal", "pointing_vec", "pol_ang"])
        t_stream_b = det_b.simulate_timestream_signal(segment, return_field=["signal"])
    else:
        t_stream_a = det_a.read_timestream(segment, return_field=["signal", "pointing_vec", "pol_ang"], noise_only=config.noise_only)
        t_stream_b = det_b.read_timestream(segment, return_field=["signal"], noise_only=config.noise_only)

    signal = 0.5*(t_stream_a["signal"] - t_stream_b["signal"])

    hitpix = cov_ut.pointing_vec_to_hitpix(time_stream_dict["pointing_vec"], config.nside_out)
    return signal, hitpix, t_stream_a["pol_ang"]

def acquire_signal_template(det, segment, template_name):
    det.config.special_tag = template_name
    t_stream = det.read_timestream(segment, return_field=["signal"])

    return t_stream[signal_type], hitpix, t_stream["pol_ang"]

#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
        
#Where all the distribution and gathering takes place
def distribute_matrix(local_full_matrix, matrix_type):
    npix = hp.nside2npix(config.nside_out)
    dim, ind_elements = cov_ut.get_dim(config.recon_pol_type)
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

def gather_map_segments(map_segment, map_type):
    npix = hp.nside2npix(config.nside_out)
    dim, ind_elements = cov_ut.get_dim(config.recon_pol_type)
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
    display_string = colored("\n#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*\n", "green")
    display_string += colored("#* BEGINNING MAP-MAKER\n", "green")
    display_string += colored("#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*\n", "green")
    display_string += "RUN TYPE : {}\n".format(run_type)
    display_string += "# OF PROCESSES: {}\n".format(size)
    display_string += "SIMULATION TAG: {}\n".format(config.sim_tag)
    display_string += "SCAN TAG: {}\n".format(config.scan_tag)
    display_string += "RECON TAG: {}\n".format(config.recon_tag)
    display_string += "SPECIAL TAG: {}\n".format(config.special_tag)
    display_string += "NSIDE OUT: {}\n".format(config.nside_out)
    if config.simulate_ts:
        display_string += "SIMULATION POLARISATION TYPE: {}\n".format(config.sim_pol_type)
        display_string += "SCAN STRATEGY : {}\n".format(config.scan_strategy_name)
        display_string += "COORDINATE SYSTEM : {}\n".format(config.coordinate_system)
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
    display_string += "# OF DETECTORS: {}\n".format(len(detector_list))
    display_string += "DETECTOR LIST: {}\n".format(', '.join(detector_list))
    display_string += "TOTAL NUMBER OF SEGMENTS : {}\n".format(num_segments)
    display_string += "TOTAL SIMULATION TIME : {}s {}d {}y\n".format(total_simulation_time, uc.convert_unit('time', total_simulation_time, 'second', 'day'), uc.convert_unit('time', total_simulation_time, 'second', 'siderial year'))
    display_string += "NOTES : {}\n".format(config.notes)
    display_string += colored("#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*\n\n", "green")
    prompt(display_string, sys.stdout)

def local_start_message():
    display_string = colored("#* RANK: {} #*#*#*#*#*#*#*#*#*#*\n".format(rank), color="yellow")
    display_string += "LOCAL # OF DETECTORS: {}\n".format(len(local_detector_list))
    display_string += "LOCAL DETECTOR LIST: {}\n".format(', '.join(local_detector_list))
    display_string += "LOCAL NUMBER OF SEGMENTS: {}\n".format(local_num_segments)
    display_string += "TOTAL LOCAL SIMULATION TIME : {}s {}d {}y\n".format(local_simulation_time, uc.convert_unit('time', local_simulation_time, 'second', 'day'), uc.convert_unit('time', local_simulation_time, 'second', 'siderial year'))
    display_string += colored("#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*\n\n", color="yellow")
    prompt(display_string)

#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
#* Main function definition. This is where the code begins when executed
#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*

if __name__=="__main__":
    start_time_main = time.time()

    # run_type can be either "mpi" or "serial"
    run_type = sys.argv[1]
    # This block just sets up the rank amd size of the simulation run
    if run_type == "mpi":
        from mpi4py import MPI          # This module is imported here to be compatible with systems not having mpi support
        comm = MPI.COMM_WORLD
        rank = comm.Get_rank()          # Rank of this particular MPI process. Ranges from 0 to (size - 1)
        size = comm.Get_size()          # Number of MPI processes
    elif run_type == "serial":
        rank = 0
        size = 1
    else:
        raise Exception("Please enter \"mpi\" OR \"serial\" as the second argument.")

    # config_file is provided as a dot separated path to the config file with genesys as the parent directory
    config_file = sys.argv[2]
    # Importing the general simulation parameters
    config = importlib.import_module(config_file).config
    # Importing the scan strategy
    scan_strategy = importlib.import_module("genesys.scan_strategy." + config.scan_strategy_config).scan_strategy
    config.__dict__.update(scan_strategy.__dict__)

    npix = 12 * config.nside_out**2
    assert (npix % size == 0), "# processes is {}. Please make # processes a factor of npix={}".format(size, npix)

    # Getting the global list of detectors and segments
    num_segments = segd.count_segments(config.band_detector_segment_dict)
    detector_list = segd.get_detector_list(config.band_detector_segment_dict)
    total_simulation_time = segd.get_total_simulation_time(config.band_detector_segment_dict)
    # Getting the local list of detectors and segments
    local_band_detector_segment_dict = segd.get_local_band_detector_segment_dict(rank, size, config.band_detector_segment_dict)
    local_num_segments = segd.count_segments(local_band_detector_segment_dict)
    local_detector_list = segd.get_detector_list(local_band_detector_segment_dict)
    local_simulation_time = segd.get_total_simulation_time(local_band_detector_segment_dict)

    # Making the data directories and prompting the global and local start messages
    if rank == 0: 
        dio.make_top_data_directories(config, dir_list=['sim_dir', 'scan_dir', 'recon_dir'])
        global_start_message()
    if run_type == "mpi":
        comm.Barrier()
    local_start_message()

    # Displaying the parameters ofthe individual detectors
    if rank == 0:
        for band_name in list(local_band_detector_segment_dict.keys()):
            for detector_name in list(local_band_detector_segment_dict[band_name].keys()):
                detector = Detector(band_name, detector_name, config, load_map=True)
                detector.display_params()
    # Waiting for all the processes to print out their respective detector parameters
    if run_type == "mpi":
        comm.Barrier()

    run_map_maker()

    if run_type == "mpi":
        comm.Barrier()
    if rank == 0:
        prompt(colored("########## MAP-MAKING DONE ##########\n", "green"))

#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
#  def write_segments(maps, map_name, recon_dir):
    #  np.save(os.path.join(recon_dir, map_name + "_segments", str(rank).zfill(4)), maps)
#
#
#  #The different sky chunks are brought together to form the final maps
#  def accumulate_segments(size):
    #  acc_start_time = time.time()
    #  dir_names = get_dir_names()
#
    #  npix = hp.nside2npix(config.nside_out)
    #  npix_segment = npix/size
    #  dim, ind_elements = cov_ut.get_dim(config.pol_type)
#
    #  # SKY MAP
    #  sky_map = np.empty((dim, npix))
    #  for i in range(size):
        #  start = i*npix_segment
        #  stop = (i+1)*npix_segment
        #  sky_map_segment = np.load(os.path.join(dir_names.recon_dir, dir_names.map_segment_dir, str(i).zfill(4) + '.npy'))
        #  sky_map[..., start:stop] = sky_map_segment
#
    #  #if config.pol_type == 'T':
    #  #    sky_map[np.isnan(sky_map)] = 0.0
    #  #else:
    #  #    sky_map[..., np.isnan(sky_map[0])] = 0.0
#
    #  hp.write_map(os.path.join(dir_names.recon_dir, "sky_map.fits"), sky_map)
    #  del sky_map
    #  del sky_map_segment
    #  shutil.rmtree(os.path.join(dir_names.recon_dir, dir_names.map_segment_dir))
    #  #*#*#*#*#*#*#*#*#*#*
#
    #  # HITMAP
    #  hitmap = np.empty(npix)
    #  for i in range(size):
        #  start = i*npix_segment
        #  stop = (i+1)*npix_segment
        #  hitmap_segment = np.load(os.path.join(dir_names.recon_dir, dir_names.hitmap_segment_dir, str(i).zfill(4) + '.npy'))
        #  hitmap[..., start:stop] = hitmap_segment
#
    #  hp.write_map(os.path.join(dir_names.recon_dir, "hitmap.fits"), hitmap)
    #  if config.pol_type == "IQU":
        #  mask_map = hitmap > 3
    #  elif config.pol_type == "QU":
        #  mask_map = hitmap > 2
    #  else:
        #  mask_map = hitmap > 1
    #  hp.write_map(os.path.join(dir_names.recon_dir, "mask.fits"), mask_map)
#
    #  del hitmap
    #  del hitmap_segment
    #  shutil.rmtree(os.path.join(dir_names.recon_dir, dir_names.hitmap_segment_dir))
    #  #*#*#*#*#*#*#*#*#*#*
#
    #  # INVERSE COVARIANCE MATRIX
    #  inverse_cov_matrix = np.empty((npix, ind_elements))
    #  for i in range(size):
        #  start = i*npix_segment
        #  stop = (i+1)*npix_segment
        #  inverse_cov_matrix_segment = np.load(os.path.join(dir_names.recon_dir, dir_names.inv_cov_matrix_segment_dir, str(i).zfill(4) + '.npy'))
        #  inverse_cov_matrix[start:stop] = inverse_cov_matrix_segment
#
    #  hp.write_map(os.path.join(dir_names.recon_dir, "inverse_covariance_maps.fits"), inverse_cov_matrix.T)
    #  del inverse_cov_matrix
    #  del inverse_cov_matrix_segment
    #  shutil.rmtree(os.path.join(dir_names.recon_dir, dir_names.inv_cov_matrix_segment_dir))
    #  #*#*#*#*#*#*#*#*#*#*
#
    #  # COVARIANCE MATRIX
    #  cov_matrix = np.empty((npix, ind_elements))
    #  for i in range(size):
        #  start = i*npix_segment
        #  stop = (i+1)*npix_segment
        #  cov_matrix_segment = np.load(os.path.join(dir_names.recon_dir, dir_names.cov_matrix_segment_dir, str(i).zfill(4) + '.npy'))
        #  if config.pol_type == "IQU":
            #  cov_matrix_segment = cov_matrix_segment.reshape((npix_segment, dim**2))[..., np.array([0,1,2,4,5,8])]
        #  elif config.pol_type =="QU":
            #  cov_matrix_segment = cov_matrix_segment.reshape((npix_segment, dim**2))[..., np.array([0,1,3])]
        #  else:
            #  cov_matrix_segment = cov_matrix_segment.reshape((npix_segment, dim**2))
        #  cov_matrix[start:stop] = cov_matrix_segment
#
    #  hp.write_map(os.path.join(dir_names.recon_dir, "covariance_maps.fits"), cov_matrix.T)
    #  del cov_matrix
    #  del cov_matrix_segment
    #  shutil.rmtree(os.path.join(dir_names.recon_dir, dir_names.cov_matrix_segment_dir))
    #  #*#*#*#*#*#*#*#*#*#*
#
    #  acc_stop_time = time.time()
    #  prompt("Time taken to accumulate segments : {}s\n".format(acc_stop_time - acc_start_time))
#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
