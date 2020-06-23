import argparse
import os
import numpy as np
import healpy as hp
import sys
from mpi4py import MPI
from genesys import load_param_file, prompt
from genesys.gen_io import GenIO
from genesys.gen_io.segment_distribution import Segment_Distributor
import genesys.gen_io.huffman as huffman

"""
The structure of the Huffman encoded HDF5 file
<channel_name>_<segment_name>.h5
    - common : (GROUP) Contains common data values for the segment and all detectors 
        - fsamp : (DSET) Sampling frequency
        - nside : (DSET) NSide of HEALPix map used for Huffman compressing pointing
        - npsi : (DSET) Discretisation over 2pi for Huffman compressing psi
        - det : (DSET) List of detector names in np.string_ format
        - polang : (DSET) List of relative orientation of polarisation axis.
            - legend : (ATTR) List of detectors
        - mbang : (DSET) List of relative orientation of main beam
            - legend : (ATTR) List of detectors
    - segment : (GROUP) The TODs from the different detector are here
        - common : (GROUP) The common TODs for the entire detector set
            - time : (DSET) The time steps.
                - unit : (ATTR) Time unit and convention
            - vsun : (DSET) Solar orbital velocity.
                - info : (ATTR) axes
                - coords : (ATTR) Coordinate system
            - psi_hwp : (DSET) The angle of the HWP at the first time sample
            - hufftree : (DSET)
            - huffsymb : (DSET)
        - <detector_name> : (GROUP) The unique TOD for the particular detector
            - tod : (DSET) The signal observed by the detector
            - flag : (DSET) Flagging of invalid data by 0. (*Huffman)
            - pixels : (DSET) Discretised pointing at a given HEALPix NSide. (*Huffman)
            - psi : (DSET) Discretised polarisation angle values. (*Huffman)
            - scalars : Array of gain and noise properties. [gain, sigma0, fknee, alpha]
                - legend : Description of the scalar elements as a ',' separated string.
"""

def huffman_compress():
    for data_block in sd.data_block_list_local:
        prompt(f"Rank {rank:^6} starting data block {data_block}", nature='info')
        for channel_name in config['channel_detector_dict'].keys():
            io_obj.open_tod_file(channel_name, data_block, 'r')
            io_obj_huff.open_tod_file(channel_name, data_block, 'w')
            segment_list = sd.get_segment_list(data_block, config['num_segments_per_data_block'])
            channel_common, channel_common_attributes = io_obj.read_channel_common(['fsamp', 'det', 'polang', 'mbang', 'coords'])
            coords = channel_common.pop('coords')
            channel_common.update({'nside': config['huffman']['nside'], 'npsi': config['huffman']['npsi']})
            io_obj_huff.write_channel_common(channel_common, channel_common_attributes)
            for segment in segment_list: 
                segment_common, segment_common_attributes = io_obj.read_segment_common(segment, ['time', 'vsun', 'psi_hwp'])
                pix_array = [[],[],[]]
                for detector_name in config['channel_detector_dict'][channel_name]:
                    tod = io_obj.read_tod(segment, detector_name, ['theta', 'phi', 'psi'])
                    # pixels
                    pixels = hp.ang2pix(config['huffman']['nside'], tod['theta'], tod['phi'])
                    delta = np.diff(pixels)
                    delta = np.insert(delta, 0, pixels[0])
                    pix_array[0].append(delta)
                    # psi
                    psi_bins = np.linspace(0, 2*np.pi, num=config['huffman']['npsi'])
                    psi_index = np.digitize(tod['psi'], psi_bins)
                    delta = np.diff(psi_index)
                    delta = np.insert(delta, 0, psi_index[0])
                    pix_array[1].append(delta)
                    # flag
                    flag = np.ones(tod['theta'].size)
                    delta = np.diff(flag)
                    delta = np.insert(delta, 0, flag[0])
                    pix_array[2].append(delta)
                h = huffman.Huffman("", config['huffman']['nside'])
                h.GenerateCode(pix_array)
                huffarray = np.append(np.append(np.array(h.node_max), h.left_nodes), h.right_nodes)
                segment_common.update({'hufftree': huffarray, 'huffsymb': h.symbols})
                io_obj_huff.write_segment_common(segment, segment_common, segment_common_attributes)
                for detector_name in config['channel_detector_dict'][channel_name]:
                    tod = io_obj.read_tod(segment, detector_name, ['signal', 'theta', 'phi', 'psi'])
                    tod['tod'] = tod.pop('signal')
                    scalars = io_obj.read_tod_scalars(segment, detector_name)
                    # signal, noise and scalars
                    io_obj_huff.write_tod(segment, detector_name, tod, tod_write_field=['tod'])
                    io_obj_huff.write_tod_scalars(segment, detector_name, scalars)
                    # flag
                    flag = np.ones(tod['tod'].size)
                    delta = np.diff(flag)
                    delta = np.insert(delta, 0, flag[0])
                    io_obj_huff.write_tod(segment, detector_name, {'flag': np.void(bytes(h.byteCode(delta)))}, tod_write_field=['flag'])
                    # pixels
                    pixels = hp.ang2pix(config['huffman']['nside'], tod['theta'], tod['phi'])
                    delta = np.diff(pixels)
                    delta = np.insert(delta, 0, pixels[0])
                    io_obj_huff.write_tod(segment, detector_name, {'pix': np.void(bytes(h.byteCode(delta)))}, tod_write_field=['pix'])
                    # psi
                    psi_bins = np.linspace(0, 2*np.pi, num=config['huffman']['npsi'])
                    psi_index = np.digitize(tod['psi'], psi_bins)
                    delta = np.diff(psi_index)
                    delta = np.insert(delta, 0, psi_index[0])
                    io_obj_huff.write_tod(segment, detector_name, {'psi': np.void(bytes(h.byteCode(delta)))}, tod_write_field=['psi'])
            io_obj.close_tod_file()
            io_obj_huff.close_tod_file()

if __name__=='__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('config_file', type=str, help='Name of config file under relative path config_files/')
    in_args = parser.parse_args()

    comm = MPI.COMM_WORLD
    size = comm.Get_size()
    rank = comm.Get_rank()

    current_dir = os.path.dirname(__file__)
    config_file = os.path.join(current_dir,'config_files', in_args.config_file)
    config = load_param_file(file_path=config_file)
    if type(config['data_block_list']) != list:
            config['data_block_list'] = list(eval(config['data_block_list']))

    sd = Segment_Distributor(config['data_block_list'], config['num_segments_per_data_block'], size, rank)
    io_obj = GenIO(config)
    io_obj_huff = GenIO(config, huffman=True)
    io_obj_huff.make_top_data_directories(['tod_dir'])

    huffman_compress()
