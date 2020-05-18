import argparse
import os
import numpy as np
import healpy as hp
from mpi4py import MPI
from genesys import load_param_file, prompt
from genesys.gen_io import GenIO
from genesys.gen_io.segment_distribution import Segment_Distributor
import genesys.gen_io.huffman as huffman

def huffman_compress():
    for data_block in sd.data_block_list_local:
        prompt(f"Rank {rank:^6} starting data block {data_block}", nature='info')
        for channel_name in config['channel_detector_dict'].keys():
            io_obj.open_tod_file(channel_name, data_block, 'r')
            io_obj_huff.open_tod_file(channel_name, data_block, 'w')
            segment_list = sd.get_segment_list(data_block, config['num_segments_per_data_block'])
            channel_commons = io_obj.read_channel_common(['fsamp', 'det', 'polang', 'mbangs', 'coords'])
            coords = channel_commons.pop('coords')
            channel_commons.update({'nside': config['huffman']['nside'], 'npsi': config['huffman']['npsi']})
            det_list = ','.join(config['channel_detector_dict'][channel_name])
            channel_common_attributes = {'polang': {'legend': det_list}, 'mbangs': {'legend': det_list}}
            io_obj_huff.write_channel_common(channel_commons, channel_common_attributes)
            for segment in segment_list: 
                segment_commons = {}
                segment_commons['time'] = io_obj.read_segment_common(segment, ['time_0'])['time_0']
                segment_commons.update({'vsun': [0.0,0.0,0.0]})
                segment_commons_attributes = {'vsun': {'info': '[x,y,z]', 'coords': coords}, 'time': {'type': 'second'}}
                pix_array = [[],[],[]]
                for detector_name in config['channel_detector_dict'][channel_name]:
                    tod = io_obj.read_detector_data(segment, detector_name, ['theta', 'phi', 'psi'])
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
                segment_commons.update({'hufftree': huffarray, 'huffsymb': h.symbols})
                io_obj_huff.write_segment_common(segment, segment_commons, segment_commons_attributes)
                for detector_name in config['channel_detector_dict'][channel_name]:
                    tod = io_obj.read_detector_data(segment, detector_name, ['noise', 'signal', 'theta', 'phi', 'psi'])
                    scalars = io_obj.read_detector_scalars(segment, detector_name, ['sigma0', 'gain', 'fknee', 'alpha'])
                    # signal, noise and scalars
                    io_obj_huff.write_detector_data(segment, detector_name, tod, scalars, tod_write_field=['signal'])
                    # flag
                    flag = np.ones(tod['signal'].size)
                    delta = np.diff(flag)
                    delta = np.insert(delta, 0, flag[0])
                    io_obj_huff.write_detector_data(segment, detector_name, {'flag': np.void(bytes(h.byteCode(delta)))}, tod_write_field=['flag'])
                    # pixels
                    pixels = hp.ang2pix(config['huffman']['nside'], tod['theta'], tod['phi'])
                    delta = np.diff(pixels)
                    delta = np.insert(delta, 0, pixels[0])
                    io_obj_huff.write_detector_data(segment, detector_name, {'pix': np.void(bytes(h.byteCode(delta)))}, tod_write_field=['pix'])
                    # psi
                    psi_bins = np.linspace(0, 2*np.pi, num=config['huffman']['npsi'])
                    psi_index = np.digitize(tod['psi'], psi_bins)
                    delta = np.diff(psi_index)
                    delta = np.insert(delta, 0, psi_index[0])
                    io_obj_huff.write_detector_data(segment, detector_name, {'psi': np.void(bytes(h.byteCode(delta)))}, tod_write_field=['psi'])
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
