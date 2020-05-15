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


#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*

#
#  prompt(f"Rank {rank} segment list: {segment_list_local}")
#  for segment in segment_list_local:
    #  segment_name = str(segment).zfill(7)
    #  prompt(f"Rank {rank} doing segment {segment}")
    #  outname = 'LFT_40GHz_' + segment_name + '.h5'
    #  out_f = h5py.File(os.path.join(out_dir, outname), 'w')
#
    #  prefix = '/common/'
    #  out_f.create_dataset(prefix + 'fsamp', data=fsamp)
    #  out_f.create_dataset(prefix + 'nside', data=nside)
    #  out_f.create_dataset(prefix + 'npsi', data=npsi)
    #  out_f.create_dataset(prefix + 'det', data=np.string_(', '.join(det_list)))
#
    #  polang = np.radians([0.0, 90.0, 45.0, 135.0, 0.0, 90.0])
    #  mbang = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
    #  out_f.create_dataset(prefix + 'polang', data=polang)
    #  out_f.create_dataset(prefix + 'mbang', data=mbang)
    #  out_f[prefix + 'polang'].attrs['legend'] = ', '.join(det_list)
    #  out_f[prefix + 'mbang'].attrs['legend'] = ', '.join(det_list)
#
#
    #  for chunk in chunk_list:
        #  prefix = '/' + str(chunk+1).zfill(6) + '/common/'
        #  start_time = (segment - 1)*segment_size + chunk*chunk_size
        #  out_f.create_dataset(prefix + 'time', data=start_time)
        #  out_f[prefix + 'time'].attrs['type'] = 'second'
#
        #  out_f.create_dataset(prefix + 'vsun', data=np.random.normal(size=3))
        #  out_f[prefix + 'vsun'].attrs['info'] = '[x,y,z]'
        #  out_f[prefix + 'vsun'].attrs['coords'] = 'galactic'
        #  pix_array = [[],[],[]]
        #  i_start = chunk*nsamp
        #  i_stop = (chunk+1)*nsamp
        #  for det in det_list:
            #  det_file = h5py.File(os.path.join(tod_dir, det, segment_name) + '.hdf5', 'r')
            #  # pixels
            #  pixels = hp.ang2pix(nside, det_file['theta'][i_start:i_stop], det_file['phi'][i_start:i_stop])
            #  delta = np.diff(pixels)
            #  delta = np.insert(delta, 0, pixels[0])
            #  pix_array[0].append(delta)
            #  # psi
            #  psi_bins = np.linspace(0, 2*np.pi, num=npsi)
            #  psi_index = np.digitize(det_file['psi'][i_start:i_stop], psi_bins)
            #  delta = np.diff(psi_index)
            #  delta = np.insert(delta, 0, psi_index[0])
            #  pix_array[1].append(delta)
            #  # flag
            #  flag = np.ones(nsamp)
            #  delta = np.diff(flag)
            #  delta = np.insert(delta, 0, flag[0])
            #  pix_array[2].append(delta)
#
        #  h = huffman.Huffman("", nside)
        #  h.GenerateCode(pix_array)
        #  huffarray = np.append(np.append(np.array(h.node_max), h.left_nodes), h.right_nodes)
        #  out_f.create_dataset(prefix + 'hufftree', data=huffarray)
        #  out_f.create_dataset(prefix + 'huffsymb', data=h.symbols)
#
        #  for det in det_list:
            #  prefix = '/' + str(chunk+1).zfill(6) + '/' + det + '/'
            #  # signal
            #  signal = det_file['signal'][i_start:i_stop]
            #  out_f.create_dataset(prefix + 'tod', data=det_file['signal'][i_start:i_stop])
            #  # flag
            #  flag = np.ones(nsamp)
            #  delta = np.diff(flag)
            #  delta = np.insert(delta, 0, flag[0])
            #  out_f.create_dataset(prefix + 'flag', data=np.void(bytes(h.byteCode(delta))))
            #  # pixels
            #  pixels = hp.ang2pix(nside, det_file['theta'][i_start:i_stop], det_file['phi'][i_start:i_stop])
            #  delta = np.diff(pixels)
            #  delta = np.insert(delta, 0, pixels[0])
            #  out_f.create_dataset(prefix + 'pix', data=np.void(bytes(h.byteCode(delta))))
            #  # psi
            #  psi_bins = np.linspace(0, 2*np.pi, num=npsi)
            #  psi_index = np.digitize(det_file['psi'][i_start:i_stop], psi_bins)
            #  delta = np.diff(psi_index)
            #  delta = np.insert(delta, 0, psi_index[0])
            #  out_f.create_dataset(prefix + 'psi', data=np.void(bytes(h.byteCode(delta))))
            #  # scalars
            #  gain = 1.0
            #  sigma0 = 10.0
            #  fknee = 0.02
            #  alpha = 1.0
            #  out_f.create_dataset(prefix + 'scalars', data=np.array([gain, sigma0, fknee, alpha]).flatten())
            #  out_f[prefix + 'scalars'].attrs['legend'] = 'gain, sigma0, fknee, alpha'
#
    #  out_f.close()
