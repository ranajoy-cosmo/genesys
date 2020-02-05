import h5py
import numpy as np
import healpy as hp
import huffman
import os

vel_file = "/mn/stornext/u3/ranajoyb/genesys/data/satellite_velocity.fits"
tod_dir = "/mn/stornext/u3/ranajoyb/genesys/output/sim_test_old/scan_test/LFT_40GHz"
det_list = ["0001a", "0001b", "0002a", "0002b", "0003a", "0003b"]
segment = "0000002"
out_dir = "/mn/stornext/u3/ranajoyb/genesys/output/sim_test_old/"
outname = 'LFT_40GHz_' + segment + '.h5'

out_f = h5py.File(os.path.join(out_dir, outname), 'w')

nside = 256
npsi = 4096
fsamp = 31
nsamp = 4076190

out_f.create_dataset('/common/fsamp', data=fsamp)
out_f.create_dataset('/common/nside', data=nside)
out_f.create_dataset('/common/npsi', data=npsi)
out_f.create_dataset('/common/det', data=', '.join(det_list))

polang = np.radians([0.0, 90.0, 45.0, 135.0, 0.0, 90.0])
mbang = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
out_f.create_dataset('/common/polang', data=polang)
out_f.create_dataset('/common/mbang', data=mbang)
out_f['/common/polang'].attrs['legend'] = ', '.join(det_list)
out_f['/common/mbang'].attrs['legend'] = ', '.join(det_list)

out_f.create_dataset('common/time', data =np.arange(nsamp) / fsamp)
out_f['/common/time'].attrs['type'] = 'second'

out_f.create_dataset('common/vsun', data=np.random.normal(size=(3,nsamp)))
out_f['/common/vsun'].attrs['info'] = '[x,y,z]'
out_f['/common/vsun'].attrs['coords'] = 'galactic'

pix_array = [[],[],[]]
for det in det_list:
    det_file = h5py.File(os.path.join(tod_dir, det, segment) + '.hdf5', 'r')
    # pixels
    pixels = hp.ang2pix(nside, det_file['theta'][:], det_file['phi'][:])
    delta = np.diff(pixels)
    delta = np.insert(delta, 0, pixels[0])
    pix_array[0].append(delta)
    # psi
    psi_bins = np.linspace(0, 2*np.pi, num=npsi)
    psi_index = np.digitize(det_file['psi'][:], psi_bins)
    delta = np.diff(psi_index)
    delta = np.insert(delta, 0, psi_index[0])
    pix_array[1].append(delta)
    # flag
    flag = np.ones(nsamp)
    delta = np.diff(flag)
    delta = np.insert(delta, 0, flag[0])
    pix_array[2].append(delta)

h = huffman.Huffman("", nside)
h.GenerateCode(pix_array)
huffarray = np.append(np.append(np.array(h.node_max), h.left_nodes), h.right_nodes)
out_f.create_dataset('/common/hufftree', data=huffarray)
out_f.create_dataset('/common/huffsymb', data=h.symbols)

for det in det_list:
    det_file = h5py.File(os.path.join(tod_dir, det, segment) + '.hdf5', 'r')
    # signal
    out_f.create_dataset(det+'/tod', data=det_file['signal'])
    # flag
    flag = np.ones(nsamp)
    delta = np.diff(flag)
    delta = np.insert(delta, 0, flag[0])
    out_f.create_dataset(det+'/flag', data=np.void(bytes(h.byteCode(delta))))
    # pixels
    pixels = hp.ang2pix(nside, det_file['theta'][:], det_file['phi'][:])
    delta = np.diff(pixels)
    delta = np.insert(delta, 0, pixels[0])
    out_f.create_dataset(det+'/pix', data=np.void(bytes(h.byteCode(delta))))
    # psi
    psi_bins = np.linspace(0, 2*np.pi, num=npsi)
    psi_index = np.digitize(det_file['psi'][:], psi_bins)
    delta = np.diff(psi_index)
    delta = np.insert(delta, 0, psi_index[0])
    out_f.create_dataset(det+'/psi', data=np.void(bytes(h.byteCode(delta))))

out_f.close()
