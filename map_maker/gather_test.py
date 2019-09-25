#!/usr/bin/env python

import numpy as np
import healpy as hp
import matplotlib.pyplot as plt
from mpi4py import MPI
comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()

nside = 128
npix = 12*128**2
dim = 1

sky_map = hp.read_map("/mn/stornext/u3/ranajoyb/genesys/maps/map_files/test_map.fits", field=(0))
sky_map = hp.ud_grade(sky_map, nside_out=128)

local_sky_map = sky_map[..., int(rank * npix/size): int((rank+1) * npix/size)].T.flatten()

print(local_sky_map.shape)

total_sky_map = None
if rank == 0:
    total_sky_map = np.empty(dim*npix, dtype=np.float)

comm.Gather(local_sky_map, total_sky_map, root=0)

if rank == 0:
    total_sky_map = total_sky_map.reshape((npix, dim)).T
    hp.mollview(total_sky_map)
    plt.show()
