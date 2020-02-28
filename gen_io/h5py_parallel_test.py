from mpi4py import MPI
import h5py

comm = MPI.COMM_WORLD

rank = comm.rank

f = h5py.File('parallel_test.hdf5', 'w', driver='mpio', comm=comm)

dset = f.create_dataset('test', (4,), dtype='i')
dset[rank] = rank

f.close()
