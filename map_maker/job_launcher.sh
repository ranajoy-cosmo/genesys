#!/bin/bash

NUM_PROC=96
export OMP_NUM_THREADS=1
CONFIG_FILE=genesys.map_maker.config_files.cmb_bharat

mpiexec -n $NUM_PROC python map_maker.py mpi $CONFIG_FILE
# python map_maker.py serial $CONFIG_FILE
