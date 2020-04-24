#!/usr/bin/bash

CONFIG_FILE=test.yaml
VERBOSITY=2
NUM_PROC=1
OUT_FILE=out.txt

mpirun -n $NUM_PROC python sim_tod.py $CONFIG_FILE mpi -v $VERBOSITY
# mpirun -n $NUM_PROC python sim_timestream.py $CONFIG_FILE mpi -v $VERBOSITY > $OUT_FILE 2>&1
