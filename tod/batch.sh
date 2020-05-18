#!/usr/bin/bash

CONFIG_FILE=test.yaml
VERBOSITY=2
NUM_PROC=48
OUT_FILE=out.txt

# mpirun -n $NUM_PROC python sim_tod.py $CONFIG_FILE mpi -v $VERBOSITY
mpirun -n $NUM_PROC python huff_compress.py $CONFIG_FILE
# python sim_tod.py $CONFIG_FILE serial -v $VERBOSITY
# mpirun -n $NUM_PROC python sim_tod.py $CONFIG_FILE mpi -v $VERBOSITY > $OUT_FILE 2>&1
# mprof run python sim_tod.py $CONFIG_FILE serial -v $VERBOSITY
