#!/usr/bin/bash

CONFIG_FILE=cmb_only_white_noise.yaml
VERBOSITY=2
NUM_PROC=122
OUT_FILE=out.txt

mpirun -n $NUM_PROC python sim_tod.py $CONFIG_FILE mpi -v $VERBOSITY  > $OUT_FILE 2>&1 &&
mpirun -n $NUM_PROC python huff_compress.py $CONFIG_FILE &&
python move_params_and_output.py $CONFIG_FILE $OUT_FILE

# python sim_tod.py $CONFIG_FILE serial -v $VERBOSITY
#
# CONFIG_FILE=cmb_with_doppler_white_noise.yaml
#
# mpirun -n $NUM_PROC python sim_tod.py $CONFIG_FILE mpi -v $VERBOSITY > $OUT_FILE 2>&1 &&
# mpirun -n $NUM_PROC python huff_compress.py $CONFIG_FILE &&
# python move_params_and_output.py $CONFIG_FILE $OUT_FILE
