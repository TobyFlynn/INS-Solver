#!/bin/bash

for i in $(seq 0 3)
do
    mpirun -n 1 ../../../bin/ins2d_mpi \
        -config config.ini \
        -input m${i}/mesh.h5 \
        -output m${i}/
done
