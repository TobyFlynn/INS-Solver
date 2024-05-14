#!/bin/bash

for n in $(seq 0 3)
do
    out_dir=m$n
    mkdir $out_dir
    scale=$(echo 2^-$n | bc -l)
    gmsh \
        mesh.geo \
        -3 \
        -format vtk \
        -algo hxt \
        -nt 3 \
        -clscale $scale\
        -o $out_dir/mesh.vtk
    mpirun -n 1 vtk2hdf5_3D \
        -file $out_dir/mesh.vtk
    mv mesh.h5 $out_dir
done
