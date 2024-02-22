#!/bin/bash

for n in $(seq 0 3)
do
    out_dir=m$n
    mkdir $out_dir
    scale=$(echo 2^-$n | bc -l)
    gmsh \
        mesh.geo \
        -2 \
        -format vtk \
        -algo front2d \
        -nt 3 \
        -clscale $scale\
        -o $out_dir/mesh.vtk
    mpirun -n 1 ${OP2_DG_TOOLKIT_DIR}/bin/vtk2hdf5_2D \
        -file $out_dir/mesh.vtk \
        -bc none
    mv mesh.h5 $out_dir
done
