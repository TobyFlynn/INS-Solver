#!/bin/bash

gmsh \
    mesh.geo \
    -2 \
    -format vtk \
    -algo del2d \
    -nt 2 

mpirun -n 1 ${OP2_DG_TOOLKIT_DIR}/bin/vtk2hdf5_2D \
    -file mesh.vtk \
    -bc none

