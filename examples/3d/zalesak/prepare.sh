#!/bin/bash

gmsh \
    mesh.geo \
    -3 \
    -format vtk \
    -algo hxt \
    -nt 2 

mpirun -n 1 ${OP2_DG_TOOLKIT_DIR}/bin/vtk2hdf5_3D \
    -file mesh.vtk 

