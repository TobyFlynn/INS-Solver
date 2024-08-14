#!/bin/bash

gmsh \
    m0.geo \
    -2 \
    -format vtk \
    -algo front2d \
    -nt 2 

mpirun -n 1 ${OP2_DG_TOOLKIT_DIR}/bin/vtk2hdf5_2D \
    -file m0.vtk \
    -bc euler_vortex

mv mesh.h5 m0.h5

gmsh \
    m1.geo \
    -2 \
    -format vtk \
    -algo front2d \
    -nt 2 

mpirun -n 1 ${OP2_DG_TOOLKIT_DIR}/bin/vtk2hdf5_2D \
    -file m1.vtk \
    -bc euler_vortex

mv mesh.h5 m1.h5

gmsh \
    m2.geo \
    -2 \
    -format vtk \
    -algo front2d \
    -nt 2 

mpirun -n 1 ${OP2_DG_TOOLKIT_DIR}/bin/vtk2hdf5_2D \
    -file m2.vtk \
    -bc euler_vortex

mv mesh.h5 m2.h5

gmsh \
    m3.geo \
    -2 \
    -format vtk \
    -algo front2d \
    -nt 2 

mpirun -n 1 ${OP2_DG_TOOLKIT_DIR}/bin/vtk2hdf5_2D \
    -file m3.vtk \
    -bc euler_vortex

mv mesh.h5 m3.h5