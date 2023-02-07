#!/bin/bash

set -e

rm -rf build
rm -rf gen

mkdir -p gen/kernels
mkdir -p gen/poisson/p_multigrid
mkdir -p gen/solvers/ls_utils
mkdir -p gen/matrices/poisson
mkdir -p gen/linear_solvers/petsc

python3 preprocessor.py 3

cd gen

# python3 $OP2_TRANSLATOR ins.cpp \
#         ins_data.cpp solver.cpp poisson/petsc/poisson.cpp \
#         ls/ls.cpp utils.cpp utils.cu ls/ls_reinit.cpp ls/ls_reinit.cu \
#         ls/ls_reinit_mpi.cpp ls/ls_reinit_mpi.cu ls/ls_reinit_mpi_naive.cpp \
#         ls/ls_reinit_mpi_naive.cu poisson/matrix/poisson_mat.cpp \
#         poisson/p_multigrid/p_multigrid.cpp kernels/

# sed -i '10i extern double reynolds;' openmp/ins_kernels.cpp

python3 $OP2_TRANSLATOR ins.cpp \
        solvers/advection_solver.cpp \
        solvers/ls_solver.cpp \
        matrices/poisson/poisson_mat.cpp \
        linear_solvers/petsc/poisson.cpp \
        solvers/mp_ins_solver.cpp \
        matrices/poisson/factor_poisson_mat.cpp \
        kernels/

cd ..

mkdir build

cd build

cmake .. \
  -DOP2_DIR=/home/u1717021/Code/PhD/OP2-Common/op2 \
  -DOPENBLAS_DIR=/opt/OpenBLAS \
  -DPETSC_DIR=/home/u1717021/Code/PhD/petsc/arch-linux-c-debug \
  -DPART_LIB_NAME=PARMETIS \
  -DPARMETIS_DIR=/home/u1717021/Code/PhD/ParMetis_Libs \
  -DOP2DGTOOLKIT_DIR=/home/u1717021/Code/PhD/OP2-DG-Toolkit/build \
  -DHDF5_DIR=/usr/local/module-software/hdf5-1.12.0-parallel \
  -DARMA_DIR=/home/u1717021/Code/PhD/armadillo-10.5.3/build \
  -DORDER=3 \
  -DBUILD_SN=ON \
  -DBUILD_CPU=ON \
  -DBUILD_MPI=ON \
  -DBUILD_GPU=ON

make -j 8
make
