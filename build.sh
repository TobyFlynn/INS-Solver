#!/bin/bash

set -e

rm -rf build
rm -rf gen

mkdir -p gen/kernels
mkdir -p gen/ls
mkdir -p gen/poisson
mkdir -p gen/io

python3 preprocessor.py 2

cd gen

python3 $OP2_TRANSLATOR ins.cpp \
        ins_data.cpp solver.cpp poisson/poisson.cpp \
        poisson/poisson_sub_mat.cpp timing.cpp ls/ls.cpp \
        poisson/poisson_cpu.cpp poisson/poisson_gpu.cpp \
        utils.cpp utils.cu ls/ls_reinit.cpp ls/ls_reinit.cu \
        io/save_solution.cpp io/save_solution_mpi.cpp kernels/

sed -i '10i extern double reynolds;' openmp/ins_kernels.cpp

cd ..

mkdir build

cd build

cmake .. \
  -DOP2_DIR=/dcs/pg20/u1717021/PhD/OP2-Common/op2 \
  -DCGNS_DIR=/dcs/pg20/u1717021/PhD/apps \
  -DOPENBLAS_DIR=/dcs/pg20/u1717021/PhD/apps \
  -DPETSC_DIR=$PETSC_DIR \
  -DPART_LIB_NAME=PARMETIS \
  -DPARMETIS_DIR=/dcs/pg20/u1717021/PhD/apps \
  -DOP2DGTOOLKIT_DIR=/dcs/pg20/u1717021/PhD/OP2-DG-Toolkit/build \
  -DARMA_DIR=/dcs/pg20/u1717021/PhD/apps \
  -DBUILD_SN=ON \
  -DBUILD_CPU=ON

make
