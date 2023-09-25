#!/bin/bash

set -e

rm -rf build

./code_gen.sh

mkdir build

cd build

#-DAMGX_DIR=/home/u1717021/Code/PhD/AMGX-install \
#-DHYPRE_DIR=/home/u1717021/Code/PhD/hypre-install \

cmake .. \
  -DOP2_DIR=/home/u1717021/Code/PhD/OP2-My-Fork/op2 \
  -DOPENBLAS_DIR=/opt/OpenBLAS \
  -DPETSC_DIR=/home/u1717021/Code/PhD/petsc/arch-linux-c-debug \
  -DPART_LIB_NAME=PARMETIS \
  -DPARMETIS_DIR=/home/u1717021/Code/PhD/ParMetis_Libs \
  -DOP2DGTOOLKIT_DIR=/home/u1717021/Code/PhD/OP2-DG-Toolkit/build \
  -DHDF5_DIR=/usr/local/module-software/hdf5-1.12.0-parallel \
  -DARMA_DIR=/home/u1717021/Code/PhD/armadillo-10.5.3/build \
  -DINIPP_DIR=/home/u1717021/Code/PhD/inipp/inipp \
  -DHYPRE_DIR=/home/u1717021/Code/PhD/hypre-cpu-sp-install \
  -DORDER=$ORDER \
  -DSOA=$SOA \
  -DBUILD_SN=ON \
  -DBUILD_CPU=OFF \
  -DBUILD_MPI=ON \
  -DBUILD_GPU=ON

make -j 16
make
