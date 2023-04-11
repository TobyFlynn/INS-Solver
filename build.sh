#!/bin/bash

set -e

rm -rf build
rm -rf gen2d
rm -rf gen3d

mkdir gen2d
cd gen2d

mkdir -p kernels
mkdir -p solvers/2d/ls_utils
mkdir -p solvers/3d/ls_utils
mkdir -p matrices/2d
mkdir -p matrices/3d/custom_kernels
mkdir -p linear_solvers/petsc_amg
mkdir -p linear_solvers/petsc_utils
mkdir -p linear_solvers/petsc_block_jacobi
mkdir -p linear_solvers/petsc_jacobi/custom_kernels
mkdir -p linear_solvers/pmultigrid/custom_kernels
mkdir -p linear_solvers/petsc_pmultigrid
mkdir -p linear_solvers/petsc_inv_mass

cd ..

mkdir gen3d
cd gen3d

mkdir -p kernels
mkdir -p solvers/2d/ls_utils
mkdir -p solvers/3d/ls_utils
mkdir -p matrices/2d
mkdir -p matrices/3d/custom_kernels
mkdir -p linear_solvers/petsc_amg
mkdir -p linear_solvers/petsc_utils
mkdir -p linear_solvers/petsc_block_jacobi
mkdir -p linear_solvers/petsc_jacobi/custom_kernels
mkdir -p linear_solvers/pmultigrid/custom_kernels
mkdir -p linear_solvers/petsc_pmultigrid
mkdir -p linear_solvers/petsc_inv_mass

cd ..

ORDER=3

python3 preprocessor.py 2 $ORDER

python3 preprocessor.py 3 $ORDER

cd gen2d

# python3 $OP2_TRANSLATOR ins.cpp \
#         ins_data.cpp solver.cpp poisson/petsc/poisson.cpp \
#         ls/ls.cpp utils.cpp utils.cu ls/ls_reinit.cpp ls/ls_reinit.cu \
#         ls/ls_reinit_mpi.cpp ls/ls_reinit_mpi.cu ls/ls_reinit_mpi_naive.cpp \
#         ls/ls_reinit_mpi_naive.cu poisson/matrix/poisson_mat.cpp \
#         poisson/p_multigrid/p_multigrid.cpp kernels/

# sed -i '10i extern double reynolds;' openmp/ins_kernels.cpp

python3 $OP2_TRANSLATOR ins2d.cpp \
        solvers/2d/advection_solver.cpp \
        solvers/2d/ls_solver.cpp \
        solvers/2d/mp_ins_solver.cpp \
        solvers/2d/ins_solver.cpp \
        solvers/2d/ce_solver.cpp \
        matrices/poisson_matrix.cpp \
        matrices/poisson_coarse_matrix.cpp \
        matrices/poisson_semi_matrix_free.cpp \
        matrices/poisson_matrix_free_diag.cpp \
        matrices/poisson_matrix_free.cpp \
        matrices/2d/poisson_matrix.cpp \
        matrices/2d/poisson_coarse_matrix.cpp \
        matrices/2d/poisson_semi_matrix_free.cpp \
        matrices/2d/poisson_matrix_free.cpp \
        matrices/2d/mm_poisson_matrix.cpp \
        matrices/2d/mm_poisson_matrix_free.cpp \
        matrices/2d/factor_poisson_matrix.cpp \
        matrices/2d/factor_poisson_coarse_matrix.cpp \
        matrices/2d/factor_poisson_semi_matrix_free.cpp \
        matrices/2d/factor_poisson_matrix_free.cpp \
        matrices/2d/factor_mm_poisson_matrix.cpp \
        matrices/2d/factor_mm_poisson_matrix_free.cpp \
        matrices/2d/cub_poisson_matrix.cpp \
        matrices/2d/cub_factor_poisson_matrix.cpp \
        matrices/2d/cub_mm_poisson_matrix.cpp \
        matrices/2d/cub_factor_mm_poisson_matrix.cpp \
        linear_solvers/petsc_block_jacobi/petsc_block_jacobi.cpp \
        linear_solvers/pmultigrid/pmultigrid.cpp \
        linear_solvers/petsc_inv_mass/petsc_inv_mass.cpp \
        linear_solvers/petsc_jacobi/petsc_jacobi.cpp \
        kernels/

sed -i "4i #include \"dg_compiler_defs.h\"" cuda/ins2d_kernels.cu
sed -i "4i #include \"dg_compiler_defs.h\"" openmp/ins2d_kernels.cpp
sed -i "4i #include \"dg_compiler_defs.h\"" seq/ins2d_seqkernels.cpp
sed -i "5i #include \"cblas.h\"" openmp/ins2d_kernels.cpp
sed -i "5i #include \"cblas.h\"" seq/ins2d_seqkernels.cpp

cd ..

cd gen3d

python3 $OP2_TRANSLATOR ins3d.cpp \
        solvers/3d/advection_solver.cpp \
        solvers/3d/ins_solver.cpp \
        solvers/3d/ls_solver.cpp \
        solvers/3d/mp_ins_solver.cpp \
        matrices/poisson_matrix.cpp \
        matrices/poisson_coarse_matrix.cpp \
        matrices/poisson_semi_matrix_free.cpp \
        matrices/poisson_matrix_free_diag.cpp \
        matrices/3d/poisson_matrix.cpp \
        matrices/3d/poisson_coarse_matrix.cpp \
        matrices/3d/poisson_semi_matrix_free.cpp \
        matrices/3d/poisson_matrix_free_mult.cpp \
        matrices/3d/poisson_matrix_free_diag.cpp \
        matrices/3d/mm_poisson_matrix.cpp \
        matrices/3d/mm_poisson_matrix_free.cpp \
        matrices/3d/factor_poisson_matrix.cpp \
        matrices/3d/factor_poisson_coarse_matrix.cpp \
        matrices/3d/factor_poisson_semi_matrix_free.cpp \
        matrices/3d/factor_poisson_matrix_free_diag.cpp \
        matrices/3d/factor_poisson_matrix_free_mult.cpp \
        matrices/3d/factor_mm_poisson_matrix.cpp \
        matrices/3d/factor_mm_poisson_semi_matrix_free.cpp \
        matrices/3d/factor_mm_poisson_matrix_free.cpp \
        matrices/3d/factor_mm_poisson_matrix_free_diag.cpp \
        linear_solvers/petsc_block_jacobi/petsc_block_jacobi.cpp \
        linear_solvers/pmultigrid/pmultigrid.cpp \
        linear_solvers/petsc_inv_mass/petsc_inv_mass.cpp \
        linear_solvers/petsc_jacobi/petsc_jacobi.cpp \
        kernels/

sed -i "4i #include \"dg_compiler_defs.h\"" cuda/ins3d_kernels.cu
sed -i "4i #include \"dg_compiler_defs.h\"" openmp/ins3d_kernels.cpp
sed -i "4i #include \"dg_compiler_defs.h\"" seq/ins3d_seqkernels.cpp
sed -i "5i #include \"liquid_whistle_consts.h\"" cuda/ins3d_kernels.cu
sed -i "5i #include \"liquid_whistle_consts.h\"" openmp/ins3d_kernels.cpp
sed -i "5i #include \"liquid_whistle_consts.h\"" seq/ins3d_seqkernels.cpp
sed -i "6i #include \"cblas.h\"" openmp/ins3d_kernels.cpp
sed -i "6i #include \"cblas.h\"" seq/ins3d_seqkernels.cpp
sed -i "73i #include \"../matrices/3d/custom_kernels/custom_fpmf_3d_mult_faces_flux.cu\"" cuda/ins3d_kernels.cu
sed -i "73i #include \"../matrices/3d/custom_kernels/custom_pmf_3d_mult_faces_flux.cu\"" cuda/ins3d_kernels.cu

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
  -DINIPP_DIR=/home/u1717021/Code/PhD/inipp/inipp \
  -DORDER=$ORDER \
  -DBUILD_SN=ON \
  -DBUILD_CPU=ON \
  -DBUILD_MPI=ON \
  -DBUILD_GPU=ON

make -j 8
make
