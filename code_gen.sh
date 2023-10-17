#!/bin/bash

set -e

dir_list=$(find src -maxdepth 3 -mindepth 1 | grep -v "\." | sed "s+src/++g")

mkdir -p code_gen/gen_2d
cd code_gen/gen_2d

for i in $dir_list
do
  mkdir -p $i
done

cd ..

mkdir gen_3d
cd gen_3d

for i in $dir_list
do
  mkdir -p $i
done

cd ../..

python3 preprocessor.py 2 $ORDER

python3 preprocessor.py 3 $ORDER

if [ $SOA = 1 ]; then
  export OP_AUTO_SOA=1
fi

cd code_gen/gen_2d

python3 $OP2_TRANSLATOR ins2d.cpp \
        solvers/2d/advection_solver.cpp \
        solvers/2d/diffusion_solver.cpp \
        solvers/2d/ls_solver.cpp \
        solvers/2d/ins_solver_base.cpp \
        solvers/2d/ins_solver.cpp \
        solvers/2d/mp_ins_solver.cpp \
        kernels/

sed -i "4i #include \"dg_compiler_defs.h\"" cuda/ins2d_kernels.cu
sed -i "4i #include \"dg_compiler_defs.h\"" hip/ins2d_kernels.cpp
sed -i "4i #include \"dg_compiler_defs.h\"" openmp/ins2d_kernels.cpp
sed -i "4i #include \"dg_compiler_defs.h\"" seq/ins2d_seqkernels.cpp
sed -i "6i #include \"dg_global_constants/dg_mat_constants_2d.h\"" cuda/ins2d_kernels.cu
sed -i "6i #include \"dg_global_constants/dg_mat_constants_2d.h\"" hip/ins2d_kernels.cpp
sed -i "6i #include \"dg_global_constants/dg_mat_constants_2d.h\"" openmp/ins2d_kernels.cpp
sed -i "6i #include \"dg_global_constants/dg_mat_constants_2d.h\"" seq/ins2d_seqkernels.cpp
sed -i "5i #include \"cblas.h\"" openmp/ins2d_kernels.cpp
sed -i "5i #include \"cblas.h\"" seq/ins2d_seqkernels.cpp

cd ..

cd gen_3d

python3 $OP2_TRANSLATOR ins3d.cpp \
        solvers/3d/advection_solver.cpp \
        solvers/3d/ins_solver_base.cpp \
        solvers/3d/ins_solver.cpp \
        solvers/3d/ls_solver.cpp \
        solvers/3d/mp_ins_solver.cpp \
        kernels/

sed -i "4i #include \"dg_compiler_defs.h\"" cuda/ins3d_kernels.cu
sed -i "4i #include \"dg_compiler_defs.h\"" hip/ins3d_kernels.cpp
sed -i "4i #include \"dg_compiler_defs.h\"" openmp/ins3d_kernels.cpp
sed -i "4i #include \"dg_compiler_defs.h\"" seq/ins3d_seqkernels.cpp
sed -i "5i #include \"liquid_whistle_consts.h\"" cuda/ins3d_kernels.cu
sed -i "5i #include \"liquid_whistle_consts.h\"" hip/ins3d_kernels.cpp
sed -i "5i #include \"liquid_whistle_consts.h\"" openmp/ins3d_kernels.cpp
sed -i "5i #include \"liquid_whistle_consts.h\"" seq/ins3d_seqkernels.cpp
sed -i "6i #include \"dg_global_constants/dg_mat_constants_3d.h\"" cuda/ins3d_kernels.cu
sed -i "6i #include \"dg_global_constants/dg_mat_constants_3d.h\"" hip/ins3d_kernels.cpp
sed -i "6i #include \"dg_global_constants/dg_mat_constants_3d.h\"" openmp/ins3d_kernels.cpp
sed -i "6i #include \"dg_global_constants/dg_mat_constants_3d.h\"" seq/ins3d_seqkernels.cpp
sed -i "7i #include \"cblas.h\"" openmp/ins3d_kernels.cpp
sed -i "7i #include \"cblas.h\"" seq/ins3d_seqkernels.cpp

cd ../..
