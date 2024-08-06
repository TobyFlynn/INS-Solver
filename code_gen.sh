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
        solvers/2d/advec_diff_solver.cpp \
        solvers/2d/ls_solver.cpp \
        solvers/2d/ins_solver_base.cpp \
        solvers/2d/ins_solver.cpp \
        solvers/2d/ins_temperature_solver.cpp \
        solvers/2d/mp_ins_solver.cpp \
        solvers/2d/ls_driver.cpp \
        measurements/2d/lift_drag_cylinder.cpp \
        measurements/2d/l2_error_vortex.cpp \
        measurements/2d/min_max_interface.cpp \
        measurements/2d/min_max_pressure.cpp \
        measurements/2d/ls_advec_error.cpp \
        measurements/2d/mass_of_phases.cpp \
        solvers/2d/slip_matrix/viscous_matrix.cpp \
        solvers/2d/slip_matrix/factor_viscous_matrix.cpp \
        solvers/2d/slip_matrix/viscous_solver.cpp \
        kernels/

cuda_line_no_2d=$(grep -n op_cuda_reduction cuda/ins2d_kernels.cu | cut -d : -f 1)
hip_line_no_2d=$(grep -n op_hip_reduction hip/ins2d_kernels.cpp | cut -d : -f 1)
openmp_line_no_2d=$(grep -n op_lib_cpp openmp/ins2d_kernels.cpp | cut -d : -f 1)
seq_line_no_2d=$(grep -n op_lib_cpp seq/ins2d_seqkernels.cpp | cut -d : -f 1)

cuda_line_no_2d=$((cuda_line_no_2d+1))
hip_line_no_2d=$((hip_line_no_2d+1))
openmp_line_no_2d=$((openmp_line_no_2d+1))
seq_line_no_2d=$((seq_line_no_2d+1))

text_gpu_2d="#include \"dg_compiler_defs.h\"\n#include \"problem_specific_2d.h\"\n#include \"dg_global_constants/dg_mat_constants_2d.h\""
text_cpu_2d="#include \"dg_compiler_defs.h\"\n#include \"problem_specific_2d.h\"\n#include \"dg_global_constants/dg_mat_constants_2d.h\"\n#include \"cblas.h\""

sed -i "${cuda_line_no_2d}i $text_gpu_2d" cuda/ins2d_kernels.cu
sed -i "${hip_line_no_2d}i $text_gpu_2d" hip/ins2d_kernels.cpp
sed -i "${openmp_line_no_2d}i $text_cpu_2d" openmp/ins2d_kernels.cpp
sed -i "${seq_line_no_2d}i $text_cpu_2d" seq/ins2d_seqkernels.cpp

cd ..

cd gen_3d

python3 $OP2_TRANSLATOR ins3d.cpp \
        solvers/3d/advection_solver.cpp \
        solvers/3d/diffusion_solver.cpp \
        solvers/3d/ins_solver_base.cpp \
        solvers/3d/ins_solver.cpp \
        solvers/3d/ls_solver.cpp \
        solvers/3d/mp_ins_solver.cpp \
        solvers/3d/ls_driver.cpp \
        measurements/3d/enstrophy.cpp \
        solvers/3d/slip_matrix/viscous_matrix.cpp \
        solvers/3d/slip_matrix/factor_viscous_matrix.cpp \
        solvers/3d/slip_matrix/viscous_solver.cpp \
        kernels/

cuda_line_no_3d=$(grep -n op_cuda_reduction cuda/ins3d_kernels.cu | cut -d : -f 1)
hip_line_no_3d=$(grep -n op_hip_reduction hip/ins3d_kernels.cpp | cut -d : -f 1)
openmp_line_no_3d=$(grep -n op_lib_cpp openmp/ins3d_kernels.cpp | cut -d : -f 1)
seq_line_no_3d=$(grep -n op_lib_cpp seq/ins3d_seqkernels.cpp | cut -d : -f 1)

cuda_line_no_3d=$((cuda_line_no_3d+1))
hip_line_no_3d=$((hip_line_no_3d+1))
openmp_line_no_3d=$((openmp_line_no_3d+1))
seq_line_no_3d=$((seq_line_no_3d+1))

text_gpu_3d="#include \"dg_compiler_defs.h\"\n#include \"problem_specific_3d.h\"\n#include \"dg_global_constants/dg_mat_constants_3d.h\""
text_cpu_3d="#include \"dg_compiler_defs.h\"\n#include \"problem_specific_3d.h\"\n#include \"dg_global_constants/dg_mat_constants_3d.h\"\n#include \"cblas.h\""

sed -i "${cuda_line_no_3d}i $text_gpu_3d" cuda/ins3d_kernels.cu
sed -i "${hip_line_no_3d}i $text_gpu_3d" hip/ins3d_kernels.cpp
sed -i "${openmp_line_no_3d}i $text_cpu_3d" openmp/ins3d_kernels.cpp
sed -i "${seq_line_no_3d}i $text_cpu_3d" seq/ins3d_seqkernels.cpp

cd ../..
