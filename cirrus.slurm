#!/bin/bash
#
#SBATCH --partition=gpu
#SBATCH --qos=gpu
#SBATCH --gres=gpu:4
#SBATCH --nodes=6
#SBATCH --time=00:15:00
#SBATCH --exclusive

module load openmpi/4.1.4-cuda-11.8
module load nvidia/nvhpc-nompi/22.11

export OMP_PLACES=cores
export OMPI_MCA_pml=ob1
export OMP_NUM_THREADS=10
export OPENBLAS_NUM_THREADS=1

ins_dir=/path/to/ins/executables
grid_dir=/path/to/grids

srun --ntasks=24 --cpus-per-task=10 --hint=nomultithread $ins_dir/ins3d_mpi_cuda \
	-config ./cirrus.ini \
	-input $grid_dir/re_1000.h5 \
	-gpudirect
