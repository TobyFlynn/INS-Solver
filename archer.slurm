#!/bin/bash

# Slurm job options (job-name, compute nodes, job time)
#SBATCH --job-name=DG-INS
#SBATCH --time=00:15:00
#SBATCH --nodes=8
#SBATCH --ntasks-per-node=128
#SBATCH --cpus-per-task=1
#SBATCH --partition=standard
#SBATCH --qos=standard

# Modules
module load PrgEnv-cray
module load cray-hdf5-parallel
module load scotch
module load parmetis
module swap craype-network-ofi craype-network-ucx
module swap cray-mpich cray-mpich-ucx
module load xpmem/2.5.2-2.4_3.30__gd0f7936.shasta
module list 

# Set the number of threads to 1
#   This prevents any threaded system libraries from automatically 
#   using threading.
export OMP_NUM_THREADS=1

ins_dir=/path/to/ins/executables
grid_dir=/path/to/grids

srun --distribution=block:block --hint=nomultithread $ins_dir/ins3d_mpi \
	-config ./archer.ini \
	-input $grid_dir/liquid_whistle.h5
