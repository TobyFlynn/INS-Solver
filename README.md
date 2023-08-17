# INS-Solver
This multiphase incompressible Navier-Stokes solver has the following dependencies:
- [OP2 on its force halo compute branch](https://github.com/OP-DSL/OP2-Common/tree/feature/force_halo_compute) built with [PT-Scotch](https://gitlab.inria.fr/scotch/scotch) and/or [ParMETIS](http://glaros.dtc.umn.edu/gkhome/metis/parmetis/download), [HDF5](https://www.hdfgroup.org/solutions/hdf5/), MPI and CUDA (only required for GPU builds).
- [OP2-DG-Toolkit]() and its dependencies

Once the dependencies have been built, the `build.sh` script can be used to build the solver. The CMake command at the end of the script must be modified to include the locations of where the dependencies have been installed, like for OP2-DG-Toolkit.

## Running the solver
The solver takes in a configuration file and a HDF5 mesh file. In addition, GPUDirect can be utilised by OP2 if the `-gpudirect` flag is passed. So an MPI+CUDA run would be started with the following command: `mpirun -n 2 ./ins_mpi_cuda -input ./mesh.h5 -config config.ini -gpudirect`

## ARCHER2
The build script used on ARCHER2 is `build-archer.sh`. An example configuration file is `archer.ini`. An example SLURM submission script is `archer.slurm`. Some of the dependencies were provided by the clusters module environment and the full module environment used was:
- cce/15.0.0
- cray-libsci/22.12.1.1
- cray-ucx/2.7.0-1
- load-epcc-module
- scotch/7.0.3
- craype-x86-rome
- craype/2.7.19
- PrgEnv-cray/8.3.3
- craype-network-ucx
- bolt/0.8
- metis/5.1.0
- xpmem/2.5.2-2.4_3.30__gd0f7936.shasta
- cray-dsmml/0.2.2
- cray-hdf5-parallel/1.12.2.1
- cray-mpich-ucx/8.1.23
- epcc-setup-env
- parmetis/4.0.3

## Cirrus
The build script used on Cirrus is `build-cirrus.sh`. An example configuration file is `cirrus.ini`. An example SLURM submission script is `cirrus.slurm`. Some of the dependencies were provided by the clusters module environment and the full module environment used was:
- git/2.37.3
- /mnt/lustre/indy2lfs/sw/modulefiles/epcc/setup-env
- openmpi/4.1.4-cuda-11.8
- cmake/3.25.2
- flex/2.6.4
- epcc/utils
- gcc/8.2.0(default)
- nvidia/nvhpc-nompi/22.11
- bison/3.6.4
