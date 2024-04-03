CONFIG_AR := ar rcs
CC = mpicc
CXX = mpicxx
NVCC = nvcc
HIPCC = hipcc

CUDA_ARCH = 70
HIP_ARCH = gfx908

OPENMP_FLAG := -fopenmp

OP2_DIR := /home/toby/Code/PhD/OP2-Fork/op2
OPENBLAS_DIR := $(shell spack location -i openblas)
ARMA_DIR := $(shell spack location -i armadillo)
HDF5_DIR := $(shell spack location -i hdf5)
PETSC_DIR := $(shell spack location -i petsc)
INIPP_DIR := $(shell spack location -i inipp)
MPI_DIR := $(shell spack location -i openmpi)
OP2_DG_TOOLKIT_DIR = /home/toby/Code/PhD/OP2-DG-Toolkit
CDT_DIR := $(shell spack location -i cdt)
#HYPRE_DIR =

PART_LIB_NAME = PTSCOTCH
#PARMETIS_DIR = /dcs/pg20/u1717021/PhD/apps
PTSCOTCH_DIR = $(shell spack location -i scotch)
#PARTITION_LIB = -L$(PARMETIS_DIR)/lib -lparmetis -lmetis -lGKlib
PARTITION_LIB = -L${PTSCOTCH_DIR}/lib -lptscotch -lscotch -lptscotcherr -lscotcherr -lptscotcherrexit -lscotcherrexit
#PARTITION_LIB = -L$(PARMETIS_DIR)/lib -lparmetis -lmetis -lGKlib -L${PTSCOTCH_DIR}/lib -lptscotch -lscotch -lptscotcherr -lscotcherr -lptscotcherrexit -lscotcherrexit

CUBLAS_LIB = -lcublas
HIPBLAS_LIB = -lhipblas

ORDER = 3
SOA = 1
BUILD_WITH_HYPRE = 0

# Probably do not need to change derived variables below this comment unless
# dependencies were installed in unusual locations

EXTRA_LIBS_CPU = 
EXTRA_LIBS_CUDA = 
EXTRA_LIBS_HIP = 

OP2_INC = -I$(OP2_DIR)/include
OPENBLAS_INC = -I$(OPENBLAS_DIR)/include
ARMA_INC = -I$(ARMA_DIR)/include
PETSC_INC = -I$(PETSC_DIR)/include
INIPP_INC = -I$(INIPP_DIR)/include
MPI_INC = -I$(MPI_DIR)/include
OP2_DG_TOOLKIT_INC = -I$(OP2_DG_TOOLKIT_DIR)/include
HYPRE_INC = -I$(HYPRE_DIR)/include
CDT_INC = -I$(CDT_DIR)/include

OPENBLAS_LIB = -L$(OPENBLAS_DIR)/lib -lopenblas
PETSC_LIB = -L$(PETSC_DIR)/lib -lpetsc
ARMA_LIB = -L$(ARMA_DIR)/lib -larmadillo
HDF5_LIB = -L$(HDF5_DIR)/lib -lhdf5
MPI_LIB = -L$(MPI_DIR)/lib -lmpi
HYPRE_LIB = -L$(HYPRE_DIR)/lib -lHYPRE

OP2_OPENMP_LIBS = -L$(OP2_DIR)/lib -lop2_openmp -lop2_hdf5
OP2_CUDA_LIBS = -L$(OP2_DIR)/lib -lop2_cuda -lop2_hdf5
OP2_HIP_LIBS = -L$(OP2_DIR)/lib -lop2_hip -lop2_hdf5
OP2_MPI_LIBS = -L$(OP2_DIR)/lib -lop2_mpi
OP2_MPI_CUDA_LIBS = -L$(OP2_DIR)/lib -lop2_mpi_cuda
OP2_MPI_HIP_LIBS = -L$(OP2_DIR)/lib -lop2_mpi_hip

OP2_TRANSLATOR = $(OP2_DIR)/../translator/c/op2.py
