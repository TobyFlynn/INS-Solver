CONFIG_AR := ar rcs
CC = mpicc
CXX = mpicxx
NVCC = nvcc
HIPCC = hipcc

CUDA_ARCH = 70
HIP_ARCH = gfx908

OPENMP_FLAG := -fopenmp

OP2_DIR = /dcs/pg20/u1717021/PhD/OP2-My-Fork/op2
OPENBLAS_DIR = /dcs/pg20/u1717021/PhD/apps
ARMA_DIR = /dcs/pg20/u1717021/PhD/apps
HDF5_DIR = /dcs/pg20/u1717021/PhD/apps
PETSC_DIR = /dcs/pg20/u1717021/PhD/petsc-install
INIPP_DIR = /dcs/pg20/u1717021/PhD/inipp/inipp
HIGHFIVE_DIR = /dcs/pg20/u1717021/PhD/HighFive/include
MPI_DIR = /dcs/pg20/u1717021/PhD/apps
OP2_DG_TOOLKIT_DIR = /dcs/pg20/u1717021/PhD/OP2-DG-Toolkit
#HYPRE_DIR =

PART_LIB_NAME = PARMETIS
PARMETIS_DIR = /dcs/pg20/u1717021/PhD/apps
#PTSCOTCH_DIR = /dcs/pg20/u1717021/PhD/apps
PARTITION_LIB = -L$(PARMETIS_DIR)/lib -lparmetis -lmetis -lGKlib
#PARTITION_LIB = -L${PTSCOTCH_DIR}/lib -lptscotch -lscotch -lptscotcherr -lscotcherr -lptscotcherrexit -lscotcherrexit
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
INIPP_INC = -I$(INIPP_DIR)
MPI_INC = -I$(MPI_DIR)/include
OP2_DG_TOOLKIT_INC = -I$(OP2_DG_TOOLKIT_DIR)/include
HYPRE_INC = -I$(HYPRE_DIR)/include

OPENBLAS_LIB = -L$(OPENBLAS_DIR)/lib -lopenblas
PETSC_LIB = -L$(PETSC_DIR)/lib -lpetsc
ARMA_LIB = -L$(ARMA_DIR)/lib64 -larmadillo
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