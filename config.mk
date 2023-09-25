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

PART_LIB_NAME = PARMETIS
PARMETIS_DIR = /dcs/pg20/u1717021/PhD/apps
#PTSCOTCH_DIR = /dcs/pg20/u1717021/PhD/apps
PARTITION_LIB = -L$(PARMETIS_DIR)/lib -lparmetis -lmetis -lGKlib
#PARTITION_LIB = -L${PTSCOTCH_DIR}/lib -lptscotch -lscotch -lptscotcherr -lscotcherr -lptscotcherrexit -lscotcherrexit
#PARTITION_LIB = -L$(PARMETIS_DIR)/lib -lparmetis -lmetis -lGKlib -L${PTSCOTCH_DIR}/lib -lptscotch -lscotch -lptscotcherr -lscotcherr -lptscotcherrexit -lscotcherrexit

CUBLAS_LIB = -lcublas
HIPBLAS_LIB = -lhipblas

ORDER = 3
SOA = true