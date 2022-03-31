# Require CMake 3.18
cmake_minimum_required(VERSION 3.18)
project(INS-Solver CXX)

# Import compiler flags for all build types
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++14 -fopenmp")
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS}")

# Default OP2 location for me
if(NOT OP2_DIR)
  set(OP2_DIR "/home/u1717021/Code/PhD/OP2-Common/op2")
endif()
include_directories(${OP2_DIR}/include)

# CGNS library for me
if(NOT CGNS_DIR)
  set(CGNS_DIR "/home/u1717021/Code/PhD/CGNS/build")
endif()
include_directories(${CGNS_DIR}/include)

# OpenBLAS library for me
if(NOT OPENBLAS_DIR)
  set(OPENBLAS_DIR "/opt/OpenBLAS")
endif()
include_directories(${OPENBLAS_DIR}/include)

# HYPRE library for me
if(NOT HYPRE_DIR)
  set(HYPRE_DIR "/home/u1717021/Code/PhD/hypre/src/hypre")
endif()
include_directories(${HYPRE_DIR}/include)

set(PETSC_DIR $ENV{PETSC_DIR})
include_directories(${PETSC_DIR}/include ${PETSC_DIR}/arch-linux-c-debug/include)

set(PARMETIS_DIR "/home/u1717021/Code/PhD/parmetis-4.0.3/build")
include_directories(${PARMETIS_DIR}/include)

if(NOT OP2DGTOOLKIT_DIR)
  set(OP2DGTOOLKIT_DIR "/home/u1717021/Code/PhD/OP2-DG-Toolkit/build")
endif()
include_directories(${OP2DGTOOLKIT_DIR}/include)

if(NOT ARMA_DIR)
  set(ARMA_DIR "/home/u1717021/Code/PhD/armadillo-10.5.3/build")
endif()
include_directories(${ARMA_DIR}/include)

enable_language(CUDA)
find_package(CUDAToolkit)
include_directories(${CUDAToolkit_INCLUDE_DIRS})
set(CMAKE_CUDA_COMPILER ${CUDAToolkit_NVCC_EXECUTABLE})
set(CMAKE_CUDA_ARCHITECTURES 75)

include_directories(${PROJECT_SOURCE_DIR}/include)

add_subdirectory(tools)
add_subdirectory(gen)
