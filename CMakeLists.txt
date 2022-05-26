# Require CMake 3.18
cmake_minimum_required(VERSION 3.18)
project(INS-Solver CXX)

# Import compiler flags for all build types
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++14 -fopenmp")
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS}")

# Default OP2 location for me
if(NOT OP2_DIR)
  set(OP2_DIR "/gpfs/users/tsf79-sxa03/OP2-Common/op2")
endif()
include_directories(${OP2_DIR}/include)

include_directories(${CGNS_DIR}/include)

# OpenBLAS library for me
find_package(OpenBLAS)
include_directories(${OpenBLAS_INCLUDE_DIRS})

set(PETSC_DIR $ENV{PETSC_DIR}/arch-linux-c-debug)
include_directories($ENV{PETSC_DIR}/include ${PETSC_DIR}/include)

if(NOT OP2DGTOOLKIT_DIR)
  set(OP2DGTOOLKIT_DIR "/gpfs/users/tsf79-sxa03/OP2-DG-Toolkit/build")
endif()
include_directories(${OP2DGTOOLKIT_DIR}/include)

if(NOT ARMA_DIR)
  set(ARMA_DIR "/gpfs/users/tsf79-sxa03/armadillo-10.6.1/build")
endif()
include_directories(${ARMA_DIR}/include)

enable_language(CUDA)
find_package(CUDAToolkit)
include_directories(${CUDAToolkit_INCLUDE_DIRS})
set(CMAKE_CUDA_COMPILER ${CUDAToolkit_NVCC_EXECUTABLE})
set(CMAKE_CUDA_ARCHITECTURES 70)

include_directories(${PROJECT_SOURCE_DIR}/include)

add_subdirectory(tools)
add_subdirectory(gen)
