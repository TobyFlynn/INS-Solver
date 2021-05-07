project(INS-Solver)

# Require CMake 3.20
cmake_minimum_required(VERSION 3.20)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Load OP2 configuration (Set OP2_DIR to the directory containing OP2Config.cmake)
find_package(OP2 REQUIRED PATHS
  ${OP2-APPS_SOURCE_DIR}/../../op2/c/build
  ${PROJECT_SOURCE_DIR}/../../../op2/c/build
  ${CMAKE_INSTALL_PREFIX}/lib/op2)

# Import compiler flags for all build types
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OP2_CXX_FLAGS}")
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OP2_C_FLAGS}")
set(CUDA_NVCC_FLAGS "${CUDA_NVCC_FLAGS} -arch=sm_70 -Xcompiler -mno-float128")

include_directories(${OP2_INCLUDE_DIRS})
add_definitions(${OP2_USER_DEFINITIONS})

include(OP2Application)

set(PETSC_DIR $ENV{PETSC_DIR})

find_package(OpenBLAS)

enable_language(CUDA)
find_package(CUDAToolkit)
include_directories(${CUDAToolkit_INCLUDE_DIRS})
set(CMAKE_CUDA_COMPILER ${CUDAToolkit_NVCC_EXECUTABLE})

add_subdirectory(tools)
add_subdirectory(src)
