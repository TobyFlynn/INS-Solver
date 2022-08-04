# Require CMake 3.18
cmake_minimum_required(VERSION 3.18)
project(INS-Solver CXX)

enable_language(C)

# Import compiler flags for all build types
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++14 -fopenmp")
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS}")

# Default OP2 location for me
if(NOT OP2_DIR)
  message(FATAL_ERROR "OP2_DIR is not specified")
endif()
include_directories(${OP2_DIR}/include)

# CGNS library for me
if(NOT CGNS_DIR)
  message(FATAL_ERROR "CGNS_DIR is not specified")
endif()
include_directories(${CGNS_DIR}/include)
set(CGNS_LIB "-L${CGNS_DIR}/lib -lcgns")

# OpenBLAS library for me
if(NOT OPENBLAS_DIR)
  message(FATAL_ERROR "OPENBLAS_DIR is not specified")
endif()
include_directories(${OPENBLAS_DIR}/include)
set(OPENBLAS_LIB "-L${OPENBLAS_DIR}/lib -lopenblas")

if(NOT PETSC_DIR)
  message(FATAL_ERROR "PETSC_DIR is not specified")
endif()
include_directories(${PETSC_DIR}/include)
set(PETSC_LIB "-L${PETSC_DIR}/lib -lpetsc")

if(BUILD_MPI)
  if(NOT PART_LIB_NAME)
    message(FATAL_ERROR "PART_LIB_NAME is not specified")
  endif()

  if(PART_LIB_NAME STREQUAL "PARMETIS")
    if(NOT PARMETIS_DIR)
      message(FATAL_ERROR "PARMETIS_DIR is not specified")
    endif()
    include_directories(${PARMETIS_DIR}/include)
    set(PART_LIB "-L${PARMETIS_DIR}/lib -lparmetis -lmetis -lGKlib")
  elseif(PART_LIB_NAME STREQUAL "PTSCOTCH")
    if(NOT PTSCOTCH_DIR)
      message(FATAL_ERROR "PTSCOTCH_DIR is not specified")
    endif()
    include_directories(${PTSCOTCH_DIR}/include)
    set(PART_LIB "-L${PTSCOTCH_DIR}/lib -lptscotch -lscotch -lptscotcherr -lscotcherr -lptscotcherrexit -lscotcherrexit")
  else()
    message(FATAL_ERROR "Unrecognised OP2 partitioning library specified in PART_LIB_NAME")
  endif()
endif()

if(NOT OP2DGTOOLKIT_DIR)
  message(FATAL_ERROR "OP2DGTOOLKIT_DIR is not specified")
endif()
include_directories(${OP2DGTOOLKIT_DIR}/include)

if(NOT ARMA_DIR)
  message(FATAL_ERROR "ARMA_DIR is not specified")
endif()
include_directories(${ARMA_DIR}/include)
set(ARMA_LIB "-L${ARMA_DIR}/lib -larmadillo")

if(BUILD_GPU)
  enable_language(CUDA)
  find_package(CUDAToolkit)
  include_directories(${CUDAToolkit_INCLUDE_DIRS})
  set(CMAKE_CUDA_COMPILER ${CUDAToolkit_NVCC_EXECUTABLE})
  set(CMAKE_CUDA_ARCHITECTURES 70 75)
  set(CMAKE_CUDA_FLAGS "-Xcompiler -mno-float128")
endif()

include_directories(${PROJECT_SOURCE_DIR}/include)

add_subdirectory(tools)
add_subdirectory(gen)
