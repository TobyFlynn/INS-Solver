set(COMMON_SRC solver_op.cpp poisson_op.cpp poisson_sub_mat_op.cpp ins_data_op.cpp timing_op.cpp shared_save_functions.cpp ls_op.cpp)
set(CPU_SRC poisson_cpu_op.cpp utils_op.cpp)
set(GPU_SRC poisson_gpu_op.cu utils_op.cu)
set(SN_SRC save_solution_op.cpp load_mesh.cpp)
set(MPI_SRC save_solution_mpi_op.cpp load_mesh_mpi.cpp mpi_helper_func.cpp)

set(GFORT_LIB /gpfs/users/shared/apps/core/gcc/7.3.0/lib64/libgfortran.so.4)

add_executable(ins_openmp ins_op.cpp ${COMMON_SRC} ${CPU_SRC} ${SN_SRC} openmp/ins_kernels.cpp)
target_link_libraries(ins_openmp -L${OP2_DIR}/lib -lop2_openmp -lop2_hdf5 libcgns.so libhdf5.so ${PETSC_DIR}/lib/libpetsc.so -L${ARMA_DIR}/lib64 -larmadillo ${OpenBLAS_LIBRARIES} ${OP2DGTOOLKIT_DIR}/lib/libop2dgtoolkit_openmp.so ${GFORT_LIB})
target_compile_definitions(ins_openmp PRIVATE DG_ORDER=3)

add_executable(ins_cuda ins_op.cpp ${COMMON_SRC} ${GPU_SRC} ${SN_SRC} cuda/ins_kernels.cu)
target_link_libraries(ins_cuda -L${OP2_DIR}/lib -lop2_cuda -lop2_hdf5 libcgns.so libhdf5.so ${PETSC_DIR}/lib/libpetsc.so -L${ARMA_DIR}/lib64 -larmadillo ${OpenBLAS_LIBRARIES} ${OP2DGTOOLKIT_DIR}/lib/libop2dgtoolkit_cuda.so ${GFORT_LIB} CUDA::cudart_static CUDA::cublas)
target_compile_definitions(ins_cuda PRIVATE BUILD_CUDA OP2_DG_CUDA INS_CUDA DG_ORDER=3)

add_executable(ins_mpi_openmp ins_op.cpp ${COMMON_SRC} ${CPU_SRC} ${MPI_SRC} openmp/ins_kernels.cpp)
target_link_libraries(ins_mpi_openmp -L${OP2_DIR}/lib -lop2_mpi -lop2_openmp -lop2_hdf5 libcgns.so libhdf5.so ${PETSC_DIR}/lib/libpetsc.so -L${ARMA_DIR}/lib64 -larmadillo ${OpenBLAS_LIBRARIES} ${OP2DGTOOLKIT_DIR}/lib/libop2dgtoolkit_mpi_openmp.so ${GFORT_LIB} libptscotch.a libscotch.a)
target_compile_definitions(ins_mpi_openmp PRIVATE INS_MPI DG_ORDER=3)

add_executable(ins_mpi_cuda ins_op.cpp ${COMMON_SRC} ${GPU_SRC} ${MPI_SRC} cuda/ins_kernels.cu)
target_link_libraries(ins_mpi_cuda -L${OP2_DIR}/lib -lop2_mpi_cuda -lop2_hdf5 libcgns.so libhdf5.so ${PETSC_DIR}/lib/libpetsc.so -L${ARMA_DIR}/lib64 -larmadillo ${OpenBLAS_INCLUDE_DIRS} ${OP2DGTOOLKIT_DIR}/lib/libop2dgtoolkit_mpi_cuda.so ${GFORT_LIB} CUDA::cudart_static CUDA::cublas libptscotch.a libscotch.a)
target_compile_definitions(ins_mpi_cuda PRIVATE BUILD_CUDA OP2_DG_CUDA INS_CUDA INS_MPI DG_ORDER=3)
