set(OPENBLAS_SRC openBLAS_op/init_gauss_grad_op.cpp openBLAS_op/init_gauss_grad_neighbour_op.cpp)
set(CUBLAS_SRC cuBLAS_op/init_gauss_grad_op.cu cuBLAS_op/init_gauss_grad_neighbour_op.cu)

set(COMMON_SRC solver_op.cpp ins_data_op.cpp timing_op.cpp ls_op.cpp poisson_op.cpp shared_save_functions.cpp)
set(CPU_SRC utils_op.cpp poisson_cpu_op.cpp)
set(GPU_SRC utils_op.cu poisson_gpu_op.cu)

set(SN_SRC save_solution_op.cpp load_mesh.cpp)
set(MPI_SRC save_solution_mpi_op.cpp load_mesh_mpi.cpp mpi_helper_func.cpp)

op2_application(ins_openmp LIBS op2_openmp op2_hdf5 SOURCES ins_op.cpp ${COMMON_SRC} ${CPU_SRC} ${SN_SRC} ${OPENBLAS_SRC} openmp/ins_kernels.cpp)
target_include_directories(ins_openmp PUBLIC ${CGNS_DIR}/include ${ARMA_DIR}/include ${PETSC_DIR}/include ${PETSC_DIR}/arch-linux-c-debug/include ${OPENBLAS_DIR}/include ${OP2DGTOOLKIT_DIR}/include ${PROJECT_SOURCE_DIR}/include)
target_link_libraries(ins_openmp ${CGNS_DIR}/lib/libcgns.so ${PETSC_DIR}/arch-linux-c-debug/lib/libpetsc.so -L${ARMA_DIR}/lib -larmadillo -L${OPENBLAS_DIR}/lib -lopenblas ${OP2DGTOOLKIT_DIR}/lib/libop2dgtoolkit_openmp.so)
target_compile_definitions(ins_openmp PRIVATE DG_ORDER=4)

op2_application(ins_cuda LIBS op2_cuda op2_hdf5 SOURCES ins_op.cpp ${COMMON_SRC} ${GPU_SRC} ${SN_SRC} ${CUBLAS_SRC} cuda/ins_kernels.cu)
target_compile_definitions(ins_cuda PUBLIC -DBUILD_CUDA)
target_include_directories(ins_cuda PUBLIC ${CGNS_DIR}/include ${ARMA_DIR}/include ${PETSC_DIR}/include ${PETSC_DIR}/arch-linux-c-debug/include ${CUDAToolkit_INCLUDE_DIRS} ${OP2DGTOOLKIT_DIR}/include ${PROJECT_SOURCE_DIR}/include)
target_link_libraries(ins_cuda ${CGNS_DIR}/lib/libcgns.so ${PETSC_DIR}/arch-linux-c-debug/lib/libpetsc.so -L${ARMA_DIR}/lib -larmadillo -L${OPENBLAS_DIR}/lib -lopenblas CUDA::cudart_static CUDA::cublas ${OP2DGTOOLKIT_DIR}/lib/libop2dgtoolkit_cuda.so)
set_property(TARGET ins_cuda PROPERTY CUDA_ARCHITECTURES OFF)
target_compile_definitions(ins_cuda PRIVATE INS_CUDA OP2_DG_CUDA DG_ORDER=4)

op2_application(ins_mpi_openmp LIBS op2_mpi op2_openmp op2_hdf5 SOURCES ins_op.cpp ${COMMON_SRC} ${CPU_SRC} ${MPI_SRC} ${OPENBLAS_SRC} openmp/ins_kernels.cpp)
target_include_directories(ins_mpi_openmp PUBLIC ${CGNS_DIR}/include ${ARMA_DIR}/include ${PETSC_DIR}/include ${PETSC_DIR}/arch-linux-c-debug/include ${OPENBLAS_DIR}/include ${OP2DGTOOLKIT_DIR}/include ${PROJECT_SOURCE_DIR}/include ${PARMETIS_DIR}/include)
target_link_libraries(ins_mpi_openmp ${CGNS_DIR}/lib/libcgns.so ${PETSC_DIR}/arch-linux-c-debug/lib/libpetsc.so ${PARMETIS_DIR}/lib/libparmetis.a -L${ARMA_DIR}/lib -larmadillo -L${OPENBLAS_DIR}/lib -lopenblas ${OP2DGTOOLKIT_DIR}/lib/libop2dgtoolkit_mpi_openmp.so)
target_compile_definitions(ins_mpi_openmp PRIVATE INS_MPI DG_ORDER=4)

op2_application(ins_mpi_cuda LIBS op2_mpi_cuda op2_hdf5 SOURCES ins_op.cpp ${COMMON_SRC} ${GPU_SRC} ${MPI_SRC} ${CUBLAS_SRC} cuda/ins_kernels.cu)
target_compile_definitions(ins_mpi_cuda PUBLIC -DBUILD_CUDA)
target_include_directories(ins_mpi_cuda PUBLIC ${CGNS_DIR}/include ${ARMA_DIR}/include ${PETSC_DIR}/include ${PETSC_DIR}/arch-linux-c-debug/include ${CUDAToolkit_INCLUDE_DIRS} ${OP2DGTOOLKIT_DIR}/include ${PROJECT_SOURCE_DIR}/include)
target_link_libraries(ins_mpi_cuda ${CGNS_DIR}/lib/libcgns.so ${PETSC_DIR}/arch-linux-c-debug/lib/libpetsc.so -L${ARMA_DIR}/lib -larmadillo -L${OPENBLAS_DIR}/lib -lopenblas CUDA::cudart_static CUDA::cublas ${OP2DGTOOLKIT_DIR}/lib/libop2dgtoolkit_mpi_cuda.so)
set_property(TARGET ins_mpi_cuda PROPERTY CUDA_ARCHITECTURES OFF)
target_compile_definitions(ins_mpi_cuda PRIVATE INS_CUDA INS_MPI OP2_DG_CUDA DG_ORDER=4)
