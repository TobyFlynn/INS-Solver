set(KERNEL_HEADERS kernels/init_grid.h kernels/set_ic.h kernels/calc_dt.h kernels/lift_drag.h kernels/init_nodes.h)
list(APPEND KERNEL_HEADERS kernels/init_cubature.h kernels/init_cubature_grad.h kernels/init_cubature_OP.h)
list(APPEND KERNEL_HEADERS kernels/gauss_tau.h kernels/gauss_tau_bc.h kernels/init_gauss.h kernels/init_gauss_grad.h kernels/init_gauss_grad_neighbour.h kernels/gauss_grad_faces.h kernels/init_gauss_grad2.h kernels/gauss_op.h kernels/gauss_gfi_faces.h kernels/gauss_reverse.h)
list(APPEND KERNEL_HEADERS kernels/div.h kernels/curl.h kernels/grad.h)
list(APPEND KERNEL_HEADERS kernels/advection_flux.h kernels/advection_faces.h kernels/advection_bc.h kernels/advection_numerical_flux.h kernels/advection_intermediate_vel.h)
list(APPEND KERNEL_HEADERS kernels/pressure_bc.h kernels/pressure_bc2.h kernels/pressure_rhs.h kernels/pressure_update_vel.h)
list(APPEND KERNEL_HEADERS kernels/viscosity_rhs.h kernels/viscosity_reset_bc.h kernels/viscosity_bc.h)
list(APPEND KERNEL_HEADERS kernels/poisson_test_init.h kernels/poisson_test_bc.h kernels/poisson_test_error.h)
list(APPEND KERNEL_HEADERS kernels/poisson_mf2.h kernels/poisson_mf2_mass.h kernels/poisson_mf2_faces.h kernels/poisson_mf2_op.h kernels/poisson_mf2_opf.h  kernels/poisson_mf2_opbf.h kernels/poisson_mf2_bc.h kernels/poisson_mf2_apply_bc)

set(OPENBLAS_SRC openBLAS_op/init_grid_op.cpp openBLAS_op/op2_gemv_op.cpp openBLAS_op/op2_gemm_op.cpp openBLAS_op/op2_gemv_batch_op.cpp openBLAS_op/op2_gemm_batch_op.cpp)
list(APPEND OPENBLAS_SRC openBLAS_op/init_gauss_op.cpp openBLAS_op/init_gauss_grad_op.cpp openBLAS_op/init_gauss_grad_neighbour_op.cpp)
#list(APPEND OPENBLAS_SRC openBLAS_op/poisson_mf2_op.cpp)

set(CUBLAS_SRC cuBLAS_op/init_grid_op.cu cuBLAS_op/op2_gemv_op.cu cuBLAS_op/op2_gemm_op.cu cuBLAS_op/op2_gemv_batch_op.cu cuBLAS_op/op2_gemm_batch_op.cu)
list(APPEND CUBLAS_SRC cuBLAS_op/init_gauss_op.cu cuBLAS_op/init_gauss_grad_op.cu cuBLAS_op/init_gauss_grad_neighbour_op.cu)
#list(APPEND CUBLAS_SRC cuBLAS_op/poisson_mf2_op.cu)

set(COMMON_SRC solver_op.cpp operators_op.cpp poisson_op.cpp poisson_m_op.cpp poisson_mf2_op.cpp ins_data_op.cpp timing_op.cpp)
set(CPU_SRC poisson_copy_op.cpp constants.cpp)
set(GPU_SRC poisson_copy_op.cu constants.cu)
set(SN_SRC save_solution_op.cpp load_mesh.cpp)
set(MPI_SRC save_solution_mpi_op.cpp load_mesh_mpi.cpp mpi_helper_func.cpp)

op2_application(ins_seq LIBS op2_seq op2_hdf5 SOURCES ins_op.cpp ${COMMON_SRC} ${CPU_SRC} ${SN_SRC} ${OPENBLAS_SRC} seq/ins_seqkernels.cpp ${KERNEL_HEADERS})

# Link CGNS Library
target_include_directories(ins_seq PUBLIC ${CGNS_DIR}/include ${PETSC_DIR}/include ${PETSC_DIR}/arch-linux-c-debug/include ${OPENBLAS_DIR}/include ${PROJECT_SOURCE_DIR}/include)
target_link_libraries(ins_seq ${CGNS_DIR}/lib/libcgns.so ${PETSC_DIR}/arch-linux-c-debug/lib/libpetsc.so -L${OPENBLAS_DIR}/lib -lopenblas)

op2_application(ins_openmp LIBS op2_openmp op2_hdf5 SOURCES ins_op.cpp ${COMMON_SRC} ${CPU_SRC} ${SN_SRC} ${OPENBLAS_SRC} openmp/ins_kernels.cpp ${KERNEL_HEADERS})

# Link CGNS Library
target_include_directories(ins_openmp PUBLIC ${CGNS_DIR}/include ${PETSC_DIR}/include ${PETSC_DIR}/arch-linux-c-debug/include ${OPENBLAS_DIR}/include ${PROJECT_SOURCE_DIR}/include)
target_link_libraries(ins_openmp ${CGNS_DIR}/lib/libcgns.so ${PETSC_DIR}/arch-linux-c-debug/lib/libpetsc.so -L${OPENBLAS_DIR}/lib -lopenblas)

op2_application(ins_cuda LIBS op2_cuda op2_hdf5 SOURCES ins_op.cpp ${COMMON_SRC} ${GPU_SRC} ${SN_SRC} ${CUBLAS_SRC} cuda/ins_kernels.cu ${KERNEL_HEADERS})
target_compile_definitions(ins_cuda PUBLIC -DBUILD_CUDA)
# Link CGNS Library
target_include_directories(ins_cuda PUBLIC ${CGNS_DIR}/include ${PETSC_DIR}/include ${PETSC_DIR}/arch-linux-c-debug/include ${CUDAToolkit_INCLUDE_DIRS} ${PROJECT_SOURCE_DIR}/include)
target_link_libraries(ins_cuda ${CGNS_DIR}/lib/libcgns.so ${PETSC_DIR}/arch-linux-c-debug/lib/libpetsc.so CUDA::cudart_static CUDA::cublas)
set_property(TARGET ins_cuda PROPERTY CUDA_ARCHITECTURES OFF)
target_compile_definitions(ins_cuda PRIVATE INS_CUDA)

op2_application(ins_mpi LIBS op2_mpi op2_hdf5 SOURCES ins_op.cpp ${COMMON_SRC} ${CPU_SRC} ${MPI_SRC} ${OPENBLAS_SRC} seq/ins_seqkernels.cpp ${KERNEL_HEADERS})
target_include_directories(ins_mpi PUBLIC ${CGNS_DIR}/include ${PETSC_DIR}/include ${PETSC_DIR}/arch-linux-c-debug/include ${OPENBLAS_DIR}/include ${PROJECT_SOURCE_DIR}/include ${PARMETIS_DIR}/include)
target_link_libraries(ins_mpi ${CGNS_DIR}/lib/libcgns.so ${PETSC_DIR}/arch-linux-c-debug/lib/libpetsc.so ${PARMETIS_DIR}/lib/libparmetis.a -L${OPENBLAS_DIR}/lib -lopenblas)
target_compile_definitions(ins_mpi PRIVATE INS_MPI)

op2_application(ins_mpi_openmp LIBS op2_mpi op2_openmp op2_hdf5 SOURCES ins_op.cpp ${COMMON_SRC} ${CPU_SRC} ${MPI_SRC} ${OPENBLAS_SRC} openmp/ins_kernels.cpp ${KERNEL_HEADERS})
target_include_directories(ins_mpi_openmp PUBLIC ${CGNS_DIR}/include ${PETSC_DIR}/include ${PETSC_DIR}/arch-linux-c-debug/include ${OPENBLAS_DIR}/include ${PROJECT_SOURCE_DIR}/include ${PARMETIS_DIR}/include)
target_link_libraries(ins_mpi_openmp ${CGNS_DIR}/lib/libcgns.so ${PETSC_DIR}/arch-linux-c-debug/lib/libpetsc.so ${PARMETIS_DIR}/lib/libparmetis.a -L${OPENBLAS_DIR}/lib -lopenblas)
target_compile_definitions(ins_mpi_openmp PRIVATE INS_MPI)

op2_application(ins_mpi_cuda LIBS op2_mpi_cuda op2_hdf5 SOURCES ins_op.cpp ${COMMON_SRC} ${GPU_SRC} ${MPI_SRC} ${CUBLAS_SRC} cuda/ins_kernels.cu ${KERNEL_HEADERS})
target_compile_definitions(ins_mpi_cuda PUBLIC -DBUILD_CUDA)
# Link CGNS Library
target_include_directories(ins_mpi_cuda PUBLIC ${CGNS_DIR}/include ${PETSC_DIR}/include ${PETSC_DIR}/arch-linux-c-debug/include ${CUDAToolkit_INCLUDE_DIRS} ${PROJECT_SOURCE_DIR}/include)
target_link_libraries(ins_mpi_cuda ${CGNS_DIR}/lib/libcgns.so ${PETSC_DIR}/arch-linux-c-debug/lib/libpetsc.so CUDA::cudart_static CUDA::cublas)
set_property(TARGET ins_mpi_cuda PROPERTY CUDA_ARCHITECTURES OFF)
target_compile_definitions(ins_mpi_cuda PRIVATE INS_CUDA INS_MPI)

op2_application(p_openmp LIBS op2_openmp op2_hdf5 SOURCES poisson_test_op.cpp ${COMMON_SRC} ${CPU_SRC} ${SN_SRC} ${OPENBLAS_SRC} openmp/ins_kernels.cpp ${KERNEL_HEADERS})
target_include_directories(p_openmp PUBLIC ${CGNS_DIR}/include ${PETSC_DIR}/include ${PETSC_DIR}/arch-linux-c-debug/include ${OPENBLAS_DIR}/include ${PROJECT_SOURCE_DIR}/include)
target_link_libraries(p_openmp ${CGNS_DIR}/lib/libcgns.so ${PETSC_DIR}/arch-linux-c-debug/lib/libpetsc.so -L${OPENBLAS_DIR}/lib -lopenblas)
target_compile_definitions(p_openmp PRIVATE POISSON_TEST)

op2_application(p_cuda LIBS op2_cuda op2_hdf5 SOURCES poisson_test_op.cpp ${COMMON_SRC} ${GPU_SRC} ${SN_SRC} ${CUBLAS_SRC} cuda/ins_kernels.cu ${KERNEL_HEADERS})
target_compile_definitions(p_cuda PUBLIC -DBUILD_CUDA)
target_include_directories(p_cuda PUBLIC ${CGNS_DIR}/include ${PETSC_DIR}/include ${PETSC_DIR}/arch-linux-c-debug/include ${PROJECT_SOURCE_DIR}/include ${CUDAToolkit_INCLUDE_DIRS})
target_link_libraries(p_cuda ${CGNS_DIR}/lib/libcgns.so ${PETSC_DIR}/arch-linux-c-debug/lib/libpetsc.so CUDA::cudart_static CUDA::cublas)
set_property(TARGET p_cuda PROPERTY CUDA_ARCHITECTURES OFF)
target_compile_definitions(p_cuda PRIVATE INS_CUDA POISSON_TEST)

op2_application(p_mpi_openmp LIBS op2_mpi op2_openmp op2_hdf5 SOURCES poisson_test_op.cpp ${COMMON_SRC} ${CPU_SRC} ${MPI_SRC} ${OPENBLAS_SRC} openmp/ins_kernels.cpp ${KERNEL_HEADERS})
target_include_directories(p_mpi_openmp PUBLIC ${CGNS_DIR}/include ${PETSC_DIR}/include ${PETSC_DIR}/arch-linux-c-debug/include ${OPENBLAS_DIR}/include ${PROJECT_SOURCE_DIR}/include ${PARMETIS_DIR}/include)
target_link_libraries(p_mpi_openmp ${CGNS_DIR}/lib/libcgns.so ${PETSC_DIR}/arch-linux-c-debug/lib/libpetsc.so -L${OPENBLAS_DIR}/lib -lopenblas)
target_compile_definitions(p_mpi_openmp PRIVATE INS_MPI POISSON_TEST)

op2_application(p_mpi_cuda LIBS op2_mpi_cuda op2_hdf5 SOURCES poisson_test_op.cpp ${COMMON_SRC} ${GPU_SRC} ${MPI_SRC} ${CUBLAS_SRC} cuda/ins_kernels.cu ${KERNEL_HEADERS})
target_compile_definitions(p_mpi_cuda PUBLIC -DBUILD_CUDA)
target_include_directories(p_mpi_cuda PUBLIC ${CGNS_DIR}/include ${PETSC_DIR}/include ${PETSC_DIR}/arch-linux-c-debug/include ${PROJECT_SOURCE_DIR}/include ${PARMETIS_DIR}/include ${CUDAToolkit_INCLUDE_DIRS})
target_link_libraries(p_mpi_cuda ${CGNS_DIR}/lib/libcgns.so ${PETSC_DIR}/arch-linux-c-debug/lib/libpetsc.so CUDA::cudart_static CUDA::cublas)
set_property(TARGET p_mpi_cuda PROPERTY CUDA_ARCHITECTURES OFF)
target_compile_definitions(p_mpi_cuda PRIVATE INS_MPI INS_CUDA POISSON_TEST)
