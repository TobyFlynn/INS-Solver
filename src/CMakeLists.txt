set(KERNEL_HEADERS kernels/init_grid.h kernels/set_ic.h kernels/calc_dt.h kernels/lift_drag.h kernels/min_max.h)
list(APPEND KERNEL_HEADERS kernels/init_cubature.h kernels/init_cubature_grad.h kernels/init_cubature_OP.h)
list(APPEND KERNEL_HEADERS kernels/gauss_tau.h kernels/gauss_tau_bc.h kernels/init_gauss.h kernels/init_gauss_grad.h kernels/init_gauss_grad_neighbour.h kernels/gauss_grad_faces.h kernels/init_gauss_grad2.h kernels/gauss_op.h kernels/gauss_gfi_faces.h kernels/gauss_reverse.h)
list(APPEND KERNEL_HEADERS kernels/div.h kernels/curl.h kernels/grad.h)
list(APPEND KERNEL_HEADERS kernels/advection_flux.h kernels/advection_faces.h kernels/advection_bc.h kernels/advection_numerical_flux.h kernels/advection_intermediate_vel.h)
list(APPEND KERNEL_HEADERS kernels/pressure_bc.h kernels/pressure_rhs.h kernels/pressure_update_vel.h)
list(APPEND KERNEL_HEADERS kernels/viscosity_rhs.h kernels/viscosity_reset_bc.h kernels/viscosity_bc.h)
list(APPEND KERNEL_HEADERS kernels/poisson_rhs_faces.h kernels/poisson_rhs_bc.h kernels/poisson_rhs_flux.h kernels/poisson_rhs_J.h kernels/poisson_rhs_qbc.h kernels/poisson_rhs_qflux.h kernels/poisson_bc.h kernels/poisson_bc_J.h kernels/poisson_bc2.h kernels/poisson_bc3.h)
list(APPEND KERNEL_HEADERS kernels/cub_grad.h kernels/cub_div.h kernels/poisson_test_init.h kernels/poisson_test_bc.h kernels/poisson_test_error.h kernels/tau.h kernels/tau_bc.h)
list(APPEND KERNEL_HEADERS kernels/poisson_mf2.h kernels/poisson_mf2_mass.h kernels/poisson_mf2_faces.h kernels/poisson_mf2_op.h kernels/poisson_mf2_opf.h  kernels/poisson_mf2_opbf.h kernels/poisson_mf2_bc.h kernels/poisson_mf2_apply_bc)

set(OPENBLAS_SRC openBLAS_op/init_grid_op.cpp openBLAS_op/op2_gemv_op.cpp openBLAS_op/viscosity_rhs_op.cpp)
list(APPEND OPENBLAS_SRC openBLAS_op/cubature_mm_op.cpp openBLAS_op/cubature_op_op.cpp)
list(APPEND OPENBLAS_SRC openBLAS_op/init_gauss_op.cpp openBLAS_op/init_gauss_grad_op.cpp openBLAS_op/gauss_op_op.cpp openBLAS_op/gauss_opf_op.cpp openBLAS_op/init_gauss_grad_neighbour_op.cpp)
list(APPEND OPENBLAS_SRC openBLAS_op/poisson_rhs_mass_op.cpp openBLAS_op/poisson_mf2_op.cpp)

set(CUBLAS_SRC cuBLAS_op/init_grid_op.cu cuBLAS_op/op2_gemv_op.cu cuBLAS_op/viscosity_rhs_op.cu)
list(APPEND CUBLAS_SRC cuBLAS_op/cubature_mm_op.cu cuBLAS_op/cubature_op_op.cu)
list(APPEND CUBLAS_SRC cuBLAS_op/init_gauss_op.cu cuBLAS_op/init_gauss_grad_op.cu cuBLAS_op/gauss_op_op.cu cuBLAS_op/gauss_opf_op.cu cuBLAS_op/init_gauss_grad_neighbour_op.cu)
list(APPEND CUBLAS_SRC cuBLAS_op/poisson_rhs_mass_op.cu cuBLAS_op/poisson_mf2_op.cu)

set(CPU_SRC operators_op.cpp poisson_op.cpp poisson_m_op.cpp poisson_mf_op.cpp poisson_mf2_op.cpp poisson_copy_op.cpp ins_data_op.cpp save_solution_op.cpp timing_op.cpp constants.cpp)
set(GPU_SRC operators_op.cpp poisson_op.cpp poisson_m_op.cpp poisson_mf_op.cpp poisson_mf2_op.cpp poisson_copy_op.cu ins_data_op.cpp save_solution_op.cpp timing_op.cpp constants.cu)

op2_application(ins_seq LIBS op2_seq op2_hdf5 SOURCES ins_op.cpp ${CPU_SRC} ${OPENBLAS_SRC} seq/ins_seqkernels.cpp ${KERNEL_HEADERS})

# Link CGNS Library
target_include_directories(ins_seq PUBLIC ${CGNS_DIR}/include ${PETSC_DIR}/include ${PETSC_DIR}/arch-linux-c-debug/include ${OPENBLAS_DIR}/include)
target_link_libraries(ins_seq ${CGNS_DIR}/lib/libcgns.so ${PETSC_DIR}/arch-linux-c-debug/lib/libpetsc.so -L${OPENBLAS_DIR}/lib -lopenblas)

op2_application(ins_openmp LIBS op2_openmp op2_hdf5 SOURCES ins_op.cpp ${CPU_SRC} ${OPENBLAS_SRC} openmp/ins_kernels.cpp ${KERNEL_HEADERS})

# Link CGNS Library
target_include_directories(ins_openmp PUBLIC ${CGNS_DIR}/include ${PETSC_DIR}/include ${PETSC_DIR}/arch-linux-c-debug/include ${OPENBLAS_DIR}/include)
target_link_libraries(ins_openmp ${CGNS_DIR}/lib/libcgns.so ${PETSC_DIR}/arch-linux-c-debug/lib/libpetsc.so -L${OPENBLAS_DIR}/lib -lopenblas)

op2_application(ins_cuda LIBS op2_cuda op2_hdf5 SOURCES ins_op.cpp ${GPU_SRC} ${CUBLAS_SRC} cuda/ins_kernels.cu ${KERNEL_HEADERS})
target_compile_definitions(ins_cuda PUBLIC -DBUILD_CUDA)
# Link CGNS Library
target_include_directories(ins_cuda PUBLIC ${CGNS_DIR}/include ${PETSC_DIR}/include ${PETSC_DIR}/arch-linux-c-debug/include ${CUDAToolkit_INCLUDE_DIRS})
target_link_libraries(ins_cuda ${CGNS_DIR}/lib/libcgns.so ${PETSC_DIR}/arch-linux-c-debug/lib/libpetsc.so CUDA::cudart_static CUDA::cublas)
set_property(TARGET ins_cuda PROPERTY CUDA_ARCHITECTURES OFF)
target_compile_definitions(ins_cuda PRIVATE INS_CUDA)

op2_application(p_openmp LIBS op2_openmp op2_hdf5 SOURCES poisson_test_op.cpp ${CPU_SRC} ${OPENBLAS_SRC} openmp/ins_kernels.cpp ${KERNEL_HEADERS})
target_include_directories(p_openmp PUBLIC ${CGNS_DIR}/include ${PETSC_DIR}/include ${PETSC_DIR}/arch-linux-c-debug/include ${OPENBLAS_DIR}/include)
target_link_libraries(p_openmp ${CGNS_DIR}/lib/libcgns.so ${PETSC_DIR}/arch-linux-c-debug/lib/libpetsc.so -L${OPENBLAS_DIR}/lib -lopenblas)

op2_application(p_cuda LIBS op2_cuda op2_hdf5 SOURCES poisson_test_op.cpp ${GPU_SRC} ${CUBLAS_SRC} cuda/ins_kernels.cu ${KERNEL_HEADERS})
target_compile_definitions(p_cuda PUBLIC -DBUILD_CUDA)
target_include_directories(p_cuda PUBLIC ${CGNS_DIR}/include ${PETSC_DIR}/include ${PETSC_DIR}/arch-linux-c-debug/include ${CUDAToolkit_INCLUDE_DIRS})
target_link_libraries(p_cuda ${CGNS_DIR}/lib/libcgns.so ${PETSC_DIR}/arch-linux-c-debug/lib/libpetsc.so CUDA::cudart_static CUDA::cublas)
set_property(TARGET p_cuda PROPERTY CUDA_ARCHITECTURES OFF)
target_compile_definitions(p_cuda PRIVATE INS_CUDA)
