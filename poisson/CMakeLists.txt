set(KERNEL_HEADERS kernels/init_grid.h kernels/set_ic1.h kernels/set_ic2.h kernels/set_tau.h kernels/set_tau_bc.h kernels/set_rhs.h kernels/calc_sol.h)
list(APPEND KERNEL_HEADERS kernels/div.h kernels/curl.h kernels/grad.h)
list(APPEND KERNEL_HEADERS kernels/pRHS_faces.h kernels/pRHS_bc.h kernels/pRHS_du.h kernels/pRHS_fluxq.h kernels/pRHS_J.h kernels/pRHS_qbc.h)
set(OPENBLAS_SRC openBLAS_op/init_grid_op.cpp openBLAS_op/div_op.cpp openBLAS_op/grad_op.cpp openBLAS_op/poisson_set_rhs_op.cpp openBLAS_op/poisson_rhs1_op.cpp openBLAS_op/poisson_rhs2_op.cpp openBLAS_op/poisson_calc_sol_op.cpp)
set(CUBLAS_SRC cuBLAS_op/init_grid_op.cu cuBLAS_op/div_op.cu cuBLAS_op/grad_op.cu cuBLAS_op/poisson_set_rhs_op.cu cuBLAS_op/poisson_rhs1_op.cu cuBLAS_op/poisson_rhs2_op.cu cuBLAS_op/poisson_calc_sol_op.cu)

op2_application(poisson_seq LIBS op2_seq op2_hdf5 SOURCES poisson_op.cpp operators_op.cpp poisson_rhs_op.cpp poisson_rhs_solve.cpp ins_data.cpp load_mesh.cpp save_solution.cpp ${OPENBLAS_SRC} seq/poisson_seqkernels.cpp ${KERNEL_HEADERS})

# Link CGNS Library
target_include_directories(poisson_seq PUBLIC ${CGNS_DIR}/include ${PETSC_DIR}/include ${PETSC_DIR}/arch-linux-c-debug/include ${OPENBLAS_DIR}/include)
target_link_libraries(poisson_seq ${CGNS_DIR}/lib/libcgns.so ${PETSC_DIR}/arch-linux-c-debug/lib/libpetsc.so -L${OPENBLAS_DIR}/lib -lopenblas)

op2_application(poisson_openmp LIBS op2_openmp op2_hdf5 SOURCES poisson_op.cpp operators_op.cpp poisson_rhs_op.cpp poisson_rhs_solve.cpp ins_data.cpp load_mesh.cpp save_solution.cpp ${OPENBLAS_SRC} openmp/poisson_kernels.cpp ${KERNEL_HEADERS})

# Link CGNS Library
target_include_directories(poisson_openmp PUBLIC ${CGNS_DIR}/include ${PETSC_DIR}/include ${PETSC_DIR}/arch-linux-c-debug/include ${OPENBLAS_DIR}/include)
target_link_libraries(poisson_openmp ${CGNS_DIR}/lib/libcgns.so ${PETSC_DIR}/arch-linux-c-debug/lib/libpetsc.so -L${OPENBLAS_DIR}/lib -lopenblas)

op2_application(poisson_cuda LIBS op2_cuda op2_hdf5 SOURCES poisson_op.cpp operators_op.cpp poisson_rhs_op.cu poisson_rhs_solve.cu ins_data.cpp load_mesh.cpp save_solution.cpp ${CUBLAS_SRC} cuda/poisson_kernels.cu ${KERNEL_HEADERS})
target_compile_definitions(poisson_cuda PUBLIC -DBUILD_CUDA)
# Link CGNS Library
target_include_directories(poisson_cuda PUBLIC ${CGNS_DIR}/include ${PETSC_DIR}/include ${PETSC_DIR}/arch-linux-c-debug/include ${CUDAToolkit_INCLUDE_DIRS})
target_link_libraries(poisson_cuda ${CGNS_DIR}/lib/libcgns.so ${PETSC_DIR}/arch-linux-c-debug/lib/libpetsc.so CUDA::cudart_static CUDA::cublas)
set_property(TARGET poisson_cuda PROPERTY CUDA_ARCHITECTURES OFF)
