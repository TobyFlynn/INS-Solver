set(KERNEL_HEADERS kernels/init_grid.h kernels/set_ic.h)
list(APPEND KERNEL_HEADERS kernels/div.h kernels/curl.h kernels/grad.h)
list(APPEND KERNEL_HEADERS kernels/advection_flux.h kernels/advection_faces.h kernels/advection_bc.h kernels/advection_numerical_flux.h kernels/advection_intermediate_vel.h)
list(APPEND KERNEL_HEADERS kernels/pressure_bc.h kernels/pressure_rhs.h)
list(APPEND KERNEL_HEADERS kernels/setup_poisson.h kernels/set_tau.h kernels/set_tau_bc.h)
set(OPENBLAS_SRC openBLAS_op/init_grid_op.cpp openBLAS_op/div_op.cpp openBLAS_op/grad_op.cpp openBLAS_op/advection_lift_op.cpp openBLAS_op/pressure_rhs_op.cpp)
set(CUBLAS_SRC cuBLAS_op/init_grid_op.cu cuBLAS_op/div_op.cu cuBLAS_op/grad_op.cu cuBLAS_op/advection_lift_op.cu cuBLAS_op/pressure_rhs_op.cu)

op2_application(ins_seq LIBS op2_seq op2_hdf5 SOURCES ins_op.cpp poisson_op.cpp ins_data.cpp load_mesh.cpp save_solution.cpp ${OPENBLAS_SRC} seq/ins_seqkernels.cpp ${KERNEL_HEADERS})

# Link CGNS Library
target_include_directories(ins_seq PUBLIC ${CGNS_DIR}/include ${OPENBLAS_DIR}/include)
target_link_libraries(ins_seq ${CGNS_DIR}/lib/libcgns.so -L${OPENBLAS_DIR}/lib -lopenblas)

op2_application(ins_openmp LIBS op2_openmp op2_hdf5 SOURCES ins_op.cpp poisson_op.cpp ins_data.cpp load_mesh.cpp save_solution.cpp ${OPENBLAS_SRC} openmp/ins_kernels.cpp ${KERNEL_HEADERS})

# Link CGNS Library
target_include_directories(ins_openmp PUBLIC ${CGNS_DIR}/include ${OPENBLAS_DIR}/include)
target_link_libraries(ins_openmp ${CGNS_DIR}/lib/libcgns.so -L${OPENBLAS_DIR}/lib -lopenblas)

op2_application(ins_cuda LIBS op2_cuda op2_hdf5 SOURCES ins_op.cpp poisson_op.cpp ins_data.cpp load_mesh.cpp save_solution.cpp ${CUBLAS_SRC} cuda/ins_kernels.cu ${KERNEL_HEADERS})
target_compile_definitions(ins_cuda PUBLIC -DBUILD_CUDA)
# Link CGNS Library
target_include_directories(ins_cuda PUBLIC ${CGNS_DIR}/include ${CUDAToolkit_INCLUDE_DIRS})
target_link_libraries(ins_cuda ${CGNS_DIR}/lib/libcgns.so CUDA::cudart_static CUDA::cublas)
set_property(TARGET ins_cuda PROPERTY CUDA_ARCHITECTURES OFF)
