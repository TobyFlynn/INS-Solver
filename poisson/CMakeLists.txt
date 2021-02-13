set(KERNEL_HEADERS kernels/init_grid.h kernels/set_ic1.h kernels/set_ic2.h)
list(APPEND KERNEL_HEADERS kernels/div.h kernels/curl.h kernels/grad.h)
set(OPENBLAS_SRC openBLAS_op/init_grid_op.cpp openBLAS_op/div_op.cpp openBLAS_op/grad_op.cpp)
set(CUBLAS_SRC cuBLAS_op/init_grid_op.cu cuBLAS_op/div_op.cu cuBLAS_op/grad_op.cu)

op2_application(poisson_seq LIBS op2_seq op2_hdf5 SOURCES poisson_op.cpp ins_data.cpp load_mesh.cpp save_solution.cpp ${OPENBLAS_SRC} seq/poisson_seqkernels.cpp ${KERNEL_HEADERS})

# Link CGNS Library
target_include_directories(poisson_seq PUBLIC ${CGNS_DIR}/include ${OPENBLAS_DIR}/include)
target_link_libraries(poisson_seq ${CGNS_DIR}/lib/libcgns.so -L${OPENBLAS_DIR}/lib -lopenblas)

op2_application(poisson_openmp LIBS op2_openmp op2_hdf5 SOURCES poisson_op.cpp ins_data.cpp load_mesh.cpp save_solution.cpp ${OPENBLAS_SRC} openmp/poisson_kernels.cpp ${KERNEL_HEADERS})

# Link CGNS Library
target_include_directories(poisson_openmp PUBLIC ${CGNS_DIR}/include ${OPENBLAS_DIR}/include)
target_link_libraries(poisson_openmp ${CGNS_DIR}/lib/libcgns.so -L${OPENBLAS_DIR}/lib -lopenblas)

op2_application(poisson_cuda LIBS op2_cuda op2_hdf5 SOURCES poisson_op.cpp ins_data.cpp load_mesh.cpp save_solution.cpp ${CUBLAS_SRC} cuda/poisson_kernels.cu ${KERNEL_HEADERS})
target_compile_definitions(poisson_cuda PUBLIC -DBUILD_CUDA)
# Link CGNS Library
target_include_directories(poisson_cuda PUBLIC ${CGNS_DIR}/include ${CUDAToolkit_INCLUDE_DIRS})
target_link_libraries(poisson_cuda ${CGNS_DIR}/lib/libcgns.so CUDA::cudart_static CUDA::cublas)
set_property(TARGET poisson_cuda PROPERTY CUDA_ARCHITECTURES OFF)
