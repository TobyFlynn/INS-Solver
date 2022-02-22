find_package(VTK)

if(NOT VTK_FOUND)
  message(FATAL_ERROR "VTK could not be found")
endif()

add_executable(vtk2CGNS vtk2CGNS.cpp)
target_include_directories(vtk2CGNS PUBLIC ${VTK_USE_FILE})
target_link_libraries(vtk2CGNS libcgns.so libhdf5.so ${VTK_LIBRARIES})

add_executable(interpolateFromCGNS interpolateFromCGNS.cpp)
target_include_directories(interpolateFromCGNS PUBLIC ${CGNS_DIR}/include)
target_link_libraries(interpolateFromCGNS libcgns.so libhdf5.so ${OP2DGTOOLKIT_DIR}/lib/libdgtoolkit.so -L${ARMA_DIR}/lib -larmadillo ${OpenBLAS_LIBRARIES})
