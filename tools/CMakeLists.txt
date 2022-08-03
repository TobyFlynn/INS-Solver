find_package(VTK)

if(NOT VTK_FOUND)
  message(FATAL_ERROR "VTK could not be found")
endif()

add_executable(vtk2CGNS vtk2CGNS.cpp)
target_include_directories(vtk2CGNS PUBLIC ${VTK_USE_FILE})
target_link_libraries(vtk2CGNS ${CGNS_LIB} ${VTK_LIBRARIES})

#add_executable(matlab2CGNS matlab2CGNS.cpp)
#target_link_libraries(matlab2CGNS ${CGNS_LIB})

#add_executable(vortexError vortexError.cpp)
#target_link_libraries(vortexError ${CGNS_LIB})

#add_executable(interpolateFromCGNS interpolateFromCGNS.cpp)
#target_link_libraries(interpolateFromCGNS ${CGNS_LIB} ${OP2DGTOOLKIT_DIR}/lib/libdgtoolkit.so -L${ARMA_DIR}/lib -larmadillo -L${OPENBLAS_DIR}/lib -lopenblas)
