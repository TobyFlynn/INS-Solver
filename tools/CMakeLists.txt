find_package(VTK)

if(NOT VTK_FOUND)
  message(FATAL_ERROR "VTK could not be found")
endif()

add_executable(vtk2CGNS vtk2CGNS.cpp)
target_include_directories(vtk2CGNS PUBLIC ${CGNS_DIR}/include ${VTK_USE_FILE})
target_link_libraries(vtk2CGNS ${CGNS_DIR}/lib/libcgns.so ${VTK_LIBRARIES})

add_executable(matlab2CGNS matlab2CGNS.cpp)
target_include_directories(matlab2CGNS PUBLIC ${CGNS_DIR}/include)
target_link_libraries(matlab2CGNS ${CGNS_DIR}/lib/libcgns.so)

add_executable(vortexError vortexError.cpp)
target_include_directories(vortexError PUBLIC ${CGNS_DIR}/include)
target_link_libraries(vortexError ${CGNS_DIR}/lib/libcgns.so)

add_executable(diff diff.cpp)
target_include_directories(diff PUBLIC ${CGNS_DIR}/include)
target_link_libraries(diff ${CGNS_DIR}/lib/libcgns.so)
