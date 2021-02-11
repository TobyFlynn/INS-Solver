find_package(VTK)

if(NOT VTK_FOUND)
  message(FATAL_ERROR "VTK could not be found")
endif()

add_executable(vtk2CGNS vtk2CGNS.cpp)

# Link CGNS Library
target_include_directories(vtk2CGNS PUBLIC ${CGNS_DIR}/include ${VTK_USE_FILE})
target_link_libraries(vtk2CGNS ${CGNS_DIR}/lib/libcgns.so ${VTK_LIBRARIES})
