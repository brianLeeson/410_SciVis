cmake_minimum_required(VERSION 3.10)

PROJECT(proj5)
SET(VTK_DIR ~/Builds/VTK-build/)

find_package(VTK REQUIRED)
include(${VTK_USE_FILE})

add_executable(proj5 proj5.cxx)

#SET(CMAKE_CXX_LINK_FLAGS "-framework OpenGL")
if(VTK_LIBRARIES)
  target_link_libraries(proj5 ${VTK_LIBRARIES})
else()
  target_link_libraries(proj5 vtkHybrid)
endif()


