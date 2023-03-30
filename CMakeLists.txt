#---------------------------------------------------------------------------------------------------
# CMake Build Configuration for the Ray Tracing Weekend Series
#
# See README.md for guidance.
#---------------------------------------------------------------------------------------------------

cmake_minimum_required ( VERSION 3.1.0 )

project ( MCPT_LY
  VERSION 3.0.0
  LANGUAGES CXX
)


# Set to c++11
set ( CMAKE_CXX_STANDARD 11 )

# Source
set ( COMMON_ALL
  src/common/rtweekend.h
  src/common/camera.h
  src/common/ray.h
  src/common/vec3.h
  main.cc
)

file(GLOB SOURCES ${CMAKE_SOURCE_DIR}/src/*.cc ${CMAKE_SOURCE_DIR}/src/*.h ${CMAKE_SOURCE_DIR}/src/*.cpp)

set ( SOURCES
  ${COMMON_ALL}
  ${SOURCES}
)

add_library(tinyobj STATIC tiny_obj_loader.cc)

# Executables
add_executable(MCPT_LY ${SOURCES})

include_directories(src/common src)
target_link_libraries(MCPT_LY tinyobj)
find_package(OpenMP)
if(OpenMP_CXX_FOUND)
    target_link_libraries(MCPT_LY OpenMP::OpenMP_CXX)
endif()