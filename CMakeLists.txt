cmake_minimum_required( VERSION 2.6 )
project( CircleDetection )

set( CMAKE_CXX_FLAGS_DEBUG OFF )
set( CMAKE_CXX_FLAGS_RELEASE ON )

include_directories( ${PROJECT_SOURCE_DIR}/circle_detection )
add_subdirectory( circle_detection )

add_executable( ${PROJECT_NAME} main.cpp )

target_link_libraries( ${PROJECT_NAME} circdet )

enable_testing()

add_test( RunTest CircleDetection )