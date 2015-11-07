cmake_minimum_required (VERSION 2.8)
enable_language (Fortran)

project (MGPoisson_project)
set (MGPoisson_project_Tutorial_VERSION_MAJOR 0)
set (MGPoisson_project_Tutorial_VERSION_MINOR 1)

# set (CMAKE_VERBOSE_MAKEFILE "On")

include_directories(${CMAKE_BINARY_DIR}/include)

add_subdirectory (libBaseTools)
add_subdirectory (libMGPoisson)
add_subdirectory (tests)
