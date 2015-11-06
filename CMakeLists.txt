cmake_minimum_required (VERSION 2.8)
enable_language (Fortran)
project (MGPoisson_project)

# set (CMAKE_VERBOSE_MAKEFILE "On")

add_subdirectory (libMGPoisson)
add_subdirectory (tests)
