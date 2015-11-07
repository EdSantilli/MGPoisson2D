# CMake basics.
cmake_minimum_required (VERSION 2.8)
enable_language (Fortran)

# Project name and version.
project (MGPoisson_project)
set (MGPoisson_project_Tutorial_VERSION_MAJOR 0)
set (MGPoisson_project_Tutorial_VERSION_MINOR 1)

# make verbosity.
# set (CMAKE_VERBOSE_MAKEFILE "On")

# The source folders that will be "make"ed.
add_subdirectory (libBaseTools)
add_subdirectory (libMGPoisson)
add_subdirectory (tests)
