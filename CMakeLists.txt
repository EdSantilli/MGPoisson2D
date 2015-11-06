cmake_minimum_required (VERSION 2.8)
enable_language (Fortran)

project (MGPoisson_project)
set (MGPoisson_project_Tutorial_VERSION_MAJOR 0)
set (MGPoisson_project_Tutorial_VERSION_MINOR 1)

# # configure a header file to pass some of the CMake settings
# # to the source code
# configure_file (
#   "${PROJECT_SOURCE_DIR}/TutorialConfig.h.in"
#   "${PROJECT_BINARY_DIR}/TutorialConfig.h"
#   )

# set (CMAKE_VERBOSE_MAKEFILE "On")

add_subdirectory (libMGPoisson)
add_subdirectory (tests)
