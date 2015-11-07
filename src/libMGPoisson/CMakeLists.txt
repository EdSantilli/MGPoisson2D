# Define the lib and its sources.
add_library (MGPoisson2D STATIC MGPoisson2D.f90 MGPoisson2D_DD.f90)

# This is where we should -I during compilation.
include_directories(${CMAKE_BINARY_DIR}/include)

# Put new modules into include dir.
set_target_properties(MGPoisson2D PROPERTIES Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/include)

# All targets that use this lib should also -I this include dir.
target_include_directories (MGPoisson2D PUBLIC ${CMAKE_BINARY_DIR}/include)

