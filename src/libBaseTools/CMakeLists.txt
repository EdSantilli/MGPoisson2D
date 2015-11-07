# Define the lib and its sources.
add_library (BaseTools STATIC precision.f90 arrayutils.f90)

# Put new modules into include dir.
set_target_properties(BaseTools PROPERTIES Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/include)

# All targets that use this lib should also -I this include dir.
target_include_directories (BaseTools PUBLIC ${CMAKE_BINARY_DIR}/include)
