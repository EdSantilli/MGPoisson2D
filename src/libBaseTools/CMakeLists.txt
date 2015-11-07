add_library (BaseTools STATIC precision.f90)
set_target_properties(BaseTools PROPERTIES Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/include)
target_include_directories (BaseTools PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
