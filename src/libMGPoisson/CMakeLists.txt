add_library (MGPoisson2D STATIC MGPoisson2D.f90 MGPoisson2D_DD.f90)
include_directories(${CMAKE_BINARY_DIR}/include)
set_target_properties(MGPoisson2D PROPERTIES Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/include)
target_include_directories (MGPoisson2D PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

