add_library (MGPoisson2D STATIC MGPoisson2D.f90 MGPoisson2D_DD.f90)

target_include_directories (MGPoisson2D PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/../build/libMGPoisson)

