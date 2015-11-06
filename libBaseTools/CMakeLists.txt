add_library (BaseTools STATIC precision.f90)

target_include_directories (BaseTools PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/../build/libBaseTools)

