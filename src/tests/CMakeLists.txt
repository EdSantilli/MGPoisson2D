# Define the executable and its source.
add_executable (test test.f90)

# This is where we should -I during compilation.
include_directories(${CMAKE_BINARY_DIR}/include)

# Link against these libs.
target_link_libraries (test LINK_PUBLIC BaseTools MGPoisson2D)

