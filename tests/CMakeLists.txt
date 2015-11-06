# Add executable called "helloDemo" that is built from the source files
# "demo.cxx" and "demo_b.cxx". The extensions are automatically found.
add_executable (test test.f90)

# Link the executable to the Hello library. Since the Hello library has
# public include directories we will use those link directories when building
# helloDemo
target_link_libraries (test LINK_PUBLIC BaseTools MGPoisson2D)

# include_directories ($(CMAKE_CURRENT_SOURCE_DIR)/../build/libMGPoisson)
