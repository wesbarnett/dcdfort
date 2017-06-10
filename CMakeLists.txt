
cmake_minimum_required ( VERSION 2.8.8 )
enable_language ( Fortran C )

list ( APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake )

project ( dcdfort NONE )

# Set project info
include ( info )

# Set package version
include ( version )

# Include C header files for compilation
include_directories( include )

# Set build type
include ( buildtype )

# Create library
add_subdirectory ( src )

# Publicize installed location to other CMake projects and pkgconfig
include ( public )

# Tests
# TODO
#enable_testing ()
#add_subdirectory ( src/tests )
