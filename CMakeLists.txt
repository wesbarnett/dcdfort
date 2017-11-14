#
# This file is part of libdcdfort
# https://github.com/wesbarnett/dcdfort
#
# Copyright (c) 2017 by James W. Barnett
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License along
# with this program; if not, write to the Free Software Foundation, Inc.,
# 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.

cmake_minimum_required ( VERSION 2.8.8 )
enable_language ( Fortran )

list ( APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake )

project ( dcdfort NONE )

####################################################################
# Check if GNU Fortran compiler is used
####################################################################
if ( NOT "${CMAKE_Fortran_COMPILER_ID}" STREQUAL "GNU" )
    message ( FATAL_ERROR "${CMAKE_Fortran_COMPILER_ID} compiler is not supported. GNU Fortran compiler is required.")
endif ()

####################################################################
# Set project info
####################################################################
set ( PROJECT_DESCRIPTION "Fortran toolkit for analyzing DCD trajectory files" )
set ( PROJECT_URL "https://github.com/wesbarnett/dcdfort" )
set ( LIBDIR "lib")
set ( INCLUDEDIR "include")
set ( TESTDIR ${CMAKE_CURRENT_SOURCE_DIR}/src/tests/ )
set ( INSTALL_LIBDIR ${CMAKE_INSTALL_PREFIX}/${LIBDIR} )
set ( INSTALL_INCLUDEDIR ${CMAKE_INSTALL_PREFIX}/${INCLUDEDIR} )

####################################################################
# Version number is stored in a file
####################################################################
file ( STRINGS "${CMAKE_SOURCE_DIR}/.VERSION" VERSION )
string ( REPLACE "." ";" VERSION_LIST ${VERSION} )
list ( GET VERSION_LIST 0 VERSION_MAJOR )
list ( GET VERSION_LIST 1 VERSION_MINOR )
list ( GET VERSION_LIST 2 VERSION_PATCH )
set ( PROJECT_VERSION "${VERSION_MAJOR}.${VERSION_MINOR}.${VERSION_PATCH}" )
message ( STATUS "CMake build configuration for ${CMAKE_PROJECT_NAME} ${PROJECT_VERSION}" )
string ( TOLOWER ${CMAKE_PROJECT_NAME} PACKAGE_NAME )
set ( PACKAGE_VERSION "${PACKAGE_NAME}-${VERSION}" )

####################################################################
# Set flags based on build type
####################################################################
if ( NOT CMAKE_BUILD_TYPE )
  set ( CMAKE_BUILD_TYPE "Release" )
endif ()

set ( CMAKE_Fortran_FLAGS_RELEASE "-g -O3 -Wall -fno-f2c" )
set ( CMAKE_Fortran_FLAGS_DEBUG "-g -O0 -Wall -Wextra -Wunreachable-code -Wimplicit-interface -Wimplicit-procedure -pedantic -fcheck=all" )

####################################################################
# Create library
####################################################################
# The main library to build
file ( GLOB PKG_SOURCES ${CMAKE_SOURCE_DIR}/src/*.f90 )

add_library ( ${CMAKE_PROJECT_NAME} SHARED ${PKG_SOURCES} ) 

# Where to put files after they are built, and what to name them
set_target_properties ( ${CMAKE_PROJECT_NAME}
  PROPERTIES
  OUTPUT_NAME ${CMAKE_PROJECT_NAME}
  PREFIX lib
  SOVERSION ${VERSION_MAJOR}.${VERSION_MINOR}
  VERSION ${VERSION}
  LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib"
  Fortran_MODULE_DIRECTORY "${CMAKE_BINARY_DIR}/include" )

# Install the compiled library
install ( TARGETS ${CMAKE_PROJECT_NAME} LIBRARY DESTINATION lib )

# Install the module files (in include directory)
install ( DIRECTORY "${CMAKE_BINARY_DIR}/include" DESTINATION "${CMAKE_INSTALL_PREFIX}" )

# Install the example program
install ( DIRECTORY "${CMAKE_SOURCE_DIR}/example" DESTINATION
    "${CMAKE_INSTALL_PREFIX}/share/${PACKAGE_NAME}/" )

# Install the license
install ( FILES "${CMAKE_SOURCE_DIR}/LICENSE" DESTINATION
    "${CMAKE_INSTALL_PREFIX}/share/licenes/${PACKAGE_NAME}/" )

# Install the docs
install ( DIRECTORY "${CMAKE_SOURCE_DIR}/docs/html" DESTINATION
    "${CMAKE_INSTALL_PREFIX}/share/docs/${PACKAGE_NAME}/" )
install ( FILES "${CMAKE_SOURCE_DIR}/README.md" DESTINATION
    "${CMAKE_INSTALL_PREFIX}/share/docs/${PACKAGE_NAME}/" )

####################################################################
# Publicize installed location to other CMake projects and pkgconfig
####################################################################
set ( EXPORT_INSTALL_DIR "${INSTALL_LIBDIR}/cmake/${PACKAGE_VERSION}" )

include ( CMakePackageConfigHelpers ) # Standard CMake module
write_basic_package_version_file(
"${CMAKE_BINARY_DIR}/${PACKAGE_NAME}-config-version.cmake"
  VERSION ${VERSION}
  COMPATIBILITY SameMajorVersion )

# install package config file
configure_package_config_file (
  "${CMAKE_SOURCE_DIR}/cmake/${CMAKE_PROJECT_NAME}-config.cmake.in"
  "${CMAKE_BINARY_DIR}/pkg/${PACKAGE_NAME}-config.cmake"
  INSTALL_DESTINATION "${EXPORT_INSTALL_DIR}"
  PATH_VARS EXPORT_INSTALL_DIR )

# Install the config and version files so that we can find this project with
# others
install ( FILES
  "${CMAKE_BINARY_DIR}/pkg/${PACKAGE_NAME}-config.cmake"
  "${CMAKE_BINARY_DIR}/${PACKAGE_NAME}-config-version.cmake"
  DESTINATION "${EXPORT_INSTALL_DIR}" )

# pkg-config
configure_file(
   "${CMAKE_CURRENT_SOURCE_DIR}/cmake/${CMAKE_PROJECT_NAME}.pc.cmake.in"
   "${CMAKE_CURRENT_BINARY_DIR}/lib${CMAKE_PROJECT_NAME}.pc"
   @ONLY )
install ( FILES
   "${CMAKE_CURRENT_BINARY_DIR}/lib${CMAKE_PROJECT_NAME}.pc"
   DESTINATION "${LIBDIR}/pkgconfig" )

####################################################################
# Tests
####################################################################
enable_testing ()
# Include the module files from the built libdcdfort.so
include_directories( ${CMAKE_BINARY_DIR}/include )

# Create a library for the tests and link it to the built libdcdfort.so
add_library ( dcdfort_tests SHARED "${CMAKE_CURRENT_SOURCE_DIR}/src/tests/dcdfort_tests.f90" )
target_link_libraries ( dcdfort_tests ${CMAKE_PROJECT_NAME} )

# Create the test programs and link them to the test library
foreach ( TESTNAME dcdfile read_next utils ndxfile ndxfile_grp speed )
    add_executable ( ${TESTNAME} ${TESTDIR}/${TESTNAME}.f90 )
    target_link_libraries ( ${TESTNAME} dcdfort_tests )
endforeach ( TESTNAME )

add_test ( "Read_dcdfile" dcdfile )
add_test ( "Read_next" read_next )
add_test ( "Utils" utils )
add_test ( "Index_file" ndxfile )
add_test ( "Index_file_groups" ndxfile_grp )

# Copy over files needed for running the tests
configure_file ( ${TESTDIR}/test.dcd ${CMAKE_CURRENT_BINARY_DIR}/test.dcd COPYONLY )
configure_file ( ${TESTDIR}/test.ndx ${CMAKE_CURRENT_BINARY_DIR}/test.ndx COPYONLY )
