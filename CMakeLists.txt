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
enable_language ( Fortran C )

list ( APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake )

project ( dcdfort NONE )

####################################################################
# Set project info
####################################################################
set ( PROJECT_DESCRIPTION "Fortran toolkit for analyzing DCD trajectory files" )
set ( PROJECT_URL "https://github.com/wesbarnett/dcdfort" )
set ( LIBDIR "lib")
set ( INCLUDEDIR "include")
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
# Go ahead and set this, even though it may never be used
####################################################################
if ( NOT CMAKE_BUILD_TYPE )
    set ( CMAKE_BUILD_TYPE "Release" )
endif ()

####################################################################
# Create library
####################################################################
add_subdirectory ( src )

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
add_subdirectory ( src/tests )
