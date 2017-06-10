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
enable_testing ()
add_subdirectory ( src/tests )
