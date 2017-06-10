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

# The main library to build
add_library ( ${CMAKE_PROJECT_NAME} SHARED dcdplugin.c trajectory.f90 utils.f90 )

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
