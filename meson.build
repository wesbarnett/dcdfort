#
# This file is part of libdcdfort
# https://github.com/wesbarnett/dcdfort
#
# Copyright (c) 2018 by James W. Barnett
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

project('dcdfort', 'Fortran', version : '1.2')
src = ['src/common.f90', 'src/dcdreader.f90', 'src/dcdwriter.f90', 'src/indexfile.f90', 'src/trajectory.f90', 'src/utils.f90'] 
main_lib = shared_library('dcdfort', src, install : true)

install_data('LICENSE', install_dir : 'share/licenes/dcdfort/')
install_data('README.md', install_dir : 'share/docs/dcdfort/')
install_subdir('docs/html', install_dir : 'share/docs/dcdfort/')
install_subdir('example', install_dir : 'share/dcdfort/')
meson.add_install_script('meson/install_modules.sh')

pkg = import('pkgconfig')
pkg.generate(
  libraries : main_lib, 
  name : meson.project_name(), 
  filebase : meson.project_name(),
  description : 'Fortran toolkit for analyzing DCD trajectory files',
  version : meson.project_version(), 
  url : 'https://github.com/wesbarnett/dcdfort')

main_dep = declare_dependency(link_with : main_lib)
test_lib = shared_library('dcdfort_tests', 'src/tests/dcdfort_tests.f90', dependencies : main_dep)
test_dep = declare_dependency(link_with : [ main_lib, test_lib ])

run_command('meson/copy_test_files.sh')

read_dcdfile = executable('read_dcdfile', 'src/tests/dcdfile.f90', dependencies : test_dep)
read_next = executable('read_next', 'src/tests/read_next.f90', dependencies : test_dep)
utils = executable('utils', 'src/tests/utils.f90', dependencies : test_dep)
ndxfile = executable('ndxfile', 'src/tests/ndxfile.f90', dependencies : test_dep)
ndxfile_grp = executable('ndxfile_grp', 'src/tests/ndxfile_grp.f90', dependencies : test_dep)

test('dcdfile read()', read_dcdfile)
test('dcdfile read_next()', read_next)
test('utils', utils)
test('ndxfile', ndxfile)
test('ndxfile groups', ndxfile_grp)
