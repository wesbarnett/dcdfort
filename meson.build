
project('dcdfort', 'Fortran')
src = ['src/common.f90', 'src/dcdreader.f90', 'src/dcdwriter.f90', 'src/indexfile.f90', 'src/trajectory.f90', 'src/utils.f90'] 
main_lib = shared_library('dcdfort', src, install : true)

main_dep = declare_dependency(link_with : main_lib)
test_lib = shared_library('dcdfort_tests', 'src/tests/dcdfort_tests.f90', dependencies : main_dep)
test_dep = declare_dependency(link_with : [ main_lib, test_lib ])

testdcd_file = '/'.join([meson.source_root(),'src/tests/test.dcd'])
testdcd_file_out = '/'.join([meson.build_root(),'test.dcd'])
run_command('cp', testdcd_file , testdcd_file_out)

testndx_file = '/'.join([meson.source_root(),'src/tests/test.ndx'])
testndx_file_out = '/'.join([meson.build_root(),'test.ndx'])
run_command('cp', testndx_file , testndx_file_out)

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
