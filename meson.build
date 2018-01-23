
project('dcdfort', 'Fortran')
src = ['src/common.f90', 'src/dcdreader.f90', 'src/dcdwriter.f90', 'src/indexfile.f90', 'src/trajectory.f90', 'src/utils.f90'] 
main_lib = shared_library('dcdfort', src, install : true)

main_dep = declare_dependency(link_with : main_lib)
test_lib = shared_library('dcdfort_tests', 'src/tests/dcdfort_tests.f90', dependencies : main_dep)
test_dep = declare_dependency(link_with : [ main_lib, test_lib ])

read_dcdfile = executable('read_dcdfile', 'src/tests/dcdfile.f90', dependencies : test_dep)
test('Read dcdfile', read_dcdfile)
