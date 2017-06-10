# Include the module files from the built libdcdfort.so
include_directories( ${CMAKE_BINARY_DIR}/include )

# Create a library for the tests and link it to the built libdcdfort.so
add_library ( dcdfort_tests SHARED "dcdfort_tests.f90" )
target_link_libraries ( dcdfort_tests ${CMAKE_PROJECT_NAME} )

# Create the test programs and link them to the test library
foreach ( TESTNAME dcdfile )
    add_executable ( ${TESTNAME} ${TESTNAME}.f90 )
    target_link_libraries ( ${TESTNAME} dcdfort_tests )
endforeach ( TESTNAME )

add_test ( "Read_dcdfile" dcdfile )

# Copy over files needed for running the tests
configure_file ( test.dcd ${CMAKE_CURRENT_BINARY_DIR}/test.dcd COPYONLY )
