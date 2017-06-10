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
