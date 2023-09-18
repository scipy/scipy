#----------------------------------------------------------------
# Generated CMake target import file for configuration "Release".
#----------------------------------------------------------------

# Commands may need to know the format version.
set(CMAKE_IMPORT_FILE_VERSION 1)

# Import target "gr2_graphite2" for configuration "Release"
set_property(TARGET gr2_graphite2 APPEND PROPERTY IMPORTED_CONFIGURATIONS RELEASE)
set_target_properties(gr2_graphite2 PROPERTIES
  IMPORTED_LINK_INTERFACE_LIBRARIES_RELEASE "c;gcc"
  IMPORTED_LOCATION_RELEASE "${_IMPORT_PREFIX}/lib/libgraphite2.so.3.2.1"
  IMPORTED_SONAME_RELEASE "libgraphite2.so.3"
  )

list(APPEND _IMPORT_CHECK_TARGETS gr2_graphite2 )
list(APPEND _IMPORT_CHECK_FILES_FOR_gr2_graphite2 "${_IMPORT_PREFIX}/lib/libgraphite2.so.3.2.1" )

# Commands beyond this point should not need to know the version.
set(CMAKE_IMPORT_FILE_VERSION)
