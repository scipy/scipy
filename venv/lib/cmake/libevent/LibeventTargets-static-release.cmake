#----------------------------------------------------------------
# Generated CMake target import file for configuration "Release".
#----------------------------------------------------------------

# Commands may need to know the format version.
set(CMAKE_IMPORT_FILE_VERSION 1)

# Import target "libevent::core" for configuration "Release"
set_property(TARGET libevent::core APPEND PROPERTY IMPORTED_CONFIGURATIONS RELEASE)
set_target_properties(libevent::core PROPERTIES
  IMPORTED_LINK_INTERFACE_LANGUAGES_RELEASE "C"
  IMPORTED_LOCATION_RELEASE "${_IMPORT_PREFIX}/lib/libevent_core.a"
  )

list(APPEND _cmake_import_check_targets libevent::core )
list(APPEND _cmake_import_check_files_for_libevent::core "${_IMPORT_PREFIX}/lib/libevent_core.a" )

# Import target "libevent::extra" for configuration "Release"
set_property(TARGET libevent::extra APPEND PROPERTY IMPORTED_CONFIGURATIONS RELEASE)
set_target_properties(libevent::extra PROPERTIES
  IMPORTED_LINK_INTERFACE_LANGUAGES_RELEASE "C"
  IMPORTED_LOCATION_RELEASE "${_IMPORT_PREFIX}/lib/libevent_extra.a"
  )

list(APPEND _cmake_import_check_targets libevent::extra )
list(APPEND _cmake_import_check_files_for_libevent::extra "${_IMPORT_PREFIX}/lib/libevent_extra.a" )

# Import target "libevent::openssl" for configuration "Release"
set_property(TARGET libevent::openssl APPEND PROPERTY IMPORTED_CONFIGURATIONS RELEASE)
set_target_properties(libevent::openssl PROPERTIES
  IMPORTED_LINK_INTERFACE_LANGUAGES_RELEASE "C"
  IMPORTED_LOCATION_RELEASE "${_IMPORT_PREFIX}/lib/libevent_openssl.a"
  )

list(APPEND _cmake_import_check_targets libevent::openssl )
list(APPEND _cmake_import_check_files_for_libevent::openssl "${_IMPORT_PREFIX}/lib/libevent_openssl.a" )

# Import target "libevent::pthreads" for configuration "Release"
set_property(TARGET libevent::pthreads APPEND PROPERTY IMPORTED_CONFIGURATIONS RELEASE)
set_target_properties(libevent::pthreads PROPERTIES
  IMPORTED_LINK_INTERFACE_LANGUAGES_RELEASE "C"
  IMPORTED_LOCATION_RELEASE "${_IMPORT_PREFIX}/lib/libevent_pthreads.a"
  )

list(APPEND _cmake_import_check_targets libevent::pthreads )
list(APPEND _cmake_import_check_files_for_libevent::pthreads "${_IMPORT_PREFIX}/lib/libevent_pthreads.a" )

# Commands beyond this point should not need to know the version.
set(CMAKE_IMPORT_FILE_VERSION)
