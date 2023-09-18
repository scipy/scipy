#----------------------------------------------------------------
# Generated CMake target import file for configuration "Release".
#----------------------------------------------------------------

# Commands may need to know the format version.
set(CMAKE_IMPORT_FILE_VERSION 1)

# Import target "Vorbis::vorbis" for configuration "Release"
set_property(TARGET Vorbis::vorbis APPEND PROPERTY IMPORTED_CONFIGURATIONS RELEASE)
set_target_properties(Vorbis::vorbis PROPERTIES
  IMPORTED_LOCATION_RELEASE "${_IMPORT_PREFIX}/lib/libvorbis.so.0"
  IMPORTED_SONAME_RELEASE "libvorbis.so.0.4.9"
  )

list(APPEND _IMPORT_CHECK_TARGETS Vorbis::vorbis )
list(APPEND _IMPORT_CHECK_FILES_FOR_Vorbis::vorbis "${_IMPORT_PREFIX}/lib/libvorbis.so.0" )

# Import target "Vorbis::vorbisenc" for configuration "Release"
set_property(TARGET Vorbis::vorbisenc APPEND PROPERTY IMPORTED_CONFIGURATIONS RELEASE)
set_target_properties(Vorbis::vorbisenc PROPERTIES
  IMPORTED_LOCATION_RELEASE "${_IMPORT_PREFIX}/lib/libvorbisenc.so.2"
  IMPORTED_SONAME_RELEASE "libvorbisenc.so.2.0.12"
  )

list(APPEND _IMPORT_CHECK_TARGETS Vorbis::vorbisenc )
list(APPEND _IMPORT_CHECK_FILES_FOR_Vorbis::vorbisenc "${_IMPORT_PREFIX}/lib/libvorbisenc.so.2" )

# Import target "Vorbis::vorbisfile" for configuration "Release"
set_property(TARGET Vorbis::vorbisfile APPEND PROPERTY IMPORTED_CONFIGURATIONS RELEASE)
set_target_properties(Vorbis::vorbisfile PROPERTIES
  IMPORTED_LOCATION_RELEASE "${_IMPORT_PREFIX}/lib/libvorbisfile.so.3"
  IMPORTED_SONAME_RELEASE "libvorbisfile.so.3.3.8"
  )

list(APPEND _IMPORT_CHECK_TARGETS Vorbis::vorbisfile )
list(APPEND _IMPORT_CHECK_FILES_FOR_Vorbis::vorbisfile "${_IMPORT_PREFIX}/lib/libvorbisfile.so.3" )

# Commands beyond this point should not need to know the version.
set(CMAKE_IMPORT_FILE_VERSION)
