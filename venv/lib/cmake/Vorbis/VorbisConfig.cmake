
####### Expanded from @PACKAGE_INIT@ by configure_package_config_file() #######
####### Any changes to this file will be overwritten by the next CMake run ####
####### The input file was VorbisConfig.cmake.in                            ########

get_filename_component(PACKAGE_PREFIX_DIR "${CMAKE_CURRENT_LIST_DIR}/../../../" ABSOLUTE)

macro(set_and_check _var _file)
  set(${_var} "${_file}")
  if(NOT EXISTS "${_file}")
    message(FATAL_ERROR "File or directory ${_file} referenced by variable ${_var} does not exist !")
  endif()
endmacro()

macro(check_required_components _NAME)
  foreach(comp ${${_NAME}_FIND_COMPONENTS})
    if(NOT ${_NAME}_${comp}_FOUND)
      if(${_NAME}_FIND_REQUIRED_${comp})
        set(${_NAME}_FOUND FALSE)
      endif()
    endif()
  endforeach()
endmacro()

####################################################################################

include(CMakeFindDependencyMacro)
find_dependency(Ogg REQUIRED)

include(${CMAKE_CURRENT_LIST_DIR}/VorbisTargets.cmake)

set(Vorbis_Vorbis_FOUND 1)
set(Vorbis_Enc_FOUND 0)
set(Vorbis_File_FOUND 0)

if(TARGET Vorbis::vorbisenc)
    set(Vorbis_Enc_FOUND TRUE)
endif()
if(TARGET Vorbis::vorbisfile)
    set(Vorbis_File_FOUND TRUE)
endif()

check_required_components(Vorbis Enc File)
