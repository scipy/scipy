# - Config file for the Libevent package
# It defines the following variables
#  LIBEVENT_FOUND            - true if libevent and all required components found on the system
#  LIBEVENT_xxx_FOUND        - true if component xxx(see available components) found on the system
#  LIBEVENT_VERSION          - libevent version in format Major.Minor.Patch
#  LIBEVENT_INCLUDE_DIRS     - directories where libevent header is located.
#  LIBEVENT_INCLUDE_DIR      - same as DIRS
#  LIBEVENT_LIBRARIES        - libevent library to link against.
#  LIBEVENT_LIBRARY          - same as LIBRARIES
#
# These variables are deprecated, don't use them.
#  LIBEVENT_STATIC_LIBRARIES - libraries to link against (archive/static)
#  LIBEVENT_SHARED_LIBRARIES - libraries to link against (shared)
#
# When you try to locate the libevent libraries, you should specify which components you want to use.
# The following table lists all available components. If none is given, all imported targets will used.
#  core        - the core functons of libevent
#  extra       - extra functions, contains http, dns and rpc
#  pthreads    - multiple threads for libevent, not exists on Windows
#  openssl     - openssl support for libevent
#
# By default, the shared libraries of libevent will be found. To find the static ones instead,
# you must set the LIBEVENT_STATIC_LINK variable to TRUE before calling find_package(Libevent ...).
# If no component provided, all components will be used.
# example:
#  set(LIBEVENT_STATIC_LINK TRUE)
#  find_package(Libevent 2.2 REQUIRED COMPONENTS core)
#  include_directories(${LIBEVENT_INCLUDE_DIRS})  # Can be omitted
#  target_link_libraries(myapp ${LIBEVENT_LIBRARIES})
#    or target_link_libraries(myapp libevent::core)
#
# find_package() can handle dependencies automatically. For example, given the 'openssl' component,
# all dependencies (libevent_core, libssl, libcrypto and openssl include directories) will be found.

set(CONFIG_FOR_INSTALL_TREE 1)

set(LIBEVENT_VERSION 2.1.12)

# IMPORTED targets from LibeventTargets.cmake
set(LIBEVENT_STATIC_LIBRARIES "core;extra;openssl;pthreads")
set(LIBEVENT_SHARED_LIBRARIES "core;extra;openssl;pthreads")

# Default to the same type as libevent was built:
if(NOT DEFINED LIBEVENT_STATIC_LINK)
    set(LIBEVENT_STATIC_LINK NOT ON)
endif()

set(CMAKE_FIND_LIBRARY_SUFFIXES_SAVE "${CMAKE_FIND_LIBRARY_SUFFIXES}")
if(${LIBEVENT_STATIC_LINK})
    set(_LIB_TYPE static)
    set(CMAKE_FIND_LIBRARY_SUFFIXES ${CMAKE_STATIC_LIBRARY_SUFFIX})
    set(_AVAILABLE_LIBS "${LIBEVENT_STATIC_LIBRARIES}")
else()
    set(_LIB_TYPE shared)
    set(CMAKE_FIND_LIBRARY_SUFFIXES ${CMAKE_SHARED_LIBRARY_SUFFIX})
    set(_AVAILABLE_LIBS "${LIBEVENT_SHARED_LIBRARIES}")
endif()

# Get the path of the current file.
get_filename_component(LIBEVENT_CMAKE_DIR "${CMAKE_CURRENT_LIST_FILE}" PATH)
get_filename_component(_INSTALL_PREFIX "${LIBEVENT_CMAKE_DIR}/../../.." ABSOLUTE)

macro(message_if_needed _flag _msg)
    if (NOT ${CMAKE_FIND_PACKAGE_NAME}_FIND_QUIETLY)
        message(${_flag} "${_msg}")
    endif()
endmacro()

macro(no_component_msg _comp)
    if(${CMAKE_FIND_PACKAGE_NAME}_FIND_REQUIRED_${_comp})
        set(pthreadlib)
        if(NOT WIN32)
            set(pthreadlib ", pthreads")
        endif()
        message(FATAL_ERROR "Your libevent library does not contain a ${_comp} component!\n"
                "The valid components are core, extra${pthreadlib} and openssl.")
    else()
        message_if_needed(WARNING "Your libevent library does not contain a ${_comp} component!")
    endif()
endmacro()

set(_EVENT_COMPONENTS)
if(${CMAKE_FIND_PACKAGE_NAME}_FIND_COMPONENTS)
    list(REMOVE_DUPLICATES ${CMAKE_FIND_PACKAGE_NAME}_FIND_COMPONENTS)
    foreach(_comp ${${CMAKE_FIND_PACKAGE_NAME}_FIND_COMPONENTS})
        list(FIND _AVAILABLE_LIBS ${_comp} _INDEX)
        if(_INDEX GREATER -1)
            list(APPEND _EVENT_COMPONENTS ${_comp})
        else()
            no_component_msg(${_comp})
        endif()
    endforeach()
else()
    set(_EVENT_COMPONENTS ${_AVAILABLE_LIBS})
endif()

set(_POSSIBLE_PKG_NAMES)
list(APPEND _POSSIBLE_PKG_NAMES ${CMAKE_FIND_PACKAGE_NAME} LIBEVENT Libevent libevent)
list(REMOVE_DUPLICATES _POSSIBLE_PKG_NAMES)

macro(set_case_insensitive_found _comp)
    foreach(name ${_POSSIBLE_PKG_NAMES})
        if("${_comp}" STREQUAL "")
            set(${name}_FOUND TRUE)
            set(${name}_NOTFOUND FALSE)
        else()
            set(${name}_${_comp}_FOUND TRUE)
            set(${name}_${_comp}_NOTFOUND FALSE)
        endif()
    endforeach()
endmacro()

if(CONFIG_FOR_INSTALL_TREE)
    ## Config for install tree ----------------------------------------
    # Find includes
    unset(_event_h CACHE)
    find_path(_event_h
              NAMES event2/event.h
              PATHS "${_INSTALL_PREFIX}/include"
              NO_DEFAULT_PATH)
    if(_event_h)
        set(LIBEVENT_INCLUDE_DIRS "${_event_h}")
        message_if_needed(STATUS "Found libevent include directory: ${_event_h}")
    else()
        message_if_needed(WARNING "Your libevent library does not contain header files!")
    endif()

    # Find libraries
    macro(find_event_lib _comp)
        unset(_event_lib CACHE)
        find_library(_event_lib
                    NAMES "event_${_comp}"
                    PATHS "${_INSTALL_PREFIX}/lib"
                    NO_DEFAULT_PATH)
        if(_event_lib)
            list(APPEND LIBEVENT_LIBRARIES "libevent::${_comp}")
            set_case_insensitive_found(${_comp})
            message_if_needed(STATUS "Found libevent component: ${_event_lib}")
        else()
            no_component_msg(${_comp})
        endif()
    endmacro()

    foreach(comp ${_EVENT_COMPONENTS})
        find_event_lib(${comp})
    endforeach()
else()
    ## Config for build tree ----------------------------------------
    set(LIBEVENT_INCLUDE_DIRS "/home/conda/feedstock_root/build_artifacts/libevent_1685725798133/work/include;/home/conda/feedstock_root/build_artifacts/libevent_1685725798133/work/build/include")
    foreach(_comp ${_EVENT_COMPONENTS})
        list(APPEND LIBEVENT_LIBRARIES "libevent::${_comp}")
        set_case_insensitive_found(${_comp})
    endforeach()
endif()

set(LIBEVENT_INCLUDE_DIR ${LIBEVENT_INCLUDE_DIRS})
if(LIBEVENT_LIBRARIES)
    set(LIBEVENT_LIBRARY ${LIBEVENT_LIBRARIES})
    if(CONFIG_FOR_INSTALL_TREE)
        message_if_needed(STATUS "Found libevent ${LIBEVENT_VERSION} in ${_INSTALL_PREFIX}")
    else()
        message_if_needed(STATUS "Found libevent ${LIBEVENT_VERSION} in ${LIBEVENT_CMAKE_DIR}")
    endif()

    # Avoid including targets more than one times
    if(NOT TARGET event_core_${_LIB_TYPE})
        # Include the project Targets file, this contains definitions for IMPORTED targets.
        include(${LIBEVENT_CMAKE_DIR}/LibeventTargets-${_LIB_TYPE}.cmake)
    endif()
else()
    if(${CMAKE_FIND_PACKAGE_NAME}_FIND_REQUIRED)
        message(FATAL_ERROR "Can not find any libraries for libevent.")
    else()
        message_if_needed(WARNING "Can not find any libraries for libevent.")
    endif()
endif()

set(CMAKE_FIND_LIBRARY_SUFFIXES "${CMAKE_FIND_LIBRARY_SUFFIXES_SAVE}")
unset(_LIB_TYPE)
unset(_AVAILABLE_LIBS)
unset(_EVENT_COMPONENTS)
unset(_POSSIBLE_PKG_NAMES)
unset(_INSTALL_PREFIX)
