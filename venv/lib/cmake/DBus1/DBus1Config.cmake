# - Config file for the DBus1 package
# It defines the following variables
#  DBus1_FOUND - Flag for indicating that DBus1 package has been found
#  DBus1_DEFINITIONS  - compile definitions for DBus1 [1]
#  DBus1_INCLUDE_DIRS - include directories for DBus1 [1]
#  DBus1_LIBRARIES    - cmake targets to link against

# [1] This variable is not required if DBus1_LIBRARIES is added
#     to a target with target_link_libraries

get_filename_component(DBus1_PKGCONFIG_DIR "${CMAKE_CURRENT_LIST_DIR}/../../pkgconfig" ABSOLUTE)
find_package(PkgConfig)
set(ENV{PKG_CONFIG_DIR})
set(ENV{PKG_CONFIG_PATH} ${DBus1_PKGCONFIG_DIR})
set(ENV{PKG_CONFIG_LIBDIR} ${DBus1_PKGCONFIG_DIR})
# for debugging
#set(ENV{PKG_CONFIG_DEBUG_SPEW} 1)
pkg_check_modules(PC_DBUS1 QUIET dbus-1)
set(DBus1_DEFINITIONS ${PC_DBUS1_CFLAGS_OTHER})

# find the real stuff and use pkgconfig variables as hints
# because cmake provides more search options
find_path(DBus1_INCLUDE_DIR dbus/dbus.h
    HINTS ${PC_DBUS1_INCLUDEDIR} ${PC_DBUS1_INCLUDE_DIRS}
    PATH_SUFFIXES dbus-1.0)
find_path(DBus1_ARCH_INCLUDE_DIR dbus/dbus-arch-deps.h
    HINTS ${PC_DBUS1_INCLUDE_DIRS}
    PATH_SUFFIXES dbus-1.0)
find_library(DBus1_LIBRARY NAMES ${PC_DBUS1_LIBRARIES}
    HINTS ${PC_DBUS1_LIBDIR} ${PC_DBUS1_LIBRARY_DIRS})

include(FindPackageHandleStandardArgs)
# handle the QUIETLY and REQUIRED arguments and set DBus1_FOUND to TRUE
# if all listed variables are TRUE
find_package_handle_standard_args(DBus1 DEFAULT_MSG
    DBus1_LIBRARY DBus1_INCLUDE_DIR DBus1_ARCH_INCLUDE_DIR)

# make the mentioned variables only visible in cmake gui with "advanced" enabled
mark_as_advanced(DBus1_INCLUDE_DIR DBus1_LIBRARY)

set(DBus1_LIBRARIES dbus-1)
set(DBus1_INCLUDE_DIRS "${DBus1_INCLUDE_DIR}" "${DBus1_ARCH_INCLUDE_DIR}")

# setup imported target
add_library(dbus-1 SHARED IMPORTED)
set_property(TARGET dbus-1 APPEND PROPERTY IMPORTED_LOCATION ${DBus1_LIBRARY})
set_property(TARGET dbus-1 APPEND PROPERTY INTERFACE_INCLUDE_DIRECTORIES ${DBus1_INCLUDE_DIRS})
set_property(TARGET dbus-1 APPEND PROPERTY INTERFACE_COMPILE_DEFINITIONS ${DBus1_DEFINITIONS})
