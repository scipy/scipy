if (CMAKE_VERSION VERSION_LESS 3.1.0)
    message(FATAL_ERROR "Qt 5 RemoteObjects module requires at least CMake version 3.1.0")
endif()

get_filename_component(_qt5RemoteObjects_install_prefix "${CMAKE_CURRENT_LIST_DIR}/../../../" ABSOLUTE)

# For backwards compatibility only. Use Qt5RemoteObjects_VERSION instead.
set(Qt5RemoteObjects_VERSION_STRING 5.15.8)

set(Qt5RemoteObjects_LIBRARIES Qt5::RemoteObjects)

macro(_qt5_RemoteObjects_check_file_exists file)
    if(NOT EXISTS "${file}" )
        message(FATAL_ERROR "The imported target \"Qt5::RemoteObjects\" references the file
   \"${file}\"
but this file does not exist.  Possible reasons include:
* The file was deleted, renamed, or moved to another location.
* An install or uninstall procedure did not complete successfully.
* The installation package was faulty and contained
   \"${CMAKE_CURRENT_LIST_FILE}\"
but not all the files it references.
")
    endif()
endmacro()


macro(_populate_RemoteObjects_target_properties Configuration LIB_LOCATION IMPLIB_LOCATION
      IsDebugAndRelease)
    set_property(TARGET Qt5::RemoteObjects APPEND PROPERTY IMPORTED_CONFIGURATIONS ${Configuration})

    set(imported_location "${_qt5RemoteObjects_install_prefix}/lib/${LIB_LOCATION}")
    _qt5_RemoteObjects_check_file_exists(${imported_location})
    set(_deps
        ${_Qt5RemoteObjects_LIB_DEPENDENCIES}
    )
    set(_static_deps
    )

    set_target_properties(Qt5::RemoteObjects PROPERTIES
        "IMPORTED_LOCATION_${Configuration}" ${imported_location}
        "IMPORTED_SONAME_${Configuration}" "libQt5RemoteObjects.so.5"
        # For backward compatibility with CMake < 2.8.12
        "IMPORTED_LINK_INTERFACE_LIBRARIES_${Configuration}" "${_deps};${_static_deps}"
    )
    set_property(TARGET Qt5::RemoteObjects APPEND PROPERTY INTERFACE_LINK_LIBRARIES
                 "${_deps}"
    )


endmacro()

if (NOT TARGET Qt5::RemoteObjects)

    set(_Qt5RemoteObjects_OWN_INCLUDE_DIRS "${_qt5RemoteObjects_install_prefix}/include/qt/" "${_qt5RemoteObjects_install_prefix}/include/qt/QtRemoteObjects")
    set(Qt5RemoteObjects_PRIVATE_INCLUDE_DIRS
        "${_qt5RemoteObjects_install_prefix}/include/qt/QtRemoteObjects/5.15.8"
        "${_qt5RemoteObjects_install_prefix}/include/qt/QtRemoteObjects/5.15.8/QtRemoteObjects"
    )
    include("${CMAKE_CURRENT_LIST_DIR}/ExtraSourceIncludes.cmake" OPTIONAL)

    foreach(_dir ${_Qt5RemoteObjects_OWN_INCLUDE_DIRS})
        _qt5_RemoteObjects_check_file_exists(${_dir})
    endforeach()

    # Only check existence of private includes if the Private component is
    # specified.
    list(FIND Qt5RemoteObjects_FIND_COMPONENTS Private _check_private)
    if (NOT _check_private STREQUAL -1)
        foreach(_dir ${Qt5RemoteObjects_PRIVATE_INCLUDE_DIRS})
            _qt5_RemoteObjects_check_file_exists(${_dir})
        endforeach()
    endif()

    set(Qt5RemoteObjects_INCLUDE_DIRS ${_Qt5RemoteObjects_OWN_INCLUDE_DIRS})

    set(Qt5RemoteObjects_DEFINITIONS -DQT_REMOTEOBJECTS_LIB)
    set(Qt5RemoteObjects_COMPILE_DEFINITIONS QT_REMOTEOBJECTS_LIB)
    set(_Qt5RemoteObjects_MODULE_DEPENDENCIES "Network;Core")


    set(Qt5RemoteObjects_OWN_PRIVATE_INCLUDE_DIRS ${Qt5RemoteObjects_PRIVATE_INCLUDE_DIRS})

    set(_Qt5RemoteObjects_FIND_DEPENDENCIES_REQUIRED)
    if (Qt5RemoteObjects_FIND_REQUIRED)
        set(_Qt5RemoteObjects_FIND_DEPENDENCIES_REQUIRED REQUIRED)
    endif()
    set(_Qt5RemoteObjects_FIND_DEPENDENCIES_QUIET)
    if (Qt5RemoteObjects_FIND_QUIETLY)
        set(_Qt5RemoteObjects_DEPENDENCIES_FIND_QUIET QUIET)
    endif()
    set(_Qt5RemoteObjects_FIND_VERSION_EXACT)
    if (Qt5RemoteObjects_FIND_VERSION_EXACT)
        set(_Qt5RemoteObjects_FIND_VERSION_EXACT EXACT)
    endif()

    set(Qt5RemoteObjects_EXECUTABLE_COMPILE_FLAGS "")

    foreach(_module_dep ${_Qt5RemoteObjects_MODULE_DEPENDENCIES})
        if (NOT Qt5${_module_dep}_FOUND)
            find_package(Qt5${_module_dep}
                5.15.8 ${_Qt5RemoteObjects_FIND_VERSION_EXACT}
                ${_Qt5RemoteObjects_DEPENDENCIES_FIND_QUIET}
                ${_Qt5RemoteObjects_FIND_DEPENDENCIES_REQUIRED}
                PATHS "${CMAKE_CURRENT_LIST_DIR}/.." NO_DEFAULT_PATH
            )
        endif()

        if (NOT Qt5${_module_dep}_FOUND)
            set(Qt5RemoteObjects_FOUND False)
            return()
        endif()

        list(APPEND Qt5RemoteObjects_INCLUDE_DIRS "${Qt5${_module_dep}_INCLUDE_DIRS}")
        list(APPEND Qt5RemoteObjects_PRIVATE_INCLUDE_DIRS "${Qt5${_module_dep}_PRIVATE_INCLUDE_DIRS}")
        list(APPEND Qt5RemoteObjects_DEFINITIONS ${Qt5${_module_dep}_DEFINITIONS})
        list(APPEND Qt5RemoteObjects_COMPILE_DEFINITIONS ${Qt5${_module_dep}_COMPILE_DEFINITIONS})
        list(APPEND Qt5RemoteObjects_EXECUTABLE_COMPILE_FLAGS ${Qt5${_module_dep}_EXECUTABLE_COMPILE_FLAGS})
    endforeach()
    list(REMOVE_DUPLICATES Qt5RemoteObjects_INCLUDE_DIRS)
    list(REMOVE_DUPLICATES Qt5RemoteObjects_PRIVATE_INCLUDE_DIRS)
    list(REMOVE_DUPLICATES Qt5RemoteObjects_DEFINITIONS)
    list(REMOVE_DUPLICATES Qt5RemoteObjects_COMPILE_DEFINITIONS)
    list(REMOVE_DUPLICATES Qt5RemoteObjects_EXECUTABLE_COMPILE_FLAGS)

    # It can happen that the same FooConfig.cmake file is included when calling find_package()
    # on some Qt component. An example of that is when using a Qt static build with auto inclusion
    # of plugins:
    #
    # Qt5WidgetsConfig.cmake -> Qt5GuiConfig.cmake -> Qt5Gui_QSvgIconPlugin.cmake ->
    # Qt5SvgConfig.cmake -> Qt5WidgetsConfig.cmake ->
    # finish processing of second Qt5WidgetsConfig.cmake ->
    # return to first Qt5WidgetsConfig.cmake ->
    # add_library cannot create imported target Qt5::Widgets.
    #
    # Make sure to return early in the original Config inclusion, because the target has already
    # been defined as part of the second inclusion.
    if(TARGET Qt5::RemoteObjects)
        return()
    endif()

    set(_Qt5RemoteObjects_LIB_DEPENDENCIES "Qt5::Network;Qt5::Core")


    add_library(Qt5::RemoteObjects SHARED IMPORTED)


    set_property(TARGET Qt5::RemoteObjects PROPERTY
      INTERFACE_INCLUDE_DIRECTORIES ${_Qt5RemoteObjects_OWN_INCLUDE_DIRS})
    set_property(TARGET Qt5::RemoteObjects PROPERTY
      INTERFACE_COMPILE_DEFINITIONS QT_REMOTEOBJECTS_LIB)

    set_property(TARGET Qt5::RemoteObjects PROPERTY INTERFACE_QT_ENABLED_FEATURES )
    set_property(TARGET Qt5::RemoteObjects PROPERTY INTERFACE_QT_DISABLED_FEATURES )

    # Qt 6 forward compatible properties.
    set_property(TARGET Qt5::RemoteObjects
                 PROPERTY QT_ENABLED_PUBLIC_FEATURES
                 )
    set_property(TARGET Qt5::RemoteObjects
                 PROPERTY QT_DISABLED_PUBLIC_FEATURES
                 )
    set_property(TARGET Qt5::RemoteObjects
                 PROPERTY QT_ENABLED_PRIVATE_FEATURES
                 )
    set_property(TARGET Qt5::RemoteObjects
                 PROPERTY QT_DISABLED_PRIVATE_FEATURES
                 )

    set_property(TARGET Qt5::RemoteObjects PROPERTY INTERFACE_QT_PLUGIN_TYPES "")

    set(_Qt5RemoteObjects_PRIVATE_DIRS_EXIST TRUE)
    foreach (_Qt5RemoteObjects_PRIVATE_DIR ${Qt5RemoteObjects_OWN_PRIVATE_INCLUDE_DIRS})
        if (NOT EXISTS ${_Qt5RemoteObjects_PRIVATE_DIR})
            set(_Qt5RemoteObjects_PRIVATE_DIRS_EXIST FALSE)
        endif()
    endforeach()

    if (_Qt5RemoteObjects_PRIVATE_DIRS_EXIST)
        add_library(Qt5::RemoteObjectsPrivate INTERFACE IMPORTED)
        set_property(TARGET Qt5::RemoteObjectsPrivate PROPERTY
            INTERFACE_INCLUDE_DIRECTORIES ${Qt5RemoteObjects_OWN_PRIVATE_INCLUDE_DIRS}
        )
        set(_Qt5RemoteObjects_PRIVATEDEPS)
        foreach(dep ${_Qt5RemoteObjects_LIB_DEPENDENCIES})
            if (TARGET ${dep}Private)
                list(APPEND _Qt5RemoteObjects_PRIVATEDEPS ${dep}Private)
            endif()
        endforeach()
        set_property(TARGET Qt5::RemoteObjectsPrivate PROPERTY
            INTERFACE_LINK_LIBRARIES Qt5::RemoteObjects ${_Qt5RemoteObjects_PRIVATEDEPS}
        )

        # Add a versionless target, for compatibility with Qt6.
        if(NOT "${QT_NO_CREATE_VERSIONLESS_TARGETS}" AND NOT TARGET Qt::RemoteObjectsPrivate)
            add_library(Qt::RemoteObjectsPrivate INTERFACE IMPORTED)
            set_target_properties(Qt::RemoteObjectsPrivate PROPERTIES
                INTERFACE_LINK_LIBRARIES "Qt5::RemoteObjectsPrivate"
            )
        endif()
    endif()

    _populate_RemoteObjects_target_properties(RELEASE "libQt5RemoteObjects.so.5.15.8" "" FALSE)




    # In Qt 5.15 the glob pattern was relaxed to also catch plugins not literally named Plugin.
    # Define QT5_STRICT_PLUGIN_GLOB or ModuleName_STRICT_PLUGIN_GLOB to revert to old behavior.
    if (QT5_STRICT_PLUGIN_GLOB OR Qt5RemoteObjects_STRICT_PLUGIN_GLOB)
        file(GLOB pluginTargets "${CMAKE_CURRENT_LIST_DIR}/Qt5RemoteObjects_*Plugin.cmake")
    else()
        file(GLOB pluginTargets "${CMAKE_CURRENT_LIST_DIR}/Qt5RemoteObjects_*.cmake")
    endif()

    macro(_populate_RemoteObjects_plugin_properties Plugin Configuration PLUGIN_LOCATION
          IsDebugAndRelease)
        set_property(TARGET Qt5::${Plugin} APPEND PROPERTY IMPORTED_CONFIGURATIONS ${Configuration})

        set(imported_location "${_qt5RemoteObjects_install_prefix}/plugins/${PLUGIN_LOCATION}")
        _qt5_RemoteObjects_check_file_exists(${imported_location})
        set_target_properties(Qt5::${Plugin} PROPERTIES
            "IMPORTED_LOCATION_${Configuration}" ${imported_location}
        )

    endmacro()

    if (pluginTargets)
        foreach(pluginTarget ${pluginTargets})
            include(${pluginTarget})
        endforeach()
    endif()

    include("${CMAKE_CURRENT_LIST_DIR}/Qt5RemoteObjectsConfigExtras.cmake")

    include("${CMAKE_CURRENT_LIST_DIR}/Qt5RemoteObjectsMacros.cmake")

    _qt5_RemoteObjects_check_file_exists("${CMAKE_CURRENT_LIST_DIR}/Qt5RemoteObjectsConfigVersion.cmake")
endif()

# Add a versionless target, for compatibility with Qt6.
if(NOT "${QT_NO_CREATE_VERSIONLESS_TARGETS}" AND TARGET Qt5::RemoteObjects AND NOT TARGET Qt::RemoteObjects)
    add_library(Qt::RemoteObjects INTERFACE IMPORTED)
    set_target_properties(Qt::RemoteObjects PROPERTIES
        INTERFACE_LINK_LIBRARIES "Qt5::RemoteObjects"
    )
endif()
