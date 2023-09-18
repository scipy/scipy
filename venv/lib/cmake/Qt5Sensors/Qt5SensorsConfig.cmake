if (CMAKE_VERSION VERSION_LESS 3.1.0)
    message(FATAL_ERROR "Qt 5 Sensors module requires at least CMake version 3.1.0")
endif()

get_filename_component(_qt5Sensors_install_prefix "${CMAKE_CURRENT_LIST_DIR}/../../../" ABSOLUTE)

# For backwards compatibility only. Use Qt5Sensors_VERSION instead.
set(Qt5Sensors_VERSION_STRING 5.15.8)

set(Qt5Sensors_LIBRARIES Qt5::Sensors)

macro(_qt5_Sensors_check_file_exists file)
    if(NOT EXISTS "${file}" )
        message(FATAL_ERROR "The imported target \"Qt5::Sensors\" references the file
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


macro(_populate_Sensors_target_properties Configuration LIB_LOCATION IMPLIB_LOCATION
      IsDebugAndRelease)
    set_property(TARGET Qt5::Sensors APPEND PROPERTY IMPORTED_CONFIGURATIONS ${Configuration})

    set(imported_location "${_qt5Sensors_install_prefix}/lib/${LIB_LOCATION}")
    _qt5_Sensors_check_file_exists(${imported_location})
    set(_deps
        ${_Qt5Sensors_LIB_DEPENDENCIES}
    )
    set(_static_deps
    )

    set_target_properties(Qt5::Sensors PROPERTIES
        "IMPORTED_LOCATION_${Configuration}" ${imported_location}
        "IMPORTED_SONAME_${Configuration}" "libQt5Sensors.so.5"
        # For backward compatibility with CMake < 2.8.12
        "IMPORTED_LINK_INTERFACE_LIBRARIES_${Configuration}" "${_deps};${_static_deps}"
    )
    set_property(TARGET Qt5::Sensors APPEND PROPERTY INTERFACE_LINK_LIBRARIES
                 "${_deps}"
    )


endmacro()

if (NOT TARGET Qt5::Sensors)

    set(_Qt5Sensors_OWN_INCLUDE_DIRS "${_qt5Sensors_install_prefix}/include/qt/" "${_qt5Sensors_install_prefix}/include/qt/QtSensors")
    set(Qt5Sensors_PRIVATE_INCLUDE_DIRS
        "${_qt5Sensors_install_prefix}/include/qt/QtSensors/5.15.8"
        "${_qt5Sensors_install_prefix}/include/qt/QtSensors/5.15.8/QtSensors"
    )
    include("${CMAKE_CURRENT_LIST_DIR}/ExtraSourceIncludes.cmake" OPTIONAL)

    foreach(_dir ${_Qt5Sensors_OWN_INCLUDE_DIRS})
        _qt5_Sensors_check_file_exists(${_dir})
    endforeach()

    # Only check existence of private includes if the Private component is
    # specified.
    list(FIND Qt5Sensors_FIND_COMPONENTS Private _check_private)
    if (NOT _check_private STREQUAL -1)
        foreach(_dir ${Qt5Sensors_PRIVATE_INCLUDE_DIRS})
            _qt5_Sensors_check_file_exists(${_dir})
        endforeach()
    endif()

    set(Qt5Sensors_INCLUDE_DIRS ${_Qt5Sensors_OWN_INCLUDE_DIRS})

    set(Qt5Sensors_DEFINITIONS -DQT_SENSORS_LIB)
    set(Qt5Sensors_COMPILE_DEFINITIONS QT_SENSORS_LIB)
    set(_Qt5Sensors_MODULE_DEPENDENCIES "Core")


    set(Qt5Sensors_OWN_PRIVATE_INCLUDE_DIRS ${Qt5Sensors_PRIVATE_INCLUDE_DIRS})

    set(_Qt5Sensors_FIND_DEPENDENCIES_REQUIRED)
    if (Qt5Sensors_FIND_REQUIRED)
        set(_Qt5Sensors_FIND_DEPENDENCIES_REQUIRED REQUIRED)
    endif()
    set(_Qt5Sensors_FIND_DEPENDENCIES_QUIET)
    if (Qt5Sensors_FIND_QUIETLY)
        set(_Qt5Sensors_DEPENDENCIES_FIND_QUIET QUIET)
    endif()
    set(_Qt5Sensors_FIND_VERSION_EXACT)
    if (Qt5Sensors_FIND_VERSION_EXACT)
        set(_Qt5Sensors_FIND_VERSION_EXACT EXACT)
    endif()

    set(Qt5Sensors_EXECUTABLE_COMPILE_FLAGS "")

    foreach(_module_dep ${_Qt5Sensors_MODULE_DEPENDENCIES})
        if (NOT Qt5${_module_dep}_FOUND)
            find_package(Qt5${_module_dep}
                5.15.8 ${_Qt5Sensors_FIND_VERSION_EXACT}
                ${_Qt5Sensors_DEPENDENCIES_FIND_QUIET}
                ${_Qt5Sensors_FIND_DEPENDENCIES_REQUIRED}
                PATHS "${CMAKE_CURRENT_LIST_DIR}/.." NO_DEFAULT_PATH
            )
        endif()

        if (NOT Qt5${_module_dep}_FOUND)
            set(Qt5Sensors_FOUND False)
            return()
        endif()

        list(APPEND Qt5Sensors_INCLUDE_DIRS "${Qt5${_module_dep}_INCLUDE_DIRS}")
        list(APPEND Qt5Sensors_PRIVATE_INCLUDE_DIRS "${Qt5${_module_dep}_PRIVATE_INCLUDE_DIRS}")
        list(APPEND Qt5Sensors_DEFINITIONS ${Qt5${_module_dep}_DEFINITIONS})
        list(APPEND Qt5Sensors_COMPILE_DEFINITIONS ${Qt5${_module_dep}_COMPILE_DEFINITIONS})
        list(APPEND Qt5Sensors_EXECUTABLE_COMPILE_FLAGS ${Qt5${_module_dep}_EXECUTABLE_COMPILE_FLAGS})
    endforeach()
    list(REMOVE_DUPLICATES Qt5Sensors_INCLUDE_DIRS)
    list(REMOVE_DUPLICATES Qt5Sensors_PRIVATE_INCLUDE_DIRS)
    list(REMOVE_DUPLICATES Qt5Sensors_DEFINITIONS)
    list(REMOVE_DUPLICATES Qt5Sensors_COMPILE_DEFINITIONS)
    list(REMOVE_DUPLICATES Qt5Sensors_EXECUTABLE_COMPILE_FLAGS)

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
    if(TARGET Qt5::Sensors)
        return()
    endif()

    set(_Qt5Sensors_LIB_DEPENDENCIES "Qt5::Core")


    add_library(Qt5::Sensors SHARED IMPORTED)


    set_property(TARGET Qt5::Sensors PROPERTY
      INTERFACE_INCLUDE_DIRECTORIES ${_Qt5Sensors_OWN_INCLUDE_DIRS})
    set_property(TARGET Qt5::Sensors PROPERTY
      INTERFACE_COMPILE_DEFINITIONS QT_SENSORS_LIB)

    set_property(TARGET Qt5::Sensors PROPERTY INTERFACE_QT_ENABLED_FEATURES )
    set_property(TARGET Qt5::Sensors PROPERTY INTERFACE_QT_DISABLED_FEATURES )

    # Qt 6 forward compatible properties.
    set_property(TARGET Qt5::Sensors
                 PROPERTY QT_ENABLED_PUBLIC_FEATURES
                 )
    set_property(TARGET Qt5::Sensors
                 PROPERTY QT_DISABLED_PUBLIC_FEATURES
                 )
    set_property(TARGET Qt5::Sensors
                 PROPERTY QT_ENABLED_PRIVATE_FEATURES
                 )
    set_property(TARGET Qt5::Sensors
                 PROPERTY QT_DISABLED_PRIVATE_FEATURES
                 sensorfw)

    set_property(TARGET Qt5::Sensors PROPERTY INTERFACE_QT_PLUGIN_TYPES "sensors;sensorgestures")

    set(_Qt5Sensors_PRIVATE_DIRS_EXIST TRUE)
    foreach (_Qt5Sensors_PRIVATE_DIR ${Qt5Sensors_OWN_PRIVATE_INCLUDE_DIRS})
        if (NOT EXISTS ${_Qt5Sensors_PRIVATE_DIR})
            set(_Qt5Sensors_PRIVATE_DIRS_EXIST FALSE)
        endif()
    endforeach()

    if (_Qt5Sensors_PRIVATE_DIRS_EXIST)
        add_library(Qt5::SensorsPrivate INTERFACE IMPORTED)
        set_property(TARGET Qt5::SensorsPrivate PROPERTY
            INTERFACE_INCLUDE_DIRECTORIES ${Qt5Sensors_OWN_PRIVATE_INCLUDE_DIRS}
        )
        set(_Qt5Sensors_PRIVATEDEPS)
        foreach(dep ${_Qt5Sensors_LIB_DEPENDENCIES})
            if (TARGET ${dep}Private)
                list(APPEND _Qt5Sensors_PRIVATEDEPS ${dep}Private)
            endif()
        endforeach()
        set_property(TARGET Qt5::SensorsPrivate PROPERTY
            INTERFACE_LINK_LIBRARIES Qt5::Sensors ${_Qt5Sensors_PRIVATEDEPS}
        )

        # Add a versionless target, for compatibility with Qt6.
        if(NOT "${QT_NO_CREATE_VERSIONLESS_TARGETS}" AND NOT TARGET Qt::SensorsPrivate)
            add_library(Qt::SensorsPrivate INTERFACE IMPORTED)
            set_target_properties(Qt::SensorsPrivate PROPERTIES
                INTERFACE_LINK_LIBRARIES "Qt5::SensorsPrivate"
            )
        endif()
    endif()

    _populate_Sensors_target_properties(RELEASE "libQt5Sensors.so.5.15.8" "" FALSE)




    # In Qt 5.15 the glob pattern was relaxed to also catch plugins not literally named Plugin.
    # Define QT5_STRICT_PLUGIN_GLOB or ModuleName_STRICT_PLUGIN_GLOB to revert to old behavior.
    if (QT5_STRICT_PLUGIN_GLOB OR Qt5Sensors_STRICT_PLUGIN_GLOB)
        file(GLOB pluginTargets "${CMAKE_CURRENT_LIST_DIR}/Qt5Sensors_*Plugin.cmake")
    else()
        file(GLOB pluginTargets "${CMAKE_CURRENT_LIST_DIR}/Qt5Sensors_*.cmake")
    endif()

    macro(_populate_Sensors_plugin_properties Plugin Configuration PLUGIN_LOCATION
          IsDebugAndRelease)
        set_property(TARGET Qt5::${Plugin} APPEND PROPERTY IMPORTED_CONFIGURATIONS ${Configuration})

        set(imported_location "${_qt5Sensors_install_prefix}/plugins/${PLUGIN_LOCATION}")
        _qt5_Sensors_check_file_exists(${imported_location})
        set_target_properties(Qt5::${Plugin} PROPERTIES
            "IMPORTED_LOCATION_${Configuration}" ${imported_location}
        )

    endmacro()

    if (pluginTargets)
        foreach(pluginTarget ${pluginTargets})
            include(${pluginTarget})
        endforeach()
    endif()



    _qt5_Sensors_check_file_exists("${CMAKE_CURRENT_LIST_DIR}/Qt5SensorsConfigVersion.cmake")
endif()

# Add a versionless target, for compatibility with Qt6.
if(NOT "${QT_NO_CREATE_VERSIONLESS_TARGETS}" AND TARGET Qt5::Sensors AND NOT TARGET Qt::Sensors)
    add_library(Qt::Sensors INTERFACE IMPORTED)
    set_target_properties(Qt::Sensors PROPERTIES
        INTERFACE_LINK_LIBRARIES "Qt5::Sensors"
    )
endif()
