if (CMAKE_VERSION VERSION_LESS 3.1.0)
    message(FATAL_ERROR "Qt 5 QuickShapes module requires at least CMake version 3.1.0")
endif()

get_filename_component(_qt5QuickShapes_install_prefix "${CMAKE_CURRENT_LIST_DIR}/../../../" ABSOLUTE)

# For backwards compatibility only. Use Qt5QuickShapes_VERSION instead.
set(Qt5QuickShapes_VERSION_STRING 5.15.8)

set(Qt5QuickShapes_LIBRARIES Qt5::QuickShapes)

macro(_qt5_QuickShapes_check_file_exists file)
    if(NOT EXISTS "${file}" )
        message(FATAL_ERROR "The imported target \"Qt5::QuickShapes\" references the file
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


macro(_populate_QuickShapes_target_properties Configuration LIB_LOCATION IMPLIB_LOCATION
      IsDebugAndRelease)
    set_property(TARGET Qt5::QuickShapes APPEND PROPERTY IMPORTED_CONFIGURATIONS ${Configuration})

    set(imported_location "${_qt5QuickShapes_install_prefix}/lib/${LIB_LOCATION}")
    _qt5_QuickShapes_check_file_exists(${imported_location})
    set(_deps
        ${_Qt5QuickShapes_LIB_DEPENDENCIES}
    )
    set(_static_deps
    )

    set_target_properties(Qt5::QuickShapes PROPERTIES
        "IMPORTED_LOCATION_${Configuration}" ${imported_location}
        "IMPORTED_SONAME_${Configuration}" "libQt5QuickShapes.so.5"
        # For backward compatibility with CMake < 2.8.12
        "IMPORTED_LINK_INTERFACE_LIBRARIES_${Configuration}" "${_deps};${_static_deps}"
    )
    set_property(TARGET Qt5::QuickShapes APPEND PROPERTY INTERFACE_LINK_LIBRARIES
                 "${_deps}"
    )


endmacro()

if (NOT TARGET Qt5::QuickShapes)

    set(_Qt5QuickShapes_OWN_INCLUDE_DIRS "${_qt5QuickShapes_install_prefix}/include/qt/" "${_qt5QuickShapes_install_prefix}/include/qt/QtQuickShapes")
    set(Qt5QuickShapes_PRIVATE_INCLUDE_DIRS
        "${_qt5QuickShapes_install_prefix}/include/qt/QtQuickShapes/5.15.8"
        "${_qt5QuickShapes_install_prefix}/include/qt/QtQuickShapes/5.15.8/QtQuickShapes"
    )
    include("${CMAKE_CURRENT_LIST_DIR}/ExtraSourceIncludes.cmake" OPTIONAL)

    foreach(_dir ${_Qt5QuickShapes_OWN_INCLUDE_DIRS})
        _qt5_QuickShapes_check_file_exists(${_dir})
    endforeach()

    # Only check existence of private includes if the Private component is
    # specified.
    list(FIND Qt5QuickShapes_FIND_COMPONENTS Private _check_private)
    if (NOT _check_private STREQUAL -1)
        foreach(_dir ${Qt5QuickShapes_PRIVATE_INCLUDE_DIRS})
            _qt5_QuickShapes_check_file_exists(${_dir})
        endforeach()
    endif()

    set(Qt5QuickShapes_INCLUDE_DIRS ${_Qt5QuickShapes_OWN_INCLUDE_DIRS})

    set(Qt5QuickShapes_DEFINITIONS -DQT_QUICKSHAPES_LIB)
    set(Qt5QuickShapes_COMPILE_DEFINITIONS QT_QUICKSHAPES_LIB)
    set(_Qt5QuickShapes_MODULE_DEPENDENCIES "Quick;Quick;Quick;Gui;Qml;Core")


    set(Qt5QuickShapes_OWN_PRIVATE_INCLUDE_DIRS ${Qt5QuickShapes_PRIVATE_INCLUDE_DIRS})

    set(_Qt5QuickShapes_FIND_DEPENDENCIES_REQUIRED)
    if (Qt5QuickShapes_FIND_REQUIRED)
        set(_Qt5QuickShapes_FIND_DEPENDENCIES_REQUIRED REQUIRED)
    endif()
    set(_Qt5QuickShapes_FIND_DEPENDENCIES_QUIET)
    if (Qt5QuickShapes_FIND_QUIETLY)
        set(_Qt5QuickShapes_DEPENDENCIES_FIND_QUIET QUIET)
    endif()
    set(_Qt5QuickShapes_FIND_VERSION_EXACT)
    if (Qt5QuickShapes_FIND_VERSION_EXACT)
        set(_Qt5QuickShapes_FIND_VERSION_EXACT EXACT)
    endif()

    set(Qt5QuickShapes_EXECUTABLE_COMPILE_FLAGS "")

    foreach(_module_dep ${_Qt5QuickShapes_MODULE_DEPENDENCIES})
        if (NOT Qt5${_module_dep}_FOUND)
            find_package(Qt5${_module_dep}
                5.15.8 ${_Qt5QuickShapes_FIND_VERSION_EXACT}
                ${_Qt5QuickShapes_DEPENDENCIES_FIND_QUIET}
                ${_Qt5QuickShapes_FIND_DEPENDENCIES_REQUIRED}
                PATHS "${CMAKE_CURRENT_LIST_DIR}/.." NO_DEFAULT_PATH
            )
        endif()

        if (NOT Qt5${_module_dep}_FOUND)
            set(Qt5QuickShapes_FOUND False)
            return()
        endif()

        list(APPEND Qt5QuickShapes_INCLUDE_DIRS "${Qt5${_module_dep}_INCLUDE_DIRS}")
        list(APPEND Qt5QuickShapes_PRIVATE_INCLUDE_DIRS "${Qt5${_module_dep}_PRIVATE_INCLUDE_DIRS}")
        list(APPEND Qt5QuickShapes_DEFINITIONS ${Qt5${_module_dep}_DEFINITIONS})
        list(APPEND Qt5QuickShapes_COMPILE_DEFINITIONS ${Qt5${_module_dep}_COMPILE_DEFINITIONS})
        list(APPEND Qt5QuickShapes_EXECUTABLE_COMPILE_FLAGS ${Qt5${_module_dep}_EXECUTABLE_COMPILE_FLAGS})
    endforeach()
    list(REMOVE_DUPLICATES Qt5QuickShapes_INCLUDE_DIRS)
    list(REMOVE_DUPLICATES Qt5QuickShapes_PRIVATE_INCLUDE_DIRS)
    list(REMOVE_DUPLICATES Qt5QuickShapes_DEFINITIONS)
    list(REMOVE_DUPLICATES Qt5QuickShapes_COMPILE_DEFINITIONS)
    list(REMOVE_DUPLICATES Qt5QuickShapes_EXECUTABLE_COMPILE_FLAGS)

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
    if(TARGET Qt5::QuickShapes)
        return()
    endif()

    set(_Qt5QuickShapes_LIB_DEPENDENCIES "Qt5::Quick;Qt5::Quick;Qt5::Quick;Qt5::Gui;Qt5::Qml;Qt5::Core")


    add_library(Qt5::QuickShapes SHARED IMPORTED)


    set_property(TARGET Qt5::QuickShapes PROPERTY
      INTERFACE_INCLUDE_DIRECTORIES ${_Qt5QuickShapes_OWN_INCLUDE_DIRS})
    set_property(TARGET Qt5::QuickShapes PROPERTY
      INTERFACE_COMPILE_DEFINITIONS QT_QUICKSHAPES_LIB)

    set_property(TARGET Qt5::QuickShapes PROPERTY INTERFACE_QT_ENABLED_FEATURES )
    set_property(TARGET Qt5::QuickShapes PROPERTY INTERFACE_QT_DISABLED_FEATURES )

    # Qt 6 forward compatible properties.
    set_property(TARGET Qt5::QuickShapes
                 PROPERTY QT_ENABLED_PUBLIC_FEATURES
                 )
    set_property(TARGET Qt5::QuickShapes
                 PROPERTY QT_DISABLED_PUBLIC_FEATURES
                 )
    set_property(TARGET Qt5::QuickShapes
                 PROPERTY QT_ENABLED_PRIVATE_FEATURES
                 )
    set_property(TARGET Qt5::QuickShapes
                 PROPERTY QT_DISABLED_PRIVATE_FEATURES
                 )

    set_property(TARGET Qt5::QuickShapes PROPERTY INTERFACE_QT_PLUGIN_TYPES "")

    set(_Qt5QuickShapes_PRIVATE_DIRS_EXIST TRUE)
    foreach (_Qt5QuickShapes_PRIVATE_DIR ${Qt5QuickShapes_OWN_PRIVATE_INCLUDE_DIRS})
        if (NOT EXISTS ${_Qt5QuickShapes_PRIVATE_DIR})
            set(_Qt5QuickShapes_PRIVATE_DIRS_EXIST FALSE)
        endif()
    endforeach()

    if (_Qt5QuickShapes_PRIVATE_DIRS_EXIST)
        add_library(Qt5::QuickShapesPrivate INTERFACE IMPORTED)
        set_property(TARGET Qt5::QuickShapesPrivate PROPERTY
            INTERFACE_INCLUDE_DIRECTORIES ${Qt5QuickShapes_OWN_PRIVATE_INCLUDE_DIRS}
        )
        set(_Qt5QuickShapes_PRIVATEDEPS)
        foreach(dep ${_Qt5QuickShapes_LIB_DEPENDENCIES})
            if (TARGET ${dep}Private)
                list(APPEND _Qt5QuickShapes_PRIVATEDEPS ${dep}Private)
            endif()
        endforeach()
        set_property(TARGET Qt5::QuickShapesPrivate PROPERTY
            INTERFACE_LINK_LIBRARIES Qt5::QuickShapes ${_Qt5QuickShapes_PRIVATEDEPS}
        )

        # Add a versionless target, for compatibility with Qt6.
        if(NOT "${QT_NO_CREATE_VERSIONLESS_TARGETS}" AND NOT TARGET Qt::QuickShapesPrivate)
            add_library(Qt::QuickShapesPrivate INTERFACE IMPORTED)
            set_target_properties(Qt::QuickShapesPrivate PROPERTIES
                INTERFACE_LINK_LIBRARIES "Qt5::QuickShapesPrivate"
            )
        endif()
    endif()

    _populate_QuickShapes_target_properties(RELEASE "libQt5QuickShapes.so.5.15.8" "" FALSE)







    _qt5_QuickShapes_check_file_exists("${CMAKE_CURRENT_LIST_DIR}/Qt5QuickShapesConfigVersion.cmake")
endif()

# Add a versionless target, for compatibility with Qt6.
if(NOT "${QT_NO_CREATE_VERSIONLESS_TARGETS}" AND TARGET Qt5::QuickShapes AND NOT TARGET Qt::QuickShapes)
    add_library(Qt::QuickShapes INTERFACE IMPORTED)
    set_target_properties(Qt::QuickShapes PROPERTIES
        INTERFACE_LINK_LIBRARIES "Qt5::QuickShapes"
    )
endif()
