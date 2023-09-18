if (CMAKE_VERSION VERSION_LESS 3.1.0)
    message(FATAL_ERROR "Qt 5 3DQuickExtras module requires at least CMake version 3.1.0")
endif()

get_filename_component(_qt53DQuickExtras_install_prefix "${CMAKE_CURRENT_LIST_DIR}/../../../" ABSOLUTE)

# For backwards compatibility only. Use Qt53DQuickExtras_VERSION instead.
set(Qt53DQuickExtras_VERSION_STRING 5.15.8)

set(Qt53DQuickExtras_LIBRARIES Qt5::3DQuickExtras)

macro(_qt5_3DQuickExtras_check_file_exists file)
    if(NOT EXISTS "${file}" )
        message(FATAL_ERROR "The imported target \"Qt5::3DQuickExtras\" references the file
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


macro(_populate_3DQuickExtras_target_properties Configuration LIB_LOCATION IMPLIB_LOCATION
      IsDebugAndRelease)
    set_property(TARGET Qt5::3DQuickExtras APPEND PROPERTY IMPORTED_CONFIGURATIONS ${Configuration})

    set(imported_location "${_qt53DQuickExtras_install_prefix}/lib/${LIB_LOCATION}")
    _qt5_3DQuickExtras_check_file_exists(${imported_location})
    set(_deps
        ${_Qt53DQuickExtras_LIB_DEPENDENCIES}
    )
    set(_static_deps
    )

    set_target_properties(Qt5::3DQuickExtras PROPERTIES
        "IMPORTED_LOCATION_${Configuration}" ${imported_location}
        "IMPORTED_SONAME_${Configuration}" "libQt53DQuickExtras.so.5"
        # For backward compatibility with CMake < 2.8.12
        "IMPORTED_LINK_INTERFACE_LIBRARIES_${Configuration}" "${_deps};${_static_deps}"
    )
    set_property(TARGET Qt5::3DQuickExtras APPEND PROPERTY INTERFACE_LINK_LIBRARIES
                 "${_deps}"
    )


endmacro()

if (NOT TARGET Qt5::3DQuickExtras)

    set(_Qt53DQuickExtras_OWN_INCLUDE_DIRS "${_qt53DQuickExtras_install_prefix}/include/qt/" "${_qt53DQuickExtras_install_prefix}/include/qt/Qt3DQuickExtras")
    set(Qt53DQuickExtras_PRIVATE_INCLUDE_DIRS
        "${_qt53DQuickExtras_install_prefix}/include/qt/Qt3DQuickExtras/5.15.8"
        "${_qt53DQuickExtras_install_prefix}/include/qt/Qt3DQuickExtras/5.15.8/Qt3DQuickExtras"
    )
    include("${CMAKE_CURRENT_LIST_DIR}/ExtraSourceIncludes.cmake" OPTIONAL)

    foreach(_dir ${_Qt53DQuickExtras_OWN_INCLUDE_DIRS})
        _qt5_3DQuickExtras_check_file_exists(${_dir})
    endforeach()

    # Only check existence of private includes if the Private component is
    # specified.
    list(FIND Qt53DQuickExtras_FIND_COMPONENTS Private _check_private)
    if (NOT _check_private STREQUAL -1)
        foreach(_dir ${Qt53DQuickExtras_PRIVATE_INCLUDE_DIRS})
            _qt5_3DQuickExtras_check_file_exists(${_dir})
        endforeach()
    endif()

    set(Qt53DQuickExtras_INCLUDE_DIRS ${_Qt53DQuickExtras_OWN_INCLUDE_DIRS})

    set(Qt53DQuickExtras_DEFINITIONS -DQT_3DQUICKEXTRAS_LIB)
    set(Qt53DQuickExtras_COMPILE_DEFINITIONS QT_3DQUICKEXTRAS_LIB)
    set(_Qt53DQuickExtras_MODULE_DEPENDENCIES "3DExtras;3DInput;3DQuick;3DRender;3DLogic;3DCore;Gui;Qml;Core")


    set(Qt53DQuickExtras_OWN_PRIVATE_INCLUDE_DIRS ${Qt53DQuickExtras_PRIVATE_INCLUDE_DIRS})

    set(_Qt53DQuickExtras_FIND_DEPENDENCIES_REQUIRED)
    if (Qt53DQuickExtras_FIND_REQUIRED)
        set(_Qt53DQuickExtras_FIND_DEPENDENCIES_REQUIRED REQUIRED)
    endif()
    set(_Qt53DQuickExtras_FIND_DEPENDENCIES_QUIET)
    if (Qt53DQuickExtras_FIND_QUIETLY)
        set(_Qt53DQuickExtras_DEPENDENCIES_FIND_QUIET QUIET)
    endif()
    set(_Qt53DQuickExtras_FIND_VERSION_EXACT)
    if (Qt53DQuickExtras_FIND_VERSION_EXACT)
        set(_Qt53DQuickExtras_FIND_VERSION_EXACT EXACT)
    endif()

    set(Qt53DQuickExtras_EXECUTABLE_COMPILE_FLAGS "")

    foreach(_module_dep ${_Qt53DQuickExtras_MODULE_DEPENDENCIES})
        if (NOT Qt5${_module_dep}_FOUND)
            find_package(Qt5${_module_dep}
                5.15.8 ${_Qt53DQuickExtras_FIND_VERSION_EXACT}
                ${_Qt53DQuickExtras_DEPENDENCIES_FIND_QUIET}
                ${_Qt53DQuickExtras_FIND_DEPENDENCIES_REQUIRED}
                PATHS "${CMAKE_CURRENT_LIST_DIR}/.." NO_DEFAULT_PATH
            )
        endif()

        if (NOT Qt5${_module_dep}_FOUND)
            set(Qt53DQuickExtras_FOUND False)
            return()
        endif()

        list(APPEND Qt53DQuickExtras_INCLUDE_DIRS "${Qt5${_module_dep}_INCLUDE_DIRS}")
        list(APPEND Qt53DQuickExtras_PRIVATE_INCLUDE_DIRS "${Qt5${_module_dep}_PRIVATE_INCLUDE_DIRS}")
        list(APPEND Qt53DQuickExtras_DEFINITIONS ${Qt5${_module_dep}_DEFINITIONS})
        list(APPEND Qt53DQuickExtras_COMPILE_DEFINITIONS ${Qt5${_module_dep}_COMPILE_DEFINITIONS})
        list(APPEND Qt53DQuickExtras_EXECUTABLE_COMPILE_FLAGS ${Qt5${_module_dep}_EXECUTABLE_COMPILE_FLAGS})
    endforeach()
    list(REMOVE_DUPLICATES Qt53DQuickExtras_INCLUDE_DIRS)
    list(REMOVE_DUPLICATES Qt53DQuickExtras_PRIVATE_INCLUDE_DIRS)
    list(REMOVE_DUPLICATES Qt53DQuickExtras_DEFINITIONS)
    list(REMOVE_DUPLICATES Qt53DQuickExtras_COMPILE_DEFINITIONS)
    list(REMOVE_DUPLICATES Qt53DQuickExtras_EXECUTABLE_COMPILE_FLAGS)

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
    if(TARGET Qt5::3DQuickExtras)
        return()
    endif()

    set(_Qt53DQuickExtras_LIB_DEPENDENCIES "Qt5::3DExtras;Qt5::3DInput;Qt5::3DQuick;Qt5::3DRender;Qt5::3DLogic;Qt5::3DCore;Qt5::Gui;Qt5::Qml;Qt5::Core")


    add_library(Qt5::3DQuickExtras SHARED IMPORTED)


    set_property(TARGET Qt5::3DQuickExtras PROPERTY
      INTERFACE_INCLUDE_DIRECTORIES ${_Qt53DQuickExtras_OWN_INCLUDE_DIRS})
    set_property(TARGET Qt5::3DQuickExtras PROPERTY
      INTERFACE_COMPILE_DEFINITIONS QT_3DQUICKEXTRAS_LIB)

    set_property(TARGET Qt5::3DQuickExtras PROPERTY INTERFACE_QT_ENABLED_FEATURES )
    set_property(TARGET Qt5::3DQuickExtras PROPERTY INTERFACE_QT_DISABLED_FEATURES )

    # Qt 6 forward compatible properties.
    set_property(TARGET Qt5::3DQuickExtras
                 PROPERTY QT_ENABLED_PUBLIC_FEATURES
                 )
    set_property(TARGET Qt5::3DQuickExtras
                 PROPERTY QT_DISABLED_PUBLIC_FEATURES
                 )
    set_property(TARGET Qt5::3DQuickExtras
                 PROPERTY QT_ENABLED_PRIVATE_FEATURES
                 )
    set_property(TARGET Qt5::3DQuickExtras
                 PROPERTY QT_DISABLED_PRIVATE_FEATURES
                 )

    set_property(TARGET Qt5::3DQuickExtras PROPERTY INTERFACE_QT_PLUGIN_TYPES "")

    set(_Qt53DQuickExtras_PRIVATE_DIRS_EXIST TRUE)
    foreach (_Qt53DQuickExtras_PRIVATE_DIR ${Qt53DQuickExtras_OWN_PRIVATE_INCLUDE_DIRS})
        if (NOT EXISTS ${_Qt53DQuickExtras_PRIVATE_DIR})
            set(_Qt53DQuickExtras_PRIVATE_DIRS_EXIST FALSE)
        endif()
    endforeach()

    if (_Qt53DQuickExtras_PRIVATE_DIRS_EXIST)
        add_library(Qt5::3DQuickExtrasPrivate INTERFACE IMPORTED)
        set_property(TARGET Qt5::3DQuickExtrasPrivate PROPERTY
            INTERFACE_INCLUDE_DIRECTORIES ${Qt53DQuickExtras_OWN_PRIVATE_INCLUDE_DIRS}
        )
        set(_Qt53DQuickExtras_PRIVATEDEPS)
        foreach(dep ${_Qt53DQuickExtras_LIB_DEPENDENCIES})
            if (TARGET ${dep}Private)
                list(APPEND _Qt53DQuickExtras_PRIVATEDEPS ${dep}Private)
            endif()
        endforeach()
        set_property(TARGET Qt5::3DQuickExtrasPrivate PROPERTY
            INTERFACE_LINK_LIBRARIES Qt5::3DQuickExtras ${_Qt53DQuickExtras_PRIVATEDEPS}
        )

        # Add a versionless target, for compatibility with Qt6.
        if(NOT "${QT_NO_CREATE_VERSIONLESS_TARGETS}" AND NOT TARGET Qt::3DQuickExtrasPrivate)
            add_library(Qt::3DQuickExtrasPrivate INTERFACE IMPORTED)
            set_target_properties(Qt::3DQuickExtrasPrivate PROPERTIES
                INTERFACE_LINK_LIBRARIES "Qt5::3DQuickExtrasPrivate"
            )
        endif()
    endif()

    _populate_3DQuickExtras_target_properties(RELEASE "libQt53DQuickExtras.so.5.15.8" "" FALSE)




    # In Qt 5.15 the glob pattern was relaxed to also catch plugins not literally named Plugin.
    # Define QT5_STRICT_PLUGIN_GLOB or ModuleName_STRICT_PLUGIN_GLOB to revert to old behavior.
    if (QT5_STRICT_PLUGIN_GLOB OR Qt53DQuickExtras_STRICT_PLUGIN_GLOB)
        file(GLOB pluginTargets "${CMAKE_CURRENT_LIST_DIR}/Qt53DQuickExtras_*Plugin.cmake")
    else()
        file(GLOB pluginTargets "${CMAKE_CURRENT_LIST_DIR}/Qt53DQuickExtras_*.cmake")
    endif()

    macro(_populate_3DQuickExtras_plugin_properties Plugin Configuration PLUGIN_LOCATION
          IsDebugAndRelease)
        set_property(TARGET Qt5::${Plugin} APPEND PROPERTY IMPORTED_CONFIGURATIONS ${Configuration})

        set(imported_location "${_qt53DQuickExtras_install_prefix}/plugins/${PLUGIN_LOCATION}")
        _qt5_3DQuickExtras_check_file_exists(${imported_location})
        set_target_properties(Qt5::${Plugin} PROPERTIES
            "IMPORTED_LOCATION_${Configuration}" ${imported_location}
        )

    endmacro()

    if (pluginTargets)
        foreach(pluginTarget ${pluginTargets})
            include(${pluginTarget})
        endforeach()
    endif()



    _qt5_3DQuickExtras_check_file_exists("${CMAKE_CURRENT_LIST_DIR}/Qt53DQuickExtrasConfigVersion.cmake")
endif()

# Add a versionless target, for compatibility with Qt6.
if(NOT "${QT_NO_CREATE_VERSIONLESS_TARGETS}" AND TARGET Qt5::3DQuickExtras AND NOT TARGET Qt::3DQuickExtras)
    add_library(Qt::3DQuickExtras INTERFACE IMPORTED)
    set_target_properties(Qt::3DQuickExtras PROPERTIES
        INTERFACE_LINK_LIBRARIES "Qt5::3DQuickExtras"
    )
endif()
