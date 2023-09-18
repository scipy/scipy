if (CMAKE_VERSION VERSION_LESS 3.1.0)
    message(FATAL_ERROR "Qt 5 Multimedia module requires at least CMake version 3.1.0")
endif()

get_filename_component(_qt5Multimedia_install_prefix "${CMAKE_CURRENT_LIST_DIR}/../../../" ABSOLUTE)

# For backwards compatibility only. Use Qt5Multimedia_VERSION instead.
set(Qt5Multimedia_VERSION_STRING 5.15.8)

set(Qt5Multimedia_LIBRARIES Qt5::Multimedia)

macro(_qt5_Multimedia_check_file_exists file)
    if(NOT EXISTS "${file}" )
        message(FATAL_ERROR "The imported target \"Qt5::Multimedia\" references the file
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


macro(_populate_Multimedia_target_properties Configuration LIB_LOCATION IMPLIB_LOCATION
      IsDebugAndRelease)
    set_property(TARGET Qt5::Multimedia APPEND PROPERTY IMPORTED_CONFIGURATIONS ${Configuration})

    set(imported_location "${_qt5Multimedia_install_prefix}/lib/${LIB_LOCATION}")
    _qt5_Multimedia_check_file_exists(${imported_location})
    set(_deps
        ${_Qt5Multimedia_LIB_DEPENDENCIES}
    )
    set(_static_deps
    )

    set_target_properties(Qt5::Multimedia PROPERTIES
        "IMPORTED_LOCATION_${Configuration}" ${imported_location}
        "IMPORTED_SONAME_${Configuration}" "libQt5Multimedia.so.5"
        # For backward compatibility with CMake < 2.8.12
        "IMPORTED_LINK_INTERFACE_LIBRARIES_${Configuration}" "${_deps};${_static_deps}"
    )
    set_property(TARGET Qt5::Multimedia APPEND PROPERTY INTERFACE_LINK_LIBRARIES
                 "${_deps}"
    )


endmacro()

if (NOT TARGET Qt5::Multimedia)

    set(_Qt5Multimedia_OWN_INCLUDE_DIRS "${_qt5Multimedia_install_prefix}/include/qt/" "${_qt5Multimedia_install_prefix}/include/qt/QtMultimedia")
    set(Qt5Multimedia_PRIVATE_INCLUDE_DIRS
        "${_qt5Multimedia_install_prefix}/include/qt/QtMultimedia/5.15.8"
        "${_qt5Multimedia_install_prefix}/include/qt/QtMultimedia/5.15.8/QtMultimedia"
    )
    include("${CMAKE_CURRENT_LIST_DIR}/ExtraSourceIncludes.cmake" OPTIONAL)

    foreach(_dir ${_Qt5Multimedia_OWN_INCLUDE_DIRS})
        _qt5_Multimedia_check_file_exists(${_dir})
    endforeach()

    # Only check existence of private includes if the Private component is
    # specified.
    list(FIND Qt5Multimedia_FIND_COMPONENTS Private _check_private)
    if (NOT _check_private STREQUAL -1)
        foreach(_dir ${Qt5Multimedia_PRIVATE_INCLUDE_DIRS})
            _qt5_Multimedia_check_file_exists(${_dir})
        endforeach()
    endif()

    set(Qt5Multimedia_INCLUDE_DIRS ${_Qt5Multimedia_OWN_INCLUDE_DIRS})

    set(Qt5Multimedia_DEFINITIONS -DQT_MULTIMEDIA_LIB)
    set(Qt5Multimedia_COMPILE_DEFINITIONS QT_MULTIMEDIA_LIB)
    set(_Qt5Multimedia_MODULE_DEPENDENCIES "Network;Gui;Core")


    set(Qt5Multimedia_OWN_PRIVATE_INCLUDE_DIRS ${Qt5Multimedia_PRIVATE_INCLUDE_DIRS})

    set(_Qt5Multimedia_FIND_DEPENDENCIES_REQUIRED)
    if (Qt5Multimedia_FIND_REQUIRED)
        set(_Qt5Multimedia_FIND_DEPENDENCIES_REQUIRED REQUIRED)
    endif()
    set(_Qt5Multimedia_FIND_DEPENDENCIES_QUIET)
    if (Qt5Multimedia_FIND_QUIETLY)
        set(_Qt5Multimedia_DEPENDENCIES_FIND_QUIET QUIET)
    endif()
    set(_Qt5Multimedia_FIND_VERSION_EXACT)
    if (Qt5Multimedia_FIND_VERSION_EXACT)
        set(_Qt5Multimedia_FIND_VERSION_EXACT EXACT)
    endif()

    set(Qt5Multimedia_EXECUTABLE_COMPILE_FLAGS "")

    foreach(_module_dep ${_Qt5Multimedia_MODULE_DEPENDENCIES})
        if (NOT Qt5${_module_dep}_FOUND)
            find_package(Qt5${_module_dep}
                5.15.8 ${_Qt5Multimedia_FIND_VERSION_EXACT}
                ${_Qt5Multimedia_DEPENDENCIES_FIND_QUIET}
                ${_Qt5Multimedia_FIND_DEPENDENCIES_REQUIRED}
                PATHS "${CMAKE_CURRENT_LIST_DIR}/.." NO_DEFAULT_PATH
            )
        endif()

        if (NOT Qt5${_module_dep}_FOUND)
            set(Qt5Multimedia_FOUND False)
            return()
        endif()

        list(APPEND Qt5Multimedia_INCLUDE_DIRS "${Qt5${_module_dep}_INCLUDE_DIRS}")
        list(APPEND Qt5Multimedia_PRIVATE_INCLUDE_DIRS "${Qt5${_module_dep}_PRIVATE_INCLUDE_DIRS}")
        list(APPEND Qt5Multimedia_DEFINITIONS ${Qt5${_module_dep}_DEFINITIONS})
        list(APPEND Qt5Multimedia_COMPILE_DEFINITIONS ${Qt5${_module_dep}_COMPILE_DEFINITIONS})
        list(APPEND Qt5Multimedia_EXECUTABLE_COMPILE_FLAGS ${Qt5${_module_dep}_EXECUTABLE_COMPILE_FLAGS})
    endforeach()
    list(REMOVE_DUPLICATES Qt5Multimedia_INCLUDE_DIRS)
    list(REMOVE_DUPLICATES Qt5Multimedia_PRIVATE_INCLUDE_DIRS)
    list(REMOVE_DUPLICATES Qt5Multimedia_DEFINITIONS)
    list(REMOVE_DUPLICATES Qt5Multimedia_COMPILE_DEFINITIONS)
    list(REMOVE_DUPLICATES Qt5Multimedia_EXECUTABLE_COMPILE_FLAGS)

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
    if(TARGET Qt5::Multimedia)
        return()
    endif()

    set(_Qt5Multimedia_LIB_DEPENDENCIES "Qt5::Network;Qt5::Gui;Qt5::Core")


    add_library(Qt5::Multimedia SHARED IMPORTED)


    set_property(TARGET Qt5::Multimedia PROPERTY
      INTERFACE_INCLUDE_DIRECTORIES ${_Qt5Multimedia_OWN_INCLUDE_DIRS})
    set_property(TARGET Qt5::Multimedia PROPERTY
      INTERFACE_COMPILE_DEFINITIONS QT_MULTIMEDIA_LIB)

    set_property(TARGET Qt5::Multimedia PROPERTY INTERFACE_QT_ENABLED_FEATURES )
    set_property(TARGET Qt5::Multimedia PROPERTY INTERFACE_QT_DISABLED_FEATURES )

    # Qt 6 forward compatible properties.
    set_property(TARGET Qt5::Multimedia
                 PROPERTY QT_ENABLED_PUBLIC_FEATURES
                 )
    set_property(TARGET Qt5::Multimedia
                 PROPERTY QT_DISABLED_PUBLIC_FEATURES
                 )
    set_property(TARGET Qt5::Multimedia
                 PROPERTY QT_ENABLED_PRIVATE_FEATURES
                 alsa;gstreamer_1_0;gstreamer;gstreamer_app;gstreamer_encodingprofiles;gstreamer_gl;linux_v4l;pulseaudio)
    set_property(TARGET Qt5::Multimedia
                 PROPERTY QT_DISABLED_PRIVATE_FEATURES
                 directshow;evr;gpu_vivante;gstreamer_0_10;gstreamer_photography;openal;resourcepolicy;wasapi;wmf;wmsdk;wshellitem)

    set_property(TARGET Qt5::Multimedia PROPERTY INTERFACE_QT_PLUGIN_TYPES "mediaservice;audio;video/bufferpool;video/gstvideorenderer;video/videonode;playlistformats;resourcepolicy")

    set(_Qt5Multimedia_PRIVATE_DIRS_EXIST TRUE)
    foreach (_Qt5Multimedia_PRIVATE_DIR ${Qt5Multimedia_OWN_PRIVATE_INCLUDE_DIRS})
        if (NOT EXISTS ${_Qt5Multimedia_PRIVATE_DIR})
            set(_Qt5Multimedia_PRIVATE_DIRS_EXIST FALSE)
        endif()
    endforeach()

    if (_Qt5Multimedia_PRIVATE_DIRS_EXIST)
        add_library(Qt5::MultimediaPrivate INTERFACE IMPORTED)
        set_property(TARGET Qt5::MultimediaPrivate PROPERTY
            INTERFACE_INCLUDE_DIRECTORIES ${Qt5Multimedia_OWN_PRIVATE_INCLUDE_DIRS}
        )
        set(_Qt5Multimedia_PRIVATEDEPS)
        foreach(dep ${_Qt5Multimedia_LIB_DEPENDENCIES})
            if (TARGET ${dep}Private)
                list(APPEND _Qt5Multimedia_PRIVATEDEPS ${dep}Private)
            endif()
        endforeach()
        set_property(TARGET Qt5::MultimediaPrivate PROPERTY
            INTERFACE_LINK_LIBRARIES Qt5::Multimedia ${_Qt5Multimedia_PRIVATEDEPS}
        )

        # Add a versionless target, for compatibility with Qt6.
        if(NOT "${QT_NO_CREATE_VERSIONLESS_TARGETS}" AND NOT TARGET Qt::MultimediaPrivate)
            add_library(Qt::MultimediaPrivate INTERFACE IMPORTED)
            set_target_properties(Qt::MultimediaPrivate PROPERTIES
                INTERFACE_LINK_LIBRARIES "Qt5::MultimediaPrivate"
            )
        endif()
    endif()

    _populate_Multimedia_target_properties(RELEASE "libQt5Multimedia.so.5.15.8" "" FALSE)




    # In Qt 5.15 the glob pattern was relaxed to also catch plugins not literally named Plugin.
    # Define QT5_STRICT_PLUGIN_GLOB or ModuleName_STRICT_PLUGIN_GLOB to revert to old behavior.
    if (QT5_STRICT_PLUGIN_GLOB OR Qt5Multimedia_STRICT_PLUGIN_GLOB)
        file(GLOB pluginTargets "${CMAKE_CURRENT_LIST_DIR}/Qt5Multimedia_*Plugin.cmake")
    else()
        file(GLOB pluginTargets "${CMAKE_CURRENT_LIST_DIR}/Qt5Multimedia_*.cmake")
    endif()

    macro(_populate_Multimedia_plugin_properties Plugin Configuration PLUGIN_LOCATION
          IsDebugAndRelease)
        set_property(TARGET Qt5::${Plugin} APPEND PROPERTY IMPORTED_CONFIGURATIONS ${Configuration})

        set(imported_location "${_qt5Multimedia_install_prefix}/plugins/${PLUGIN_LOCATION}")
        _qt5_Multimedia_check_file_exists(${imported_location})
        set_target_properties(Qt5::${Plugin} PROPERTIES
            "IMPORTED_LOCATION_${Configuration}" ${imported_location}
        )

    endmacro()

    if (pluginTargets)
        foreach(pluginTarget ${pluginTargets})
            include(${pluginTarget})
        endforeach()
    endif()



    _qt5_Multimedia_check_file_exists("${CMAKE_CURRENT_LIST_DIR}/Qt5MultimediaConfigVersion.cmake")
endif()

# Add a versionless target, for compatibility with Qt6.
if(NOT "${QT_NO_CREATE_VERSIONLESS_TARGETS}" AND TARGET Qt5::Multimedia AND NOT TARGET Qt::Multimedia)
    add_library(Qt::Multimedia INTERFACE IMPORTED)
    set_target_properties(Qt::Multimedia PROPERTIES
        INTERFACE_LINK_LIBRARIES "Qt5::Multimedia"
    )
endif()
