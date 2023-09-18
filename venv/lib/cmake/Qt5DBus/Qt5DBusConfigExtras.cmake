
if (NOT TARGET Qt5::qdbuscpp2xml)
    add_executable(Qt5::qdbuscpp2xml IMPORTED)

    set(imported_location "${_qt5DBus_install_prefix}/bin/qdbuscpp2xml")
    _qt5_DBus_check_file_exists(${imported_location})

    set_target_properties(Qt5::qdbuscpp2xml PROPERTIES
        IMPORTED_LOCATION ${imported_location}
    )
endif()

if (NOT TARGET Qt5::qdbusxml2cpp)
    add_executable(Qt5::qdbusxml2cpp IMPORTED)

    set(imported_location "${_qt5DBus_install_prefix}/bin/qdbusxml2cpp")
    _qt5_DBus_check_file_exists(${imported_location})

    set_target_properties(Qt5::qdbusxml2cpp PROPERTIES
        IMPORTED_LOCATION ${imported_location}
    )
endif()

set(Qt5DBus_QDBUSCPP2XML_EXECUTABLE Qt5::qdbuscpp2xml)
set(Qt5DBus_QDBUSXML2CPP_EXECUTABLE Qt5::qdbusxml2cpp)

# Create versionless tool targets.
foreach(__qt_tool qdbuscpp2xml qdbusxml2cpp)
    if(NOT "${QT_NO_CREATE_VERSIONLESS_TARGETS}" AND NOT TARGET Qt::${__qt_tool}
       AND TARGET Qt5::${__qt_tool})
        add_executable(Qt::${__qt_tool} IMPORTED)
        get_target_property(__qt_imported_location Qt5::${__qt_tool} IMPORTED_LOCATION)
        set_target_properties(Qt::${__qt_tool}
                              PROPERTIES IMPORTED_LOCATION "${__qt_imported_location}")
    endif()
endforeach()
