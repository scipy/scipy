
add_library(Qt5::QGeoPositionInfoSourceFactoryGeoclue MODULE IMPORTED)


_populate_Positioning_plugin_properties(QGeoPositionInfoSourceFactoryGeoclue RELEASE "position/libqtposition_geoclue.so" FALSE)

list(APPEND Qt5Positioning_PLUGINS Qt5::QGeoPositionInfoSourceFactoryGeoclue)
set_property(TARGET Qt5::Positioning APPEND PROPERTY QT_ALL_PLUGINS_position Qt5::QGeoPositionInfoSourceFactoryGeoclue)
set_property(TARGET Qt5::QGeoPositionInfoSourceFactoryGeoclue PROPERTY QT_PLUGIN_TYPE "position")
set_property(TARGET Qt5::QGeoPositionInfoSourceFactoryGeoclue PROPERTY QT_PLUGIN_EXTENDS "")
set_property(TARGET Qt5::QGeoPositionInfoSourceFactoryGeoclue PROPERTY QT_PLUGIN_CLASS_NAME "QGeoPositionInfoSourceFactoryGeoclue")
