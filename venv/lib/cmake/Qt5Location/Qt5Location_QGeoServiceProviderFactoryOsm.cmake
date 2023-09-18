
add_library(Qt5::QGeoServiceProviderFactoryOsm MODULE IMPORTED)


_populate_Location_plugin_properties(QGeoServiceProviderFactoryOsm RELEASE "geoservices/libqtgeoservices_osm.so" FALSE)

list(APPEND Qt5Location_PLUGINS Qt5::QGeoServiceProviderFactoryOsm)
set_property(TARGET Qt5::Location APPEND PROPERTY QT_ALL_PLUGINS_geoservices Qt5::QGeoServiceProviderFactoryOsm)
set_property(TARGET Qt5::QGeoServiceProviderFactoryOsm PROPERTY QT_PLUGIN_TYPE "geoservices")
set_property(TARGET Qt5::QGeoServiceProviderFactoryOsm PROPERTY QT_PLUGIN_EXTENDS "")
set_property(TARGET Qt5::QGeoServiceProviderFactoryOsm PROPERTY QT_PLUGIN_CLASS_NAME "QGeoServiceProviderFactoryOsm")
