
add_library(Qt5::QGeoServiceProviderFactoryMapboxGL MODULE IMPORTED)


_populate_Location_plugin_properties(QGeoServiceProviderFactoryMapboxGL RELEASE "geoservices/libqtgeoservices_mapboxgl.so" FALSE)

list(APPEND Qt5Location_PLUGINS Qt5::QGeoServiceProviderFactoryMapboxGL)
set_property(TARGET Qt5::Location APPEND PROPERTY QT_ALL_PLUGINS_geoservices Qt5::QGeoServiceProviderFactoryMapboxGL)
set_property(TARGET Qt5::QGeoServiceProviderFactoryMapboxGL PROPERTY QT_PLUGIN_TYPE "geoservices")
set_property(TARGET Qt5::QGeoServiceProviderFactoryMapboxGL PROPERTY QT_PLUGIN_EXTENDS "")
set_property(TARGET Qt5::QGeoServiceProviderFactoryMapboxGL PROPERTY QT_PLUGIN_CLASS_NAME "QGeoServiceProviderFactoryMapboxGL")
