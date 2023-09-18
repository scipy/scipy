
add_library(Qt5::QGeoServiceProviderFactoryMapbox MODULE IMPORTED)


_populate_Location_plugin_properties(QGeoServiceProviderFactoryMapbox RELEASE "geoservices/libqtgeoservices_mapbox.so" FALSE)

list(APPEND Qt5Location_PLUGINS Qt5::QGeoServiceProviderFactoryMapbox)
set_property(TARGET Qt5::Location APPEND PROPERTY QT_ALL_PLUGINS_geoservices Qt5::QGeoServiceProviderFactoryMapbox)
set_property(TARGET Qt5::QGeoServiceProviderFactoryMapbox PROPERTY QT_PLUGIN_TYPE "geoservices")
set_property(TARGET Qt5::QGeoServiceProviderFactoryMapbox PROPERTY QT_PLUGIN_EXTENDS "")
set_property(TARGET Qt5::QGeoServiceProviderFactoryMapbox PROPERTY QT_PLUGIN_CLASS_NAME "QGeoServiceProviderFactoryMapbox")
