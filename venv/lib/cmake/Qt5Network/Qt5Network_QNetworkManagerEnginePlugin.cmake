
add_library(Qt5::QNetworkManagerEnginePlugin MODULE IMPORTED)


_populate_Network_plugin_properties(QNetworkManagerEnginePlugin RELEASE "bearer/libqnmbearer.so" FALSE)

list(APPEND Qt5Network_PLUGINS Qt5::QNetworkManagerEnginePlugin)
set_property(TARGET Qt5::Network APPEND PROPERTY QT_ALL_PLUGINS_bearer Qt5::QNetworkManagerEnginePlugin)
set_property(TARGET Qt5::QNetworkManagerEnginePlugin PROPERTY QT_PLUGIN_TYPE "bearer")
set_property(TARGET Qt5::QNetworkManagerEnginePlugin PROPERTY QT_PLUGIN_EXTENDS "")
set_property(TARGET Qt5::QNetworkManagerEnginePlugin PROPERTY QT_PLUGIN_CLASS_NAME "QNetworkManagerEnginePlugin")
