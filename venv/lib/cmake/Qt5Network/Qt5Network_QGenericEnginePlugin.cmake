
add_library(Qt5::QGenericEnginePlugin MODULE IMPORTED)


_populate_Network_plugin_properties(QGenericEnginePlugin RELEASE "bearer/libqgenericbearer.so" FALSE)

list(APPEND Qt5Network_PLUGINS Qt5::QGenericEnginePlugin)
set_property(TARGET Qt5::Network APPEND PROPERTY QT_ALL_PLUGINS_bearer Qt5::QGenericEnginePlugin)
set_property(TARGET Qt5::QGenericEnginePlugin PROPERTY QT_PLUGIN_TYPE "bearer")
set_property(TARGET Qt5::QGenericEnginePlugin PROPERTY QT_PLUGIN_EXTENDS "")
set_property(TARGET Qt5::QGenericEnginePlugin PROPERTY QT_PLUGIN_CLASS_NAME "QGenericEnginePlugin")
