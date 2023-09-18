
add_library(Qt5::QConnmanEnginePlugin MODULE IMPORTED)


_populate_Network_plugin_properties(QConnmanEnginePlugin RELEASE "bearer/libqconnmanbearer.so" FALSE)

list(APPEND Qt5Network_PLUGINS Qt5::QConnmanEnginePlugin)
set_property(TARGET Qt5::Network APPEND PROPERTY QT_ALL_PLUGINS_bearer Qt5::QConnmanEnginePlugin)
set_property(TARGET Qt5::QConnmanEnginePlugin PROPERTY QT_PLUGIN_TYPE "bearer")
set_property(TARGET Qt5::QConnmanEnginePlugin PROPERTY QT_PLUGIN_EXTENDS "")
set_property(TARGET Qt5::QConnmanEnginePlugin PROPERTY QT_PLUGIN_CLASS_NAME "QConnmanEnginePlugin")
