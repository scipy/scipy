
add_library(Qt5::TinyCanBusPlugin MODULE IMPORTED)


_populate_SerialBus_plugin_properties(TinyCanBusPlugin RELEASE "canbus/libqttinycanbus.so" FALSE)

list(APPEND Qt5SerialBus_PLUGINS Qt5::TinyCanBusPlugin)
set_property(TARGET Qt5::SerialBus APPEND PROPERTY QT_ALL_PLUGINS_canbus Qt5::TinyCanBusPlugin)
set_property(TARGET Qt5::TinyCanBusPlugin PROPERTY QT_PLUGIN_TYPE "canbus")
set_property(TARGET Qt5::TinyCanBusPlugin PROPERTY QT_PLUGIN_EXTENDS "")
set_property(TARGET Qt5::TinyCanBusPlugin PROPERTY QT_PLUGIN_CLASS_NAME "TinyCanBusPlugin")
