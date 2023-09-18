
add_library(Qt5::VirtualCanBusPlugin MODULE IMPORTED)


_populate_SerialBus_plugin_properties(VirtualCanBusPlugin RELEASE "canbus/libqtvirtualcanbus.so" FALSE)

list(APPEND Qt5SerialBus_PLUGINS Qt5::VirtualCanBusPlugin)
set_property(TARGET Qt5::SerialBus APPEND PROPERTY QT_ALL_PLUGINS_canbus Qt5::VirtualCanBusPlugin)
set_property(TARGET Qt5::VirtualCanBusPlugin PROPERTY QT_PLUGIN_TYPE "canbus")
set_property(TARGET Qt5::VirtualCanBusPlugin PROPERTY QT_PLUGIN_EXTENDS "")
set_property(TARGET Qt5::VirtualCanBusPlugin PROPERTY QT_PLUGIN_CLASS_NAME "VirtualCanBusPlugin")
