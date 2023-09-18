
add_library(Qt5::SocketCanBusPlugin MODULE IMPORTED)


_populate_SerialBus_plugin_properties(SocketCanBusPlugin RELEASE "canbus/libqtsocketcanbus.so" FALSE)

list(APPEND Qt5SerialBus_PLUGINS Qt5::SocketCanBusPlugin)
set_property(TARGET Qt5::SerialBus APPEND PROPERTY QT_ALL_PLUGINS_canbus Qt5::SocketCanBusPlugin)
set_property(TARGET Qt5::SocketCanBusPlugin PROPERTY QT_PLUGIN_TYPE "canbus")
set_property(TARGET Qt5::SocketCanBusPlugin PROPERTY QT_PLUGIN_EXTENDS "")
set_property(TARGET Qt5::SocketCanBusPlugin PROPERTY QT_PLUGIN_CLASS_NAME "SocketCanBusPlugin")
