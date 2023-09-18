
add_library(Qt5::LinuxSensorPlugin MODULE IMPORTED)


_populate_Sensors_plugin_properties(LinuxSensorPlugin RELEASE "sensors/libqtsensors_linuxsys.so" FALSE)

list(APPEND Qt5Sensors_PLUGINS Qt5::LinuxSensorPlugin)
set_property(TARGET Qt5::Sensors APPEND PROPERTY QT_ALL_PLUGINS_sensors Qt5::LinuxSensorPlugin)
set_property(TARGET Qt5::LinuxSensorPlugin PROPERTY QT_PLUGIN_TYPE "sensors")
set_property(TARGET Qt5::LinuxSensorPlugin PROPERTY QT_PLUGIN_EXTENDS "")
set_property(TARGET Qt5::LinuxSensorPlugin PROPERTY QT_PLUGIN_CLASS_NAME "LinuxSensorPlugin")
