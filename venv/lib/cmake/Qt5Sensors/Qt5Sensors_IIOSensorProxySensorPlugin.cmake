
add_library(Qt5::IIOSensorProxySensorPlugin MODULE IMPORTED)


_populate_Sensors_plugin_properties(IIOSensorProxySensorPlugin RELEASE "sensors/libqtsensors_iio-sensor-proxy.so" FALSE)

list(APPEND Qt5Sensors_PLUGINS Qt5::IIOSensorProxySensorPlugin)
set_property(TARGET Qt5::Sensors APPEND PROPERTY QT_ALL_PLUGINS_sensors Qt5::IIOSensorProxySensorPlugin)
set_property(TARGET Qt5::IIOSensorProxySensorPlugin PROPERTY QT_PLUGIN_TYPE "sensors")
set_property(TARGET Qt5::IIOSensorProxySensorPlugin PROPERTY QT_PLUGIN_EXTENDS "")
set_property(TARGET Qt5::IIOSensorProxySensorPlugin PROPERTY QT_PLUGIN_CLASS_NAME "IIOSensorProxySensorPlugin")
