
add_library(Qt5::QShakeSensorGesturePlugin MODULE IMPORTED)


_populate_Sensors_plugin_properties(QShakeSensorGesturePlugin RELEASE "sensorgestures/libqtsensorgestures_shakeplugin.so" FALSE)

list(APPEND Qt5Sensors_PLUGINS Qt5::QShakeSensorGesturePlugin)
set_property(TARGET Qt5::Sensors APPEND PROPERTY QT_ALL_PLUGINS_sensorgestures Qt5::QShakeSensorGesturePlugin)
set_property(TARGET Qt5::QShakeSensorGesturePlugin PROPERTY QT_PLUGIN_TYPE "sensorgestures")
set_property(TARGET Qt5::QShakeSensorGesturePlugin PROPERTY QT_PLUGIN_EXTENDS "-")
set_property(TARGET Qt5::QShakeSensorGesturePlugin PROPERTY QT_PLUGIN_CLASS_NAME "QShakeSensorGesturePlugin")
