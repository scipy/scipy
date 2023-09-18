
add_library(Qt5::QEvdevGamepadBackendPlugin MODULE IMPORTED)


_populate_Gamepad_plugin_properties(QEvdevGamepadBackendPlugin RELEASE "gamepads/libevdevgamepad.so" FALSE)

list(APPEND Qt5Gamepad_PLUGINS Qt5::QEvdevGamepadBackendPlugin)
set_property(TARGET Qt5::Gamepad APPEND PROPERTY QT_ALL_PLUGINS_gamepads Qt5::QEvdevGamepadBackendPlugin)
set_property(TARGET Qt5::QEvdevGamepadBackendPlugin PROPERTY QT_PLUGIN_TYPE "gamepads")
set_property(TARGET Qt5::QEvdevGamepadBackendPlugin PROPERTY QT_PLUGIN_EXTENDS "")
set_property(TARGET Qt5::QEvdevGamepadBackendPlugin PROPERTY QT_PLUGIN_CLASS_NAME "QEvdevGamepadBackendPlugin")
