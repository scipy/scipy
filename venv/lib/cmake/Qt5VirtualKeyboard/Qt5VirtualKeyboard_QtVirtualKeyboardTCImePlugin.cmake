
add_library(Qt5::QtVirtualKeyboardTCImePlugin MODULE IMPORTED)


_populate_VirtualKeyboard_plugin_properties(QtVirtualKeyboardTCImePlugin RELEASE "virtualkeyboard/libqtvirtualkeyboard_tcime.so" FALSE)

list(APPEND Qt5VirtualKeyboard_PLUGINS Qt5::QtVirtualKeyboardTCImePlugin)
set_property(TARGET Qt5::VirtualKeyboard APPEND PROPERTY QT_ALL_PLUGINS_virtualkeyboard Qt5::QtVirtualKeyboardTCImePlugin)
set_property(TARGET Qt5::QtVirtualKeyboardTCImePlugin PROPERTY QT_PLUGIN_TYPE "virtualkeyboard")
set_property(TARGET Qt5::QtVirtualKeyboardTCImePlugin PROPERTY QT_PLUGIN_EXTENDS "")
set_property(TARGET Qt5::QtVirtualKeyboardTCImePlugin PROPERTY QT_PLUGIN_CLASS_NAME "QtVirtualKeyboardTCImePlugin")
