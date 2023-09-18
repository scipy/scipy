
add_library(Qt5::QtVirtualKeyboardHangulPlugin MODULE IMPORTED)


_populate_VirtualKeyboard_plugin_properties(QtVirtualKeyboardHangulPlugin RELEASE "virtualkeyboard/libqtvirtualkeyboard_hangul.so" FALSE)

list(APPEND Qt5VirtualKeyboard_PLUGINS Qt5::QtVirtualKeyboardHangulPlugin)
set_property(TARGET Qt5::VirtualKeyboard APPEND PROPERTY QT_ALL_PLUGINS_virtualkeyboard Qt5::QtVirtualKeyboardHangulPlugin)
set_property(TARGET Qt5::QtVirtualKeyboardHangulPlugin PROPERTY QT_PLUGIN_TYPE "virtualkeyboard")
set_property(TARGET Qt5::QtVirtualKeyboardHangulPlugin PROPERTY QT_PLUGIN_EXTENDS "")
set_property(TARGET Qt5::QtVirtualKeyboardHangulPlugin PROPERTY QT_PLUGIN_CLASS_NAME "QtVirtualKeyboardHangulPlugin")
