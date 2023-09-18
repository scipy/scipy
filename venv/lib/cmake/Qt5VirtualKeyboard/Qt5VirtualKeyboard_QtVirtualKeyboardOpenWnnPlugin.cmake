
add_library(Qt5::QtVirtualKeyboardOpenWnnPlugin MODULE IMPORTED)


_populate_VirtualKeyboard_plugin_properties(QtVirtualKeyboardOpenWnnPlugin RELEASE "virtualkeyboard/libqtvirtualkeyboard_openwnn.so" FALSE)

list(APPEND Qt5VirtualKeyboard_PLUGINS Qt5::QtVirtualKeyboardOpenWnnPlugin)
set_property(TARGET Qt5::VirtualKeyboard APPEND PROPERTY QT_ALL_PLUGINS_virtualkeyboard Qt5::QtVirtualKeyboardOpenWnnPlugin)
set_property(TARGET Qt5::QtVirtualKeyboardOpenWnnPlugin PROPERTY QT_PLUGIN_TYPE "virtualkeyboard")
set_property(TARGET Qt5::QtVirtualKeyboardOpenWnnPlugin PROPERTY QT_PLUGIN_EXTENDS "")
set_property(TARGET Qt5::QtVirtualKeyboardOpenWnnPlugin PROPERTY QT_PLUGIN_CLASS_NAME "QtVirtualKeyboardOpenWnnPlugin")
