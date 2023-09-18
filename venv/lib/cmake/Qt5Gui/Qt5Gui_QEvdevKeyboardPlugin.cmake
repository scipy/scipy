
add_library(Qt5::QEvdevKeyboardPlugin MODULE IMPORTED)


_populate_Gui_plugin_properties(QEvdevKeyboardPlugin RELEASE "generic/libqevdevkeyboardplugin.so" FALSE)

list(APPEND Qt5Gui_PLUGINS Qt5::QEvdevKeyboardPlugin)
set_property(TARGET Qt5::Gui APPEND PROPERTY QT_ALL_PLUGINS_generic Qt5::QEvdevKeyboardPlugin)
set_property(TARGET Qt5::QEvdevKeyboardPlugin PROPERTY QT_PLUGIN_TYPE "generic")
set_property(TARGET Qt5::QEvdevKeyboardPlugin PROPERTY QT_PLUGIN_EXTENDS "-")
set_property(TARGET Qt5::QEvdevKeyboardPlugin PROPERTY QT_PLUGIN_CLASS_NAME "QEvdevKeyboardPlugin")
