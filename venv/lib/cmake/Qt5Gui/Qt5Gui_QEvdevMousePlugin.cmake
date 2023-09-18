
add_library(Qt5::QEvdevMousePlugin MODULE IMPORTED)


_populate_Gui_plugin_properties(QEvdevMousePlugin RELEASE "generic/libqevdevmouseplugin.so" FALSE)

list(APPEND Qt5Gui_PLUGINS Qt5::QEvdevMousePlugin)
set_property(TARGET Qt5::Gui APPEND PROPERTY QT_ALL_PLUGINS_generic Qt5::QEvdevMousePlugin)
set_property(TARGET Qt5::QEvdevMousePlugin PROPERTY QT_PLUGIN_TYPE "generic")
set_property(TARGET Qt5::QEvdevMousePlugin PROPERTY QT_PLUGIN_EXTENDS "-")
set_property(TARGET Qt5::QEvdevMousePlugin PROPERTY QT_PLUGIN_CLASS_NAME "QEvdevMousePlugin")
