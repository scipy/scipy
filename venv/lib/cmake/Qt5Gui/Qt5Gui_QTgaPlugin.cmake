
add_library(Qt5::QTgaPlugin MODULE IMPORTED)


_populate_Gui_plugin_properties(QTgaPlugin RELEASE "imageformats/libqtga.so" FALSE)

list(APPEND Qt5Gui_PLUGINS Qt5::QTgaPlugin)
set_property(TARGET Qt5::Gui APPEND PROPERTY QT_ALL_PLUGINS_imageformats Qt5::QTgaPlugin)
set_property(TARGET Qt5::QTgaPlugin PROPERTY QT_PLUGIN_TYPE "imageformats")
set_property(TARGET Qt5::QTgaPlugin PROPERTY QT_PLUGIN_EXTENDS "")
set_property(TARGET Qt5::QTgaPlugin PROPERTY QT_PLUGIN_CLASS_NAME "QTgaPlugin")
