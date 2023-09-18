
add_library(Qt5::QTiffPlugin MODULE IMPORTED)


_populate_Gui_plugin_properties(QTiffPlugin RELEASE "imageformats/libqtiff.so" FALSE)

list(APPEND Qt5Gui_PLUGINS Qt5::QTiffPlugin)
set_property(TARGET Qt5::Gui APPEND PROPERTY QT_ALL_PLUGINS_imageformats Qt5::QTiffPlugin)
set_property(TARGET Qt5::QTiffPlugin PROPERTY QT_PLUGIN_TYPE "imageformats")
set_property(TARGET Qt5::QTiffPlugin PROPERTY QT_PLUGIN_EXTENDS "")
set_property(TARGET Qt5::QTiffPlugin PROPERTY QT_PLUGIN_CLASS_NAME "QTiffPlugin")
