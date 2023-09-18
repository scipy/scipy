
add_library(Qt5::QSvgPlugin MODULE IMPORTED)


_populate_Gui_plugin_properties(QSvgPlugin RELEASE "imageformats/libqsvg.so" FALSE)

list(APPEND Qt5Gui_PLUGINS Qt5::QSvgPlugin)
set_property(TARGET Qt5::Gui APPEND PROPERTY QT_ALL_PLUGINS_imageformats Qt5::QSvgPlugin)
set_property(TARGET Qt5::QSvgPlugin PROPERTY QT_PLUGIN_TYPE "imageformats")
set_property(TARGET Qt5::QSvgPlugin PROPERTY QT_PLUGIN_EXTENDS "Qt::Svg")
set_property(TARGET Qt5::QSvgPlugin PROPERTY QT_PLUGIN_CLASS_NAME "QSvgPlugin")
