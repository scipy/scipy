
add_library(Qt5::QWebpPlugin MODULE IMPORTED)


_populate_Gui_plugin_properties(QWebpPlugin RELEASE "imageformats/libqwebp.so" FALSE)

list(APPEND Qt5Gui_PLUGINS Qt5::QWebpPlugin)
set_property(TARGET Qt5::Gui APPEND PROPERTY QT_ALL_PLUGINS_imageformats Qt5::QWebpPlugin)
set_property(TARGET Qt5::QWebpPlugin PROPERTY QT_PLUGIN_TYPE "imageformats")
set_property(TARGET Qt5::QWebpPlugin PROPERTY QT_PLUGIN_EXTENDS "")
set_property(TARGET Qt5::QWebpPlugin PROPERTY QT_PLUGIN_CLASS_NAME "QWebpPlugin")
