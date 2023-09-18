
add_library(Qt5::QSvgIconPlugin MODULE IMPORTED)


_populate_Gui_plugin_properties(QSvgIconPlugin RELEASE "iconengines/libqsvgicon.so" FALSE)

list(APPEND Qt5Gui_PLUGINS Qt5::QSvgIconPlugin)
set_property(TARGET Qt5::Gui APPEND PROPERTY QT_ALL_PLUGINS_iconengines Qt5::QSvgIconPlugin)
set_property(TARGET Qt5::QSvgIconPlugin PROPERTY QT_PLUGIN_TYPE "iconengines")
set_property(TARGET Qt5::QSvgIconPlugin PROPERTY QT_PLUGIN_EXTENDS "Qt::Svg")
set_property(TARGET Qt5::QSvgIconPlugin PROPERTY QT_PLUGIN_CLASS_NAME "QSvgIconPlugin")
