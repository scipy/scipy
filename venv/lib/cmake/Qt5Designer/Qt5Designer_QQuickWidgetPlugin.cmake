
add_library(Qt5::QQuickWidgetPlugin MODULE IMPORTED)


_populate_Designer_plugin_properties(QQuickWidgetPlugin RELEASE "designer/libqquickwidget.so" FALSE)

list(APPEND Qt5Designer_PLUGINS Qt5::QQuickWidgetPlugin)
set_property(TARGET Qt5::Designer APPEND PROPERTY QT_ALL_PLUGINS_designer Qt5::QQuickWidgetPlugin)
set_property(TARGET Qt5::QQuickWidgetPlugin PROPERTY QT_PLUGIN_TYPE "designer")
set_property(TARGET Qt5::QQuickWidgetPlugin PROPERTY QT_PLUGIN_EXTENDS "")
set_property(TARGET Qt5::QQuickWidgetPlugin PROPERTY QT_PLUGIN_CLASS_NAME "QQuickWidgetPlugin")
