
add_library(Qt5::QXcbEglIntegrationPlugin MODULE IMPORTED)


_populate_Gui_plugin_properties(QXcbEglIntegrationPlugin RELEASE "xcbglintegrations/libqxcb-egl-integration.so" FALSE)

list(APPEND Qt5Gui_PLUGINS Qt5::QXcbEglIntegrationPlugin)
set_property(TARGET Qt5::Gui APPEND PROPERTY QT_ALL_PLUGINS_xcbglintegrations Qt5::QXcbEglIntegrationPlugin)
set_property(TARGET Qt5::QXcbEglIntegrationPlugin PROPERTY QT_PLUGIN_TYPE "xcbglintegrations")
set_property(TARGET Qt5::QXcbEglIntegrationPlugin PROPERTY QT_PLUGIN_EXTENDS "")
set_property(TARGET Qt5::QXcbEglIntegrationPlugin PROPERTY QT_PLUGIN_CLASS_NAME "QXcbEglIntegrationPlugin")
