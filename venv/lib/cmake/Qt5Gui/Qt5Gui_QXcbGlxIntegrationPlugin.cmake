
add_library(Qt5::QXcbGlxIntegrationPlugin MODULE IMPORTED)


_populate_Gui_plugin_properties(QXcbGlxIntegrationPlugin RELEASE "xcbglintegrations/libqxcb-glx-integration.so" FALSE)

list(APPEND Qt5Gui_PLUGINS Qt5::QXcbGlxIntegrationPlugin)
set_property(TARGET Qt5::Gui APPEND PROPERTY QT_ALL_PLUGINS_xcbglintegrations Qt5::QXcbGlxIntegrationPlugin)
set_property(TARGET Qt5::QXcbGlxIntegrationPlugin PROPERTY QT_PLUGIN_TYPE "xcbglintegrations")
set_property(TARGET Qt5::QXcbGlxIntegrationPlugin PROPERTY QT_PLUGIN_EXTENDS "")
set_property(TARGET Qt5::QXcbGlxIntegrationPlugin PROPERTY QT_PLUGIN_CLASS_NAME "QXcbGlxIntegrationPlugin")
