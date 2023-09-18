
add_library(Qt5::QEglFSIntegrationPlugin MODULE IMPORTED)


_populate_Gui_plugin_properties(QEglFSIntegrationPlugin RELEASE "platforms/libqeglfs.so" FALSE)

list(APPEND Qt5Gui_PLUGINS Qt5::QEglFSIntegrationPlugin)
set_property(TARGET Qt5::Gui APPEND PROPERTY QT_ALL_PLUGINS_platforms Qt5::QEglFSIntegrationPlugin)
set_property(TARGET Qt5::QEglFSIntegrationPlugin PROPERTY QT_PLUGIN_TYPE "platforms")
set_property(TARGET Qt5::QEglFSIntegrationPlugin PROPERTY QT_PLUGIN_EXTENDS "-")
set_property(TARGET Qt5::QEglFSIntegrationPlugin PROPERTY QT_PLUGIN_CLASS_NAME "QEglFSIntegrationPlugin")
