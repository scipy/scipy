
add_library(Qt5::QOffscreenIntegrationPlugin MODULE IMPORTED)


_populate_Gui_plugin_properties(QOffscreenIntegrationPlugin RELEASE "platforms/libqoffscreen.so" FALSE)

list(APPEND Qt5Gui_PLUGINS Qt5::QOffscreenIntegrationPlugin)
set_property(TARGET Qt5::Gui APPEND PROPERTY QT_ALL_PLUGINS_platforms Qt5::QOffscreenIntegrationPlugin)
set_property(TARGET Qt5::QOffscreenIntegrationPlugin PROPERTY QT_PLUGIN_TYPE "platforms")
set_property(TARGET Qt5::QOffscreenIntegrationPlugin PROPERTY QT_PLUGIN_EXTENDS "-")
set_property(TARGET Qt5::QOffscreenIntegrationPlugin PROPERTY QT_PLUGIN_CLASS_NAME "QOffscreenIntegrationPlugin")
