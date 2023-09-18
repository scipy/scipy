
add_library(Qt5::QXcbIntegrationPlugin MODULE IMPORTED)


_populate_Gui_plugin_properties(QXcbIntegrationPlugin RELEASE "platforms/libqxcb.so" FALSE)

list(APPEND Qt5Gui_PLUGINS Qt5::QXcbIntegrationPlugin)
set_property(TARGET Qt5::Gui APPEND PROPERTY QT_ALL_PLUGINS_platforms Qt5::QXcbIntegrationPlugin)
set_property(TARGET Qt5::QXcbIntegrationPlugin PROPERTY QT_PLUGIN_TYPE "platforms")
set_property(TARGET Qt5::QXcbIntegrationPlugin PROPERTY QT_PLUGIN_EXTENDS "-")
set_property(TARGET Qt5::QXcbIntegrationPlugin PROPERTY QT_PLUGIN_CLASS_NAME "QXcbIntegrationPlugin")
