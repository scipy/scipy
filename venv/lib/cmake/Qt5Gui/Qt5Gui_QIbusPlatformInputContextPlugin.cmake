
add_library(Qt5::QIbusPlatformInputContextPlugin MODULE IMPORTED)


_populate_Gui_plugin_properties(QIbusPlatformInputContextPlugin RELEASE "platforminputcontexts/libibusplatforminputcontextplugin.so" FALSE)

list(APPEND Qt5Gui_PLUGINS Qt5::QIbusPlatformInputContextPlugin)
set_property(TARGET Qt5::Gui APPEND PROPERTY QT_ALL_PLUGINS_platforminputcontexts Qt5::QIbusPlatformInputContextPlugin)
set_property(TARGET Qt5::QIbusPlatformInputContextPlugin PROPERTY QT_PLUGIN_TYPE "platforminputcontexts")
set_property(TARGET Qt5::QIbusPlatformInputContextPlugin PROPERTY QT_PLUGIN_EXTENDS "-")
set_property(TARGET Qt5::QIbusPlatformInputContextPlugin PROPERTY QT_PLUGIN_CLASS_NAME "QIbusPlatformInputContextPlugin")
