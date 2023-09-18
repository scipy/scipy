
add_library(Qt5::QXdgDesktopPortalThemePlugin MODULE IMPORTED)


_populate_Gui_plugin_properties(QXdgDesktopPortalThemePlugin RELEASE "platformthemes/libqxdgdesktopportal.so" FALSE)

list(APPEND Qt5Gui_PLUGINS Qt5::QXdgDesktopPortalThemePlugin)
set_property(TARGET Qt5::Gui APPEND PROPERTY QT_ALL_PLUGINS_platformthemes Qt5::QXdgDesktopPortalThemePlugin)
set_property(TARGET Qt5::QXdgDesktopPortalThemePlugin PROPERTY QT_PLUGIN_TYPE "platformthemes")
set_property(TARGET Qt5::QXdgDesktopPortalThemePlugin PROPERTY QT_PLUGIN_EXTENDS "-")
set_property(TARGET Qt5::QXdgDesktopPortalThemePlugin PROPERTY QT_PLUGIN_CLASS_NAME "QXdgDesktopPortalThemePlugin")
