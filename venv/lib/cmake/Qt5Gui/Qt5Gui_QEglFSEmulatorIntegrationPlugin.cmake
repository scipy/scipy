
add_library(Qt5::QEglFSEmulatorIntegrationPlugin MODULE IMPORTED)


_populate_Gui_plugin_properties(QEglFSEmulatorIntegrationPlugin RELEASE "egldeviceintegrations/libqeglfs-emu-integration.so" FALSE)

list(APPEND Qt5Gui_PLUGINS Qt5::QEglFSEmulatorIntegrationPlugin)
set_property(TARGET Qt5::Gui APPEND PROPERTY QT_ALL_PLUGINS_egldeviceintegrations Qt5::QEglFSEmulatorIntegrationPlugin)
set_property(TARGET Qt5::QEglFSEmulatorIntegrationPlugin PROPERTY QT_PLUGIN_TYPE "egldeviceintegrations")
set_property(TARGET Qt5::QEglFSEmulatorIntegrationPlugin PROPERTY QT_PLUGIN_EXTENDS "")
set_property(TARGET Qt5::QEglFSEmulatorIntegrationPlugin PROPERTY QT_PLUGIN_CLASS_NAME "QEglFSEmulatorIntegrationPlugin")
