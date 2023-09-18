
add_library(Qt5::QEglFSKmsEglDeviceIntegrationPlugin MODULE IMPORTED)


_populate_Gui_plugin_properties(QEglFSKmsEglDeviceIntegrationPlugin RELEASE "egldeviceintegrations/libqeglfs-kms-egldevice-integration.so" FALSE)

list(APPEND Qt5Gui_PLUGINS Qt5::QEglFSKmsEglDeviceIntegrationPlugin)
set_property(TARGET Qt5::Gui APPEND PROPERTY QT_ALL_PLUGINS_egldeviceintegrations Qt5::QEglFSKmsEglDeviceIntegrationPlugin)
set_property(TARGET Qt5::QEglFSKmsEglDeviceIntegrationPlugin PROPERTY QT_PLUGIN_TYPE "egldeviceintegrations")
set_property(TARGET Qt5::QEglFSKmsEglDeviceIntegrationPlugin PROPERTY QT_PLUGIN_EXTENDS "")
set_property(TARGET Qt5::QEglFSKmsEglDeviceIntegrationPlugin PROPERTY QT_PLUGIN_CLASS_NAME "QEglFSKmsEglDeviceIntegrationPlugin")
