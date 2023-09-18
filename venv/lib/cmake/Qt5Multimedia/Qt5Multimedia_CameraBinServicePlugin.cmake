
add_library(Qt5::CameraBinServicePlugin MODULE IMPORTED)


_populate_Multimedia_plugin_properties(CameraBinServicePlugin RELEASE "mediaservice/libgstcamerabin.so" FALSE)

list(APPEND Qt5Multimedia_PLUGINS Qt5::CameraBinServicePlugin)
set_property(TARGET Qt5::Multimedia APPEND PROPERTY QT_ALL_PLUGINS_mediaservice Qt5::CameraBinServicePlugin)
set_property(TARGET Qt5::CameraBinServicePlugin PROPERTY QT_PLUGIN_TYPE "mediaservice")
set_property(TARGET Qt5::CameraBinServicePlugin PROPERTY QT_PLUGIN_EXTENDS "")
set_property(TARGET Qt5::CameraBinServicePlugin PROPERTY QT_PLUGIN_CLASS_NAME "CameraBinServicePlugin")
