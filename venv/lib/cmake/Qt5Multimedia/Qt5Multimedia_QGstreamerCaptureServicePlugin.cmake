
add_library(Qt5::QGstreamerCaptureServicePlugin MODULE IMPORTED)


_populate_Multimedia_plugin_properties(QGstreamerCaptureServicePlugin RELEASE "mediaservice/libgstmediacapture.so" FALSE)

list(APPEND Qt5Multimedia_PLUGINS Qt5::QGstreamerCaptureServicePlugin)
set_property(TARGET Qt5::Multimedia APPEND PROPERTY QT_ALL_PLUGINS_mediaservice Qt5::QGstreamerCaptureServicePlugin)
set_property(TARGET Qt5::QGstreamerCaptureServicePlugin PROPERTY QT_PLUGIN_TYPE "mediaservice")
set_property(TARGET Qt5::QGstreamerCaptureServicePlugin PROPERTY QT_PLUGIN_EXTENDS "")
set_property(TARGET Qt5::QGstreamerCaptureServicePlugin PROPERTY QT_PLUGIN_CLASS_NAME "QGstreamerCaptureServicePlugin")
