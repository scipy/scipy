
add_library(Qt5::QGstreamerPlayerServicePlugin MODULE IMPORTED)


_populate_Multimedia_plugin_properties(QGstreamerPlayerServicePlugin RELEASE "mediaservice/libgstmediaplayer.so" FALSE)

list(APPEND Qt5Multimedia_PLUGINS Qt5::QGstreamerPlayerServicePlugin)
set_property(TARGET Qt5::Multimedia APPEND PROPERTY QT_ALL_PLUGINS_mediaservice Qt5::QGstreamerPlayerServicePlugin)
set_property(TARGET Qt5::QGstreamerPlayerServicePlugin PROPERTY QT_PLUGIN_TYPE "mediaservice")
set_property(TARGET Qt5::QGstreamerPlayerServicePlugin PROPERTY QT_PLUGIN_EXTENDS "")
set_property(TARGET Qt5::QGstreamerPlayerServicePlugin PROPERTY QT_PLUGIN_CLASS_NAME "QGstreamerPlayerServicePlugin")
