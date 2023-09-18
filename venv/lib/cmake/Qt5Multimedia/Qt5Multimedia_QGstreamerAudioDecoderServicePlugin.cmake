
add_library(Qt5::QGstreamerAudioDecoderServicePlugin MODULE IMPORTED)


_populate_Multimedia_plugin_properties(QGstreamerAudioDecoderServicePlugin RELEASE "mediaservice/libgstaudiodecoder.so" FALSE)

list(APPEND Qt5Multimedia_PLUGINS Qt5::QGstreamerAudioDecoderServicePlugin)
set_property(TARGET Qt5::Multimedia APPEND PROPERTY QT_ALL_PLUGINS_mediaservice Qt5::QGstreamerAudioDecoderServicePlugin)
set_property(TARGET Qt5::QGstreamerAudioDecoderServicePlugin PROPERTY QT_PLUGIN_TYPE "mediaservice")
set_property(TARGET Qt5::QGstreamerAudioDecoderServicePlugin PROPERTY QT_PLUGIN_EXTENDS "")
set_property(TARGET Qt5::QGstreamerAudioDecoderServicePlugin PROPERTY QT_PLUGIN_CLASS_NAME "QGstreamerAudioDecoderServicePlugin")
