
add_library(Qt5::QPulseAudioPlugin MODULE IMPORTED)


_populate_Multimedia_plugin_properties(QPulseAudioPlugin RELEASE "audio/libqtmedia_pulse.so" FALSE)

list(APPEND Qt5Multimedia_PLUGINS Qt5::QPulseAudioPlugin)
set_property(TARGET Qt5::Multimedia APPEND PROPERTY QT_ALL_PLUGINS_audio Qt5::QPulseAudioPlugin)
set_property(TARGET Qt5::QPulseAudioPlugin PROPERTY QT_PLUGIN_TYPE "audio")
set_property(TARGET Qt5::QPulseAudioPlugin PROPERTY QT_PLUGIN_EXTENDS "")
set_property(TARGET Qt5::QPulseAudioPlugin PROPERTY QT_PLUGIN_CLASS_NAME "QPulseAudioPlugin")
