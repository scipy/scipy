
add_library(Qt5::QAlsaPlugin MODULE IMPORTED)


_populate_Multimedia_plugin_properties(QAlsaPlugin RELEASE "audio/libqtaudio_alsa.so" FALSE)

list(APPEND Qt5Multimedia_PLUGINS Qt5::QAlsaPlugin)
set_property(TARGET Qt5::Multimedia APPEND PROPERTY QT_ALL_PLUGINS_audio Qt5::QAlsaPlugin)
set_property(TARGET Qt5::QAlsaPlugin PROPERTY QT_PLUGIN_TYPE "audio")
set_property(TARGET Qt5::QAlsaPlugin PROPERTY QT_PLUGIN_EXTENDS "")
set_property(TARGET Qt5::QAlsaPlugin PROPERTY QT_PLUGIN_CLASS_NAME "QAlsaPlugin")
