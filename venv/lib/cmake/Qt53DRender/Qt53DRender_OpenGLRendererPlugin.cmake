
add_library(Qt5::OpenGLRendererPlugin MODULE IMPORTED)


_populate_3DRender_plugin_properties(OpenGLRendererPlugin RELEASE "renderers/libopenglrenderer.so" FALSE)

list(APPEND Qt53DRender_PLUGINS Qt5::OpenGLRendererPlugin)
set_property(TARGET Qt5::3DRender APPEND PROPERTY QT_ALL_PLUGINS_renderers Qt5::OpenGLRendererPlugin)
set_property(TARGET Qt5::OpenGLRendererPlugin PROPERTY QT_PLUGIN_TYPE "renderers")
set_property(TARGET Qt5::OpenGLRendererPlugin PROPERTY QT_PLUGIN_EXTENDS "")
set_property(TARGET Qt5::OpenGLRendererPlugin PROPERTY QT_PLUGIN_CLASS_NAME "OpenGLRendererPlugin")
