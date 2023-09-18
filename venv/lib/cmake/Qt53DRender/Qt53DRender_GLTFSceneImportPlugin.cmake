
add_library(Qt5::GLTFSceneImportPlugin MODULE IMPORTED)


_populate_3DRender_plugin_properties(GLTFSceneImportPlugin RELEASE "sceneparsers/libgltfsceneimport.so" FALSE)

list(APPEND Qt53DRender_PLUGINS Qt5::GLTFSceneImportPlugin)
set_property(TARGET Qt5::3DRender APPEND PROPERTY QT_ALL_PLUGINS_sceneparsers Qt5::GLTFSceneImportPlugin)
set_property(TARGET Qt5::GLTFSceneImportPlugin PROPERTY QT_PLUGIN_TYPE "sceneparsers")
set_property(TARGET Qt5::GLTFSceneImportPlugin PROPERTY QT_PLUGIN_EXTENDS "")
set_property(TARGET Qt5::GLTFSceneImportPlugin PROPERTY QT_PLUGIN_CLASS_NAME "GLTFSceneImportPlugin")
