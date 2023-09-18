
add_library(Qt5::AssimpSceneImportPlugin MODULE IMPORTED)


_populate_3DRender_plugin_properties(AssimpSceneImportPlugin RELEASE "sceneparsers/libassimpsceneimport.so" FALSE)

list(APPEND Qt53DRender_PLUGINS Qt5::AssimpSceneImportPlugin)
set_property(TARGET Qt5::3DRender APPEND PROPERTY QT_ALL_PLUGINS_sceneparsers Qt5::AssimpSceneImportPlugin)
set_property(TARGET Qt5::AssimpSceneImportPlugin PROPERTY QT_PLUGIN_TYPE "sceneparsers")
set_property(TARGET Qt5::AssimpSceneImportPlugin PROPERTY QT_PLUGIN_EXTENDS "")
set_property(TARGET Qt5::AssimpSceneImportPlugin PROPERTY QT_PLUGIN_CLASS_NAME "AssimpSceneImportPlugin")
