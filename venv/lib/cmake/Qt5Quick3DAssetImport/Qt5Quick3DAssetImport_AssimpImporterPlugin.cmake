
add_library(Qt5::AssimpImporterPlugin MODULE IMPORTED)


_populate_Quick3DAssetImport_plugin_properties(AssimpImporterPlugin RELEASE "assetimporters/libassimp.so" FALSE)

list(APPEND Qt5Quick3DAssetImport_PLUGINS Qt5::AssimpImporterPlugin)
set_property(TARGET Qt5::Quick3DAssetImport APPEND PROPERTY QT_ALL_PLUGINS_assetimporters Qt5::AssimpImporterPlugin)
set_property(TARGET Qt5::AssimpImporterPlugin PROPERTY QT_PLUGIN_TYPE "assetimporters")
set_property(TARGET Qt5::AssimpImporterPlugin PROPERTY QT_PLUGIN_EXTENDS "")
set_property(TARGET Qt5::AssimpImporterPlugin PROPERTY QT_PLUGIN_CLASS_NAME "AssimpImporterPlugin")
