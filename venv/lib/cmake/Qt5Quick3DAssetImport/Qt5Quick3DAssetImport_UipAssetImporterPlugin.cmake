
add_library(Qt5::UipAssetImporterPlugin MODULE IMPORTED)


_populate_Quick3DAssetImport_plugin_properties(UipAssetImporterPlugin RELEASE "assetimporters/libuip.so" FALSE)

list(APPEND Qt5Quick3DAssetImport_PLUGINS Qt5::UipAssetImporterPlugin)
set_property(TARGET Qt5::Quick3DAssetImport APPEND PROPERTY QT_ALL_PLUGINS_assetimporters Qt5::UipAssetImporterPlugin)
set_property(TARGET Qt5::UipAssetImporterPlugin PROPERTY QT_PLUGIN_TYPE "assetimporters")
set_property(TARGET Qt5::UipAssetImporterPlugin PROPERTY QT_PLUGIN_EXTENDS "")
set_property(TARGET Qt5::UipAssetImporterPlugin PROPERTY QT_PLUGIN_CLASS_NAME "UipAssetImporterPlugin")
