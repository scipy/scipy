
add_library(Qt5::QCupsPrinterSupportPlugin MODULE IMPORTED)


_populate_PrintSupport_plugin_properties(QCupsPrinterSupportPlugin RELEASE "printsupport/libcupsprintersupport.so" FALSE)

list(APPEND Qt5PrintSupport_PLUGINS Qt5::QCupsPrinterSupportPlugin)
set_property(TARGET Qt5::PrintSupport APPEND PROPERTY QT_ALL_PLUGINS_printsupport Qt5::QCupsPrinterSupportPlugin)
set_property(TARGET Qt5::QCupsPrinterSupportPlugin PROPERTY QT_PLUGIN_TYPE "printsupport")
set_property(TARGET Qt5::QCupsPrinterSupportPlugin PROPERTY QT_PLUGIN_EXTENDS "")
set_property(TARGET Qt5::QCupsPrinterSupportPlugin PROPERTY QT_PLUGIN_CLASS_NAME "QCupsPrinterSupportPlugin")
