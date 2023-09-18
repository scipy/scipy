
add_library(Qt5::QPSQLDriverPlugin MODULE IMPORTED)


_populate_Sql_plugin_properties(QPSQLDriverPlugin RELEASE "sqldrivers/libqsqlpsql.so" FALSE)

list(APPEND Qt5Sql_PLUGINS Qt5::QPSQLDriverPlugin)
set_property(TARGET Qt5::Sql APPEND PROPERTY QT_ALL_PLUGINS_sqldrivers Qt5::QPSQLDriverPlugin)
set_property(TARGET Qt5::QPSQLDriverPlugin PROPERTY QT_PLUGIN_TYPE "sqldrivers")
set_property(TARGET Qt5::QPSQLDriverPlugin PROPERTY QT_PLUGIN_EXTENDS "")
set_property(TARGET Qt5::QPSQLDriverPlugin PROPERTY QT_PLUGIN_CLASS_NAME "QPSQLDriverPlugin")
