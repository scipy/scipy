
add_library(Qt5::QMYSQLDriverPlugin MODULE IMPORTED)


_populate_Sql_plugin_properties(QMYSQLDriverPlugin RELEASE "sqldrivers/libqsqlmysql.so" FALSE)

list(APPEND Qt5Sql_PLUGINS Qt5::QMYSQLDriverPlugin)
set_property(TARGET Qt5::Sql APPEND PROPERTY QT_ALL_PLUGINS_sqldrivers Qt5::QMYSQLDriverPlugin)
set_property(TARGET Qt5::QMYSQLDriverPlugin PROPERTY QT_PLUGIN_TYPE "sqldrivers")
set_property(TARGET Qt5::QMYSQLDriverPlugin PROPERTY QT_PLUGIN_EXTENDS "")
set_property(TARGET Qt5::QMYSQLDriverPlugin PROPERTY QT_PLUGIN_CLASS_NAME "QMYSQLDriverPlugin")
