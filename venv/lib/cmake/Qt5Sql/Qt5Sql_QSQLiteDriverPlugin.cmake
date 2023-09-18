
add_library(Qt5::QSQLiteDriverPlugin MODULE IMPORTED)


_populate_Sql_plugin_properties(QSQLiteDriverPlugin RELEASE "sqldrivers/libqsqlite.so" FALSE)

list(APPEND Qt5Sql_PLUGINS Qt5::QSQLiteDriverPlugin)
set_property(TARGET Qt5::Sql APPEND PROPERTY QT_ALL_PLUGINS_sqldrivers Qt5::QSQLiteDriverPlugin)
set_property(TARGET Qt5::QSQLiteDriverPlugin PROPERTY QT_PLUGIN_TYPE "sqldrivers")
set_property(TARGET Qt5::QSQLiteDriverPlugin PROPERTY QT_PLUGIN_EXTENDS "")
set_property(TARGET Qt5::QSQLiteDriverPlugin PROPERTY QT_PLUGIN_CLASS_NAME "QSQLiteDriverPlugin")
