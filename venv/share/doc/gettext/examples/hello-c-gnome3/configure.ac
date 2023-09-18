dnl Example for use of GNU gettext.
dnl This file is in the public domain.
dnl
dnl Configuration file - processed by autoconf.

AC_INIT([hello-c-gnome3], [0])
AC_CONFIG_SRCDIR([hello.c])
AM_INIT_AUTOMAKE([1.11])

AC_PROG_CC
AC_CHECK_HEADERS([unistd.h])
AM_GNU_GETTEXT([external])
AM_GNU_GETTEXT_VERSION([0.21.1])

dnl Check GNOME specific stuff.
dnl
dnl If you have full GNOME development environment installed on your
dnl system, you should be able to use the following macros:
dnl
dnl   AM_PATH_GLIB_2_0
dnl   PKG_CHECK_MODULES([GTK], [gtk+-3.0 >= 3.10], ,
dnl     [AC_MSG_ERROR([can't find gtk+-3.0])])
dnl
dnl Here, in gettext-tools/examples, we do the checks manually for
dnl simplicity.

AC_PATH_PROG([GLIB_COMPILE_RESOURCES], [glib-compile-resources])
AS_IF([test -z "$GLIB_COMPILE_RESOURCES"], [
  AC_MSG_ERROR([can't find glib-compile-resources])
])

AC_PATH_PROG([GLIB_COMPILE_SCHEMAS], [glib-compile-schemas])
AS_IF([test -z "$GLIB_COMPILE_SCHEMAS"], [
  AC_MSG_ERROR([can't find glib-compile-schemas])
])

AC_PATH_PROG([PKG_CONFIG], [pkg-config])
AS_IF([test -z "$PKG_CONFIG"], [
  AC_MSG_ERROR([can't find pkg-config])
])

GTK="gtk+-3.0 >= 3.10"
AS_IF(["$PKG_CONFIG" --exists "$GTK"], , [
  AC_MSG_ERROR([can't find $GTK])
])

GTK_CFLAGS=`"$PKG_CONFIG" --cflags "$GTK"`
AC_SUBST([GTK_CFLAGS])

GTK_LIBS=`"$PKG_CONFIG" --libs "$GTK"`
AC_SUBST([GTK_LIBS])

AC_CONFIG_FILES([Makefile])
AC_CONFIG_FILES([m4/Makefile])
AC_CONFIG_FILES([po/Makefile.in])
AC_OUTPUT
