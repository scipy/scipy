dnl Example for use of GNU gettext.
dnl This file is in the public domain.
dnl
dnl Configuration file - processed by autoconf.

AC_INIT([hello-c++-gnome], [0], , [hello-c++-gnome])
AC_CONFIG_SRCDIR([hello.cc])
AM_INIT_AUTOMAKE([1.11])

AC_PROG_CXX

GNOME_INIT
GTKMM_CFLAGS=`gtkmm-config --cflags`
AC_SUBST([GTKMM_CFLAGS])
GTKMM_LIBS=`gtkmm-config --libs`
AC_SUBST([GTKMM_LIBS])

AC_CHECK_HEADERS([unistd.h])
AM_GNU_GETTEXT([external])
AM_GNU_GETTEXT_VERSION([0.21.1])

AC_CONFIG_FILES([Makefile])
AC_CONFIG_FILES([m4/Makefile])
AC_CONFIG_FILES([po/Makefile.in])
AC_OUTPUT
