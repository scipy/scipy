dnl Example for use of GNU gettext.
dnl This file is in the public domain.
dnl
dnl Configuration file - processed by autoconf.

AC_INIT([hello-c++-wxwidgets], [0], , [hello-c++-wxwidgets])
AC_CONFIG_SRCDIR([hello.cc])
AM_INIT_AUTOMAKE([1.11])

AC_PROG_CXX
AC_CHECK_HEADERS([unistd.h])

AM_OPTIONS_WXCONFIG
AM_PATH_WXCONFIG([2.6.0], [], [AC_MSG_ERROR([wxWidgets not found])])
AC_LIB_LINKFLAGS_FROM_LIBS([WX_LDADD], [$WX_LIBS])

dnl Support for the po directory.
AM_PO_SUBDIRS

AC_CONFIG_FILES([Makefile])
AC_CONFIG_FILES([m4/Makefile])
AC_CONFIG_FILES([po/Makefile], [AM_POSTPROCESS_PO_MAKEFILE])
AC_OUTPUT
