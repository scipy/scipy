dnl Example for use of GNU gettext.
dnl This file is in the public domain.
dnl
dnl Configuration file - processed by autoconf.

AC_INIT([hello-ycp], [0])
AC_CONFIG_SRCDIR([hello.ycp])
AM_INIT_AUTOMAKE([1.11])

dnl Check for availability of YaST's engine.
test -f /usr/lib/YaST2/bin/y2base || {
  echo "*** Essential program y2base not found" 1>&2
  exit 1
}

dnl The installed program must know where to find its message catalogs.
dnl Unfortunately, YaST hard codes the message catalog directory.
localedir=/usr/share/YaST2/locale
AC_SUBST([localedir])

dnl Support for the po directory.
AM_PO_SUBDIRS

AC_CONFIG_FILES([Makefile])
AC_CONFIG_FILES([m4/Makefile])
AC_CONFIG_FILES([po/Makefile], [AM_POSTPROCESS_PO_MAKEFILE])
AC_OUTPUT
