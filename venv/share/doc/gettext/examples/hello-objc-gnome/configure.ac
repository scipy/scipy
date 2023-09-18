dnl Example for use of GNU gettext.
dnl This file is in the public domain.
dnl
dnl Configuration file - processed by autoconf.

AC_INIT([hello-objc-gnome], [0])
AC_CONFIG_SRCDIR([hello.m])
AM_INIT_AUTOMAKE([1.11])

AC_PROG_CC
AC_SUBST([OBJC],["$CC"])
AC_SUBST([OBJCFLAGS],["$CFLAGS"])

GNOME_INIT
OBGNOME_INCLUDES=
obgnome_prefix=`obgnome-config --prefix`
if test -n "$obgnome_prefix"; then
  OBGNOME_INCLUDES="-I $obgnome_prefix/include"
fi
AC_SUBST([OBGNOME_INCLUDES])

AC_CHECK_HEADERS([unistd.h])
AM_GNU_GETTEXT([external])
AM_GNU_GETTEXT_VERSION([0.21.1])

AC_CONFIG_FILES([Makefile])
AC_CONFIG_FILES([m4/Makefile])
AC_CONFIG_FILES([po/Makefile.in])
AC_OUTPUT
