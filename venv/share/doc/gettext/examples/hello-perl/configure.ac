dnl Example for use of GNU gettext.
dnl This file is in the public domain.
dnl
dnl Configuration file - processed by autoconf.

AC_INIT([hello-perl], [0])
AC_CONFIG_SRCDIR([hello-1.pl.in])
AM_INIT_AUTOMAKE([1.11])

dnl Check for availability of the Perl interpreter.
AC_PATH_PROG([PERL], [perl])
if test -z "$PERL"; then
  echo "*** Essential program perl not found" 1>&2
  exit 1
fi
AC_SUBST([PERL])

dnl The installed program must know where to find its message catalogs.
dnl Unfortunately, prefix gets only finally determined at the end of configure.
if test "X$prefix" = "XNONE"; then
  final_prefix="$ac_default_prefix"
else
  final_prefix="$prefix"
fi
save_prefix="$prefix"
prefix="$final_prefix"
eval "datarootdir=\"${datarootdir}\""
eval "localedir=\"${datadir}/locale\""
prefix="$save_prefix"
AC_SUBST([localedir])

dnl Support for the po directory.
AM_PO_SUBDIRS

AC_CONFIG_FILES([Makefile])
AC_CONFIG_FILES([hello-1:hello-1.pl.in], [chmod a+x hello-1])
AC_CONFIG_FILES([hello-2:hello-2.pl.in], [chmod a+x hello-2])
AC_CONFIG_FILES([m4/Makefile])
AC_CONFIG_FILES([po/Makefile], [AM_POSTPROCESS_PO_MAKEFILE])
AC_OUTPUT
