dnl Example for use of GNU gettext.
dnl This file is in the public domain.
dnl
dnl Configuration file - processed by autoconf.

AC_INIT([hello-tcl-tk], [0])
AC_CONFIG_SRCDIR([hello.tcl])
AM_INIT_AUTOMAKE([1.11])

dnl Check for availability of the Tcl/Tk interpreter.
AC_PATH_PROG([WISH], [wish])
if test -z "$WISH"; then
  echo "*** Essential program wish not found" 1>&2
  exit 1
fi
AC_SUBST([WISH])

dnl The installed program must know where to find its message catalogs.
dnl Unfortunately, prefix gets only finally determined at the end of configure.
if test "X$prefix" = "XNONE"; then
  final_prefix="$ac_default_prefix"
else
  final_prefix="$prefix"
fi
save_prefix="$prefix"
prefix="$final_prefix"
save_datarootdir="$datarootdir"
eval "datarootdir=\"${datarootdir}\""
save_datadir="$datadir"
eval "datadir=\"${datadir}\""
pkgdatadir="${datadir}/${PACKAGE}"
datadir="$save_datadir"
datarootdir="$save_datarootdir"
prefix="$save_prefix"
AC_SUBST([pkgdatadir])

dnl Support for the po directory.
AM_PO_SUBDIRS

AC_CONFIG_FILES([Makefile])
AC_CONFIG_FILES([hello:hello.tcl], [chmod a+x hello])
AC_CONFIG_FILES([m4/Makefile])
AC_CONFIG_FILES([po/Makefile], [AM_POSTPROCESS_PO_MAKEFILE])
AC_OUTPUT
