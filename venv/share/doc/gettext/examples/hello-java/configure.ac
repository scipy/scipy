dnl Example for use of GNU gettext.
dnl This file is in the public domain.
dnl
dnl Configuration file - processed by autoconf.

AC_INIT([hello-java], [0])
AC_CONFIG_SRCDIR([Hello.java])
AM_INIT_AUTOMAKE([1.11])

dnl Check whether we can execute Java programs.
gt_JAVAEXEC
dnl Check whether we can build Java programs.
gt_JAVACOMP([1.5])
AC_CHECK_PROG([JAR], [jar], [jar])
if test -n "$HAVE_JAVACOMP" && test -n "$JAR"; then
  BUILDJAVA=yes
else
  BUILDJAVA=no
fi
AC_SUBST([BUILDJAVA])
if test -n "$HAVE_JAVAEXEC" && test "$BUILDJAVA" = yes; then
  TESTJAVA=yes
else
  TESTJAVA=no
fi
AC_SUBST([TESTJAVA])

dnl Checks for compiler output filename suffixes.
AC_OBJEXT
AC_EXEEXT

dnl Checks for needed libraries.
AM_PATH_PROG_WITH_TEST([GETTEXT_WITH_LIBINTL_JAR], [gettext],
  [{ basedir=`echo "$ac_dir" | sed -e 's,/bin$,,'`; test -r "$basedir"/share/gettext/libintl.jar; }])
if test -z "$GETTEXT_WITH_LIBINTL_JAR"; then
  echo "Required library libintl.jar not found." 1>&2
  exit 1
fi
changequote(,)dnl
basedir=`echo "$GETTEXT_WITH_LIBINTL_JAR" | sed -e 's,/[^/]*$,,' | sed -e 's,/bin$,,'`
changequote([, ])dnl
LIBINTL_JAR="$basedir"/share/gettext/libintl.jar
AC_SUBST([LIBINTL_JAR])

dnl Support for the po directory.
AM_PO_SUBDIRS

AC_CONFIG_FILES([Makefile])
AC_CONFIG_FILES([javacomp.sh])
AC_CONFIG_FILES([javaexec.sh])
AC_CONFIG_FILES([m4/Makefile])
AC_CONFIG_FILES([po/Makefile], [AM_POSTPROCESS_PO_MAKEFILE])
AC_OUTPUT
