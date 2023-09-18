dnl Example for use of GNU gettext.
dnl This file is in the public domain.
dnl
dnl Configuration file - processed by autoconf.

AC_INIT([hello-csharp-forms], [0])
AC_CONFIG_SRCDIR([hello.cs])
AM_INIT_AUTOMAKE([1.11])

dnl Check whether we can build C# programs at all.
gt_CSHARPCOMP
if test -n "$HAVE_CSHARPCOMP"; then
  BUILDCSHARP=yes
else
  BUILDCSHARP=no
fi
AC_SUBST([BUILDCSHARP])

dnl Check whether we can execute C# programs.
gt_CSHARPEXEC([csharpexec-test.exe], [$srcdir/m4])
if test -n "$HAVE_CSHARPEXEC" && test $BUILDCSHARP = yes; then
  TESTCSHARP=yes
else
  TESTCSHARP=no
fi
AC_SUBST([TESTCSHARP])

dnl Checks for needed libraries.
AM_PATH_PROG_WITH_TEST([GETTEXT_WITH_GNU_GETTEXT_DLL], [gettext],
  [{ basedir=`echo "$ac_dir" | sed -e 's,/bin$,,'`; test -r "$basedir"/lib/GNU.Gettext.dll; }])
if test -z "$GETTEXT_WITH_GNU_GETTEXT_DLL"; then
  echo "Required library GNU.Gettext.dll not found." 1>&2
  exit 1
fi
changequote(,)dnl
basedir=`echo "$GETTEXT_WITH_GNU_GETTEXT_DLL" | sed -e 's,/[^/]*$,,' | sed -e 's,/bin$,,'`
changequote([, ])dnl
GNU_GETTEXT_DLL="$basedir"/lib/GNU.Gettext.dll
GNU_GETTEXT_LDADD="-L $basedir/lib"
GNU_GETTEXT_LIBS="-l GNU.Gettext"
AC_SUBST([GNU_GETTEXT_LDADD])
AC_SUBST([GNU_GETTEXT_LIBS])

dnl Support for the po directory.
AM_PO_SUBDIRS

AC_CONFIG_FILES([Makefile])
AC_CONFIG_FILES([csharpcomp.sh])
AC_CONFIG_FILES([csharpexec.sh])
AC_CONFIG_FILES([m4/Makefile])
AC_CONFIG_FILES([po/Makefile], [AM_POSTPROCESS_PO_MAKEFILE])
AC_OUTPUT
