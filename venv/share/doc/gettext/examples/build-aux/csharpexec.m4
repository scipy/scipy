# csharpexec.m4 serial 8
dnl Copyright (C) 2003-2005, 2009-2022 Free Software Foundation, Inc.
dnl This file is free software; the Free Software Foundation
dnl gives unlimited permission to copy and/or distribute it,
dnl with or without modifications, as long as this notice is preserved.

# Prerequisites of csharpexec.sh.
# Checks for a C# execution engine.
# gt_CSHARPEXEC or gt_CSHARPEXEC(testexecutable, its-directory)
# Sets at most one of HAVE_MONO, HAVE_CLIX.
# Sets HAVE_CSHARPEXEC to nonempty if csharpexec.sh will work.
AC_DEFUN([gt_CSHARPEXEC],
[
  AC_REQUIRE([gt_CSHARP_CHOICE])
  AC_REQUIRE([AC_CANONICAL_HOST]) dnl for cross-compiles
  AC_MSG_CHECKING([for C[#] program execution engine])
  AC_EGREP_CPP([yes], [
#if defined _WIN32 || defined __EMX__ || defined __DJGPP__
  yes
#endif
], MONO_PATH_SEPARATOR=';', MONO_PATH_SEPARATOR=':')
  HAVE_CSHARPEXEC=1
  pushdef([AC_MSG_CHECKING],[:])dnl
  pushdef([AC_CHECKING],[:])dnl
  pushdef([AC_MSG_RESULT],[:])dnl
  AC_CHECK_PROG([HAVE_MONO_IN_PATH], [mono], [yes])
  AC_CHECK_PROG([HAVE_CLIX_IN_PATH], [clix], [yes])
  popdef([AC_MSG_RESULT])dnl
  popdef([AC_CHECKING])dnl
  popdef([AC_MSG_CHECKING])dnl
  for impl in "$CSHARP_CHOICE" mono no; do
    case "$impl" in
      mono)
        if test -n "$HAVE_MONO_IN_PATH" \
           && mono --version >/dev/null 2>/dev/null \
           m4_if([$1], , , [&& mono $2/$1 >/dev/null 2>/dev/null]); then
          HAVE_MONO=1
          ac_result="mono"
          break
        fi
        ;;
      sscli)
        if test -n "$HAVE_CLIX_IN_PATH" \
           m4_if([$1], , , [&& clix $2/$1 >/dev/null 2>/dev/null]); then
          HAVE_CLIX=1
          case $host_os in
            cygwin* | mingw* | pw32*)
              CLIX_PATH_VAR=PATH
              ;;
            darwin* | rhapsody*)
              CLIX_PATH_VAR=DYLD_LIBRARY_PATH
              ;;
            *)
              CLIX_PATH_VAR=LD_LIBRARY_PATH
              ;;
          esac
          eval CLIX_PATH=\"\$CLIX_PATH_VAR\"
          ac_result="clix"
          break
        fi
        ;;
      no)
        HAVE_CSHARPEXEC=
        ac_result="no"
        break
        ;;
    esac
  done
  AC_MSG_RESULT([$ac_result])
  AC_SUBST([MONO_PATH])
  AC_SUBST([MONO_PATH_SEPARATOR])
  AC_SUBST([CLIX_PATH_VAR])
  AC_SUBST([CLIX_PATH])
  AC_SUBST([HAVE_MONO])
  AC_SUBST([HAVE_CLIX])
])
