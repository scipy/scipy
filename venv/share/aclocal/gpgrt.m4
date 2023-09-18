# gpgrt.m4 - autoconf macro to detect libgpgrt
# Copyright (C) 2002, 2003, 2004, 2011, 2014, 2017, 2018 g10 Code GmbH
#
# This file is free software; as a special exception the author gives
# unlimited permission to copy and/or distribute it, with or without
# modifications, as long as this notice is preserved.
#
# This file is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY, to the extent permitted by law; without even the
# implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
# SPDX-License-Identifier: FSFULLR
#
# Last-changed: 2018-11-13
# Note: This is a kind of duplicate of gpg-error.m4 which uses the
# future name of libgpg-error to prepare for a smooth migration in
# some distant time.

dnl AM_PATH_GPGRT([MINIMUM-VERSION,
dnl               [ACTION-IF-FOUND [, ACTION-IF-NOT-FOUND ]]])
dnl
dnl Test for libgpgrt and define GPGRT_CFLAGS, GPGRT_LIBS,
dnl GPGRT_MT_CFLAGS, and GPGRT_MT_LIBS.  The _MT_ variants are
dnl used for programs requiring real multi thread support.
dnl
AC_DEFUN([AM_PATH_GPGRT],
[ AC_REQUIRE([AC_CANONICAL_HOST])
  if test "$prefix" = NONE ; then
    prefix_option_expanded=/usr/local
  else
    prefix_option_expanded="$prefix"
  fi
  if test "$exec_prefix" = NONE ; then
    exec_prefix_option_expanded=$prefix_option_expanded
  else
    exec_prefix_option_expanded=$(prefix=$prefix_option_expanded eval echo $exec_prefix)
  fi
  libdir_option_expanded=$(prefix=$prefix_option_expanded exec_prefix=$exec_prefix_option_expanded eval echo $libdir)

  if test -f $libdir_option_expanded/pkgconfig/gpg-error.pc; then
    gpgrt_libdir=$libdir_option_expanded
  else
    if crt1_path=$(${CC:-cc} -print-file-name=crt1.o 2>/dev/null); then
      if possible_libdir=$(cd ${crt1_path%/*} && pwd 2>/dev/null); then
        if test -f $possible_libdir/pkgconfig/gpg-error.pc; then
          gpgrt_libdir=$possible_libdir
        fi
      fi
    fi
  fi

  if test -n "$gpgrt_libdir"; then
    AC_PATH_PROG(GPGRT_CONFIG, gpgrt-config, no)
    if test "$GPGRT_CONFIG" != "no"; then
      GPGRT_CONFIG="$GPGRT_CONFIG --libdir=$gpgrt_libdir"
    fi
  fi
  min_gpgrt_version=ifelse([$1], ,1.33,$1)
  AC_MSG_CHECKING(for GPG Runtime - version >= $min_gpgrt_version)
  ok=no
  if test x"$GPGRT_CONFIG" != x -a "$GPGRT_CONFIG" != "no" ; then
    req_major=`echo $min_gpgrt_version | \
               sed 's/\([[0-9]]*\)\.\([[0-9]]*\)/\1/'`
    req_minor=`echo $min_gpgrt_version | \
               sed 's/\([[0-9]]*\)\.\([[0-9]]*\)/\2/'`
    gpgrt_config_version=`$GPGRT_CONFIG --version`
    major=`echo $gpgrt_config_version | \
               sed 's/\([[0-9]]*\)\.\([[0-9]]*\).*/\1/'`
    minor=`echo $gpgrt_config_version | \
               sed 's/\([[0-9]]*\)\.\([[0-9]]*\).*/\2/'`
    if test "$major" -gt "$req_major"; then
        ok=yes
    else
        if test "$major" -eq "$req_major"; then
            if test "$minor" -ge "$req_minor"; then
               ok=yes
            fi
        fi
    fi
  fi
  if test $ok = yes; then
    GPGRT_CFLAGS=`$GPGRT_CONFIG --cflags`
    GPGRT_LIBS=`$GPGRT_CONFIG --libs`
    GPGRT_MT_CFLAGS=`$GPGRT_CONFIG --variable=mtcflags 2>/dev/null`
    GPGRT_MT_CFLAGS="$GPGRT_CFLAGS${GPGRT_CFLAGS:+ }$GPGRT_MT_CFLAGS"
    GPGRT_MT_LIBS=`$GPGRT_CONFIG --variable=mtlibs 2>/dev/null`
    GPGRT_MT_LIBS="$GPGRT_LIBS${GPGRT_LIBS:+ }$GPGRT_MT_LIBS"
    AC_MSG_RESULT([yes ($gpgrt_config_version)])
    ifelse([$2], , :, [$2])
    gpgrt_config_host=`$GPGRT_CONFIG --variable=host 2>/dev/null || echo none`
    if test x"$gpgrt_config_host" != xnone ; then
      if test x"$gpgrt_config_host" != x"$host" ; then
  AC_MSG_WARN([[
***
*** The config script "$GPGRT_CONFIG" is for $gpgrt_config_host
*** and thus may not match the used host $host.
***]])
        gpg_config_script_warn="$gpg_config_script_warn libgpgrt"
      fi
    fi
  else
    GPGRT_CFLAGS=""
    GPGRT_LIBS=""
    GPGRT_MT_CFLAGS=""
    GPGRT_MT_LIBS=""
    AC_MSG_RESULT(no)
    ifelse([$3], , :, [$3])
  fi
  AC_SUBST(GPGRT_CFLAGS)
  AC_SUBST(GPGRT_LIBS)
  AC_SUBST(GPGRT_MT_CFLAGS)
  AC_SUBST(GPGRT_MT_LIBS)
])
