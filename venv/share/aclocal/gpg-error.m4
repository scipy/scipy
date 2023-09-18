# gpg-error.m4 - autoconf macro to detect libgpg-error.
# Copyright (C) 2002, 2003, 2004, 2011, 2014, 2018, 2020, 2021, 2022
#               g10 Code GmbH
#
# This file is free software; as a special exception the author gives
# unlimited permission to copy and/or distribute it, with or without
# modifications, as long as this notice is preserved.
#
# This file is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY, to the extent permitted by law; without even the
# implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
#
# Last-changed: 2023-04-01

dnl
dnl Find gpg-error-config, for backward compatibility
dnl
dnl _AM_PATH_POSSIBLE_GPG_ERROR_CONFIG
AC_DEFUN([_AM_PATH_POSSIBLE_GPG_ERROR_CONFIG],[dnl
  gpg_error_config_prefix=""
  dnl --with-libgpg-error-prefix=PFX is the preferred name for this option,
  dnl since that is consistent with how our three siblings use the directory/
  dnl package name in --with-$dir_name-prefix=PFX.
  AC_ARG_WITH(libgpg-error-prefix,
              AS_HELP_STRING([--with-libgpg-error-prefix=PFX],
                             [prefix where GPG Error is installed (optional)]),
              [gpg_error_config_prefix="$withval"])

  dnl Accept --with-gpg-error-prefix and make it work the same as
  dnl --with-libgpg-error-prefix above, for backwards compatibility,
  dnl but do not document this old, inconsistently-named option.
  AC_ARG_WITH(gpg-error-prefix,,
              [gpg_error_config_prefix="$withval"])

  if test x"${GPG_ERROR_CONFIG}" = x ; then
     if test x"${gpg_error_config_prefix}" != x ; then
        GPG_ERROR_CONFIG="${gpg_error_config_prefix}/bin/gpg-error-config"
     else
       case "${SYSROOT}" in
         /*)
           if test -x "${SYSROOT}/bin/gpg-error-config" ; then
             GPG_ERROR_CONFIG="${SYSROOT}/bin/gpg-error-config"
           fi
           ;;
         '')
           ;;
          *)
           AC_MSG_WARN([Ignoring \$SYSROOT as it is not an absolute path.])
           ;;
       esac
     fi
  fi

  AC_PATH_PROG(GPG_ERROR_CONFIG, gpg-error-config, no)
])

dnl
dnl Find gpgrt-config, which uses .pc file
dnl (minimum pkg-config functionality, supporting cross build)
dnl
dnl _AM_PATH_GPGRT_CONFIG
AC_DEFUN([_AM_PATH_GPGRT_CONFIG],[dnl
  AC_PATH_PROG(GPGRT_CONFIG, gpgrt-config, no, [$prefix/bin:$PATH])
  if test "$GPGRT_CONFIG" != "no"; then
    # Determine gpgrt_libdir
    #
    # Get the prefix of gpgrt-config assuming it's something like:
    #   <PREFIX>/bin/gpgrt-config
    gpgrt_prefix=${GPGRT_CONFIG%/*/*}
    possible_libdir1=${gpgrt_prefix}/lib
    # Determine by using system libdir-format with CC, it's like:
    #   Normal style: /usr/lib
    #   GNU cross style: /usr/<triplet>/lib
    #   Debian style: /usr/lib/<multiarch-name>
    #   Fedora/openSUSE style: /usr/lib, /usr/lib32 or /usr/lib64
    # It is assumed that CC is specified to the one of host on cross build.
    if libdir_candidates=$(${CC:-cc} -print-search-dirs | \
          sed -n -e "/^libraries/{s/libraries: =//;s/:/\\
/g;p;}"); then
      # From the output of -print-search-dirs, select valid pkgconfig dirs.
      libdir_candidates=$(for dir in $libdir_candidates; do
        if p=$(cd $dir 2>/dev/null && pwd); then
          test -d "$p/pkgconfig" && echo $p;
        fi
      done)

      for possible_libdir0 in $libdir_candidates; do
        # possible_libdir0:
        #   Fallback candidate, the one of system-installed (by $CC)
        #   (/usr/<triplet>/lib, /usr/lib/<multiarch-name> or /usr/lib32)
        # possible_libdir1:
        #   Another candidate, user-locally-installed
        #   (<gpgrt_prefix>/lib)
        # possible_libdir2
        #   Most preferred
        #   (<gpgrt_prefix>/<triplet>/lib,
        #    <gpgrt_prefix>/lib/<multiarch-name> or <gpgrt_prefix>/lib32)
        if test "${possible_libdir0##*/}" = "lib"; then
          possible_prefix0=${possible_libdir0%/lib}
          possible_prefix0_triplet=${possible_prefix0##*/}
          if test -z "$possible_prefix0_triplet"; then
            continue
          fi
          possible_libdir2=${gpgrt_prefix}/$possible_prefix0_triplet/lib
        else
          possible_prefix0=${possible_libdir0%%/lib*}
          possible_libdir2=${gpgrt_prefix}${possible_libdir0#$possible_prefix0}
        fi
        if test -f ${possible_libdir2}/pkgconfig/gpg-error.pc; then
          gpgrt_libdir=${possible_libdir2}
        elif test -f ${possible_libdir1}/pkgconfig/gpg-error.pc; then
          gpgrt_libdir=${possible_libdir1}
        elif test -f ${possible_libdir0}/pkgconfig/gpg-error.pc; then
          gpgrt_libdir=${possible_libdir0}
        fi
        if test -n "$gpgrt_libdir"; then break; fi
      done
    fi
    if test -z "$gpgrt_libdir"; then
      # No valid pkgconfig dir in any of the system directories, fallback
      gpgrt_libdir=${possible_libdir1}
    fi
  else
    unset GPGRT_CONFIG
  fi

  if test -n "$gpgrt_libdir"; then
    GPGRT_CONFIG="$GPGRT_CONFIG --libdir=$gpgrt_libdir"
    if $GPGRT_CONFIG gpg-error >/dev/null 2>&1; then
      GPG_ERROR_CONFIG="$GPGRT_CONFIG gpg-error"
      AC_MSG_NOTICE([Use gpgrt-config with $gpgrt_libdir as gpg-error-config])
      gpg_error_config_version=`$GPG_ERROR_CONFIG --modversion`
    else
      gpg_error_config_version=`$GPG_ERROR_CONFIG --version`
      unset GPGRT_CONFIG
    fi
  elif test "$GPG_ERROR_CONFIG" != "no"; then
    gpg_error_config_version=`$GPG_ERROR_CONFIG --version`
    unset GPGRT_CONFIG
  fi
])

dnl AM_PATH_GPG_ERROR([MINIMUM-VERSION,
dnl                   [ACTION-IF-FOUND [, ACTION-IF-NOT-FOUND ]]])
dnl
dnl Test for libgpg-error and define GPG_ERROR_CFLAGS, GPG_ERROR_LIBS,
dnl GPG_ERROR_MT_CFLAGS, and GPG_ERROR_MT_LIBS.  The _MT_ variants are
dnl used for programs requireing real multi thread support.
dnl
dnl If a prefix option is not used, the config script is first
dnl searched in $SYSROOT/bin and then along $PATH.  If the used
dnl config script does not match the host specification the script
dnl is added to the gpg_config_script_warn variable.
dnl
AC_DEFUN([AM_PATH_GPG_ERROR],[dnl
AC_REQUIRE([AC_CANONICAL_HOST])dnl
AC_REQUIRE([_AM_PATH_POSSIBLE_GPG_ERROR_CONFIG])dnl
AC_REQUIRE([_AM_PATH_GPGRT_CONFIG])dnl
  min_gpg_error_version=ifelse([$1], ,1.33,$1)
  ok=no
  if test "$GPG_ERROR_CONFIG" != "no"; then
    req_major=`echo $min_gpg_error_version | \
               sed 's/\([[0-9]]*\)\.\([[0-9]]*\)/\1/'`
    req_minor=`echo $min_gpg_error_version | \
               sed 's/\([[0-9]]*\)\.\([[0-9]]*\)/\2/'`
    major=`echo $gpg_error_config_version | \
               sed 's/\([[0-9]]*\)\.\([[0-9]]*\).*/\1/'`
    minor=`echo $gpg_error_config_version | \
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
  AC_MSG_CHECKING(for GPG Error - version >= $min_gpg_error_version)
  if test $ok = yes; then
    GPG_ERROR_CFLAGS=`$GPG_ERROR_CONFIG --cflags`
    GPG_ERROR_LIBS=`$GPG_ERROR_CONFIG --libs`
    if test -z "$GPGRT_CONFIG"; then
      GPG_ERROR_MT_CFLAGS=`$GPG_ERROR_CONFIG --mt --cflags 2>/dev/null`
      GPG_ERROR_MT_LIBS=`$GPG_ERROR_CONFIG --mt --libs 2>/dev/null`
    else
      GPG_ERROR_MT_CFLAGS=`$GPG_ERROR_CONFIG --variable=mtcflags 2>/dev/null`
      GPG_ERROR_MT_CFLAGS="$GPG_ERROR_CFLAGS${GPG_ERROR_CFLAGS:+ }$GPG_ERROR_MT_CFLAGS"
      GPG_ERROR_MT_LIBS=`$GPG_ERROR_CONFIG --variable=mtlibs 2>/dev/null`
      GPG_ERROR_MT_LIBS="$GPG_ERROR_LIBS${GPG_ERROR_LIBS:+ }$GPG_ERROR_MT_LIBS"
    fi
    AC_MSG_RESULT([yes ($gpg_error_config_version)])
    ifelse([$2], , :, [$2])
    if test -z "$GPGRT_CONFIG"; then
      gpg_error_config_host=`$GPG_ERROR_CONFIG --host 2>/dev/null || echo none`
    else
      gpg_error_config_host=`$GPG_ERROR_CONFIG --variable=host 2>/dev/null || echo none`
    fi
    if test x"$gpg_error_config_host" != xnone ; then
      if test x"$gpg_error_config_host" != x"$host" ; then
  AC_MSG_WARN([[
***
*** The config script "$GPG_ERROR_CONFIG" was
*** built for $gpg_error_config_host and thus may not match the
*** used host $host.
*** You may want to use the configure option --with-libgpg-error-prefix
*** to specify a matching config script or use \$SYSROOT.
***]])
        gpg_config_script_warn="$gpg_config_script_warn libgpg-error"
      fi
    fi
  else
    GPG_ERROR_CFLAGS=""
    GPG_ERROR_LIBS=""
    GPG_ERROR_MT_CFLAGS=""
    GPG_ERROR_MT_LIBS=""
    AC_MSG_RESULT(no)
    ifelse([$3], , :, [$3])
  fi
  AC_SUBST(GPG_ERROR_CFLAGS)
  AC_SUBST(GPG_ERROR_LIBS)
  AC_SUBST(GPG_ERROR_MT_CFLAGS)
  AC_SUBST(GPG_ERROR_MT_LIBS)
])
