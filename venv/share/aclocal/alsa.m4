dnl Configure Paths for Alsa
dnl Some modifications by Richard Boulton <richard-alsa@tartarus.org>
dnl Christopher Lansdown <lansdoct@cs.alfred.edu>
dnl Jaroslav Kysela <perex@perex.cz>
dnl Last modification: $Id: alsa.m4,v 1.24 2004/09/15 18:48:07 tiwai Exp $
dnl
dnl AM_PATH_ALSA([MINIMUM-VERSION [, ACTION-IF-FOUND [, ACTION-IF-NOT-FOUND]]])
dnl Test for libasound, and define ALSA_CFLAGS, ALSA_LIBS and
dnl ALSA_TOPOLOGY_LIBS as appropriate.
dnl
dnl enables arguments --with-alsa-prefix=
dnl                   --with-alsa-inc-prefix=
dnl                   --disable-alsatest
dnl
dnl For backwards compatibility, if ACTION_IF_NOT_FOUND is not specified,
dnl and the alsa libraries are not found, a fatal AC_MSG_ERROR() will result.
dnl

AC_DEFUN([AM_PATH_ALSA],
[dnl Save the original CFLAGS, LDFLAGS, and LIBS
alsa_save_CFLAGS="$CFLAGS"
alsa_save_LDFLAGS="$LDFLAGS"
alsa_save_LIBS="$LIBS"
alsa_found=yes
alsa_topology_found=no

dnl
dnl Get the cflags and libraries for alsa
dnl
AC_ARG_WITH(alsa-prefix,
  AS_HELP_STRING([--with-alsa-prefix=PFX], [Prefix where Alsa library is installed(optional)]),
  [alsa_prefix="$withval"], [alsa_prefix=""])

AC_ARG_WITH(alsa-inc-prefix,
  AS_HELP_STRING([--with-alsa-inc-prefix=PFX], [Prefix where include libraries are (optional)]),
  [alsa_inc_prefix="$withval"], [alsa_inc_prefix=""])

AC_ARG_ENABLE(alsa-topology,
  AS_HELP_STRING([--enable-alsatopology], [Force to use the Alsa topology library]),
  [enable_atopology="$enableval"],
  [enable_atopology=no])

AC_ARG_ENABLE(alsatest,
  AS_HELP_STRING([--disable-alsatest], [Do not try to compile and run a test Alsa program]),
  [enable_alsatest="$enableval"],
  [enable_alsatest=yes])

dnl Add any special include directories
AC_MSG_CHECKING(for ALSA CFLAGS)
if test "$alsa_inc_prefix" != "" ; then
	ALSA_CFLAGS="$ALSA_CFLAGS -I$alsa_inc_prefix"
	CFLAGS="$CFLAGS -I$alsa_inc_prefix"
fi
AC_MSG_RESULT($ALSA_CFLAGS)

AC_CHECK_LIB(c, dlopen, LIBDL="", [AC_CHECK_LIB(dl, dlopen, LIBDL="-ldl")])

dnl add any special lib dirs
AC_MSG_CHECKING(for ALSA LDFLAGS)
if test "$alsa_prefix" != "" ; then
	ALSA_LIBS="$ALSA_LIBS -L$alsa_prefix"
	LDFLAGS="$LDFLAGS $ALSA_LIBS"
fi

dnl add the alsa library
ALSA_LIBS="$ALSA_LIBS -lasound -lm $LIBDL -lpthread"
LIBS="$ALSA_LIBS $LIBS"
AC_MSG_RESULT($ALSA_LIBS)

dnl Check for a working version of libasound that is of the right version.
if test "x$enable_alsatest" = "xyes"; then

AC_MSG_CHECKING([required libasound headers version])
min_alsa_version=ifelse([$1], , 0.1.1, $1)
no_alsa=""
    alsa_min_major_version=`echo $min_alsa_version | \
           sed 's/\([[0-9]]*\).\([[0-9]]*\).\([[0-9]]*\)/\1/'`
    alsa_min_minor_version=`echo $min_alsa_version | \
           sed 's/\([[0-9]]*\).\([[0-9]]*\).\([[0-9]]*\)/\2/'`
    alsa_min_micro_version=`echo $min_alsa_version | \
           sed 's/\([[0-9]]*\).\([[0-9]]*\).\([[0-9]]*\)/\3/'`
AC_MSG_RESULT($alsa_min_major_version.$alsa_min_minor_version.$alsa_min_micro_version)

AC_LANG_PUSH([C])
AC_MSG_CHECKING([for libasound headers version >= $alsa_min_major_version.$alsa_min_minor_version.$alsa_min_micro_version ($min_alsa_version)])
AC_COMPILE_IFELSE([AC_LANG_PROGRAM([[
#include <alsa/asoundlib.h>
]], [[
/* ensure backward compatibility */
#if !defined(SND_LIB_MAJOR) && defined(SOUNDLIB_VERSION_MAJOR)
#define SND_LIB_MAJOR SOUNDLIB_VERSION_MAJOR
#endif
#if !defined(SND_LIB_MINOR) && defined(SOUNDLIB_VERSION_MINOR)
#define SND_LIB_MINOR SOUNDLIB_VERSION_MINOR
#endif
#if !defined(SND_LIB_SUBMINOR) && defined(SOUNDLIB_VERSION_SUBMINOR)
#define SND_LIB_SUBMINOR SOUNDLIB_VERSION_SUBMINOR
#endif

#  if(SND_LIB_MAJOR > $alsa_min_major_version)
  exit(0);
#  else
#    if(SND_LIB_MAJOR < $alsa_min_major_version)
#       error not present
#    endif

#   if(SND_LIB_MINOR > $alsa_min_minor_version)
  exit(0);
#   else
#     if(SND_LIB_MINOR < $alsa_min_minor_version)
#          error not present
#      endif

#      if(SND_LIB_SUBMINOR < $alsa_min_micro_version)
#        error not present
#      endif
#    endif
#  endif
exit(0);
]])],
  [AC_MSG_RESULT(found.)],
  [AC_MSG_RESULT(not present.)
   ifelse([$3], , [AC_MSG_ERROR(Sufficiently new version of libasound not found.)])
   alsa_found=no]
)
AC_LANG_POP([C])

AC_LANG_PUSH([C])
AC_MSG_CHECKING([for libatopology (sound headers version > 1.1.9)])
AC_COMPILE_IFELSE([AC_LANG_PROGRAM([[
#include <alsa/asoundlib.h>
#include <alsa/topology.h>
]], [[
/* ensure backward compatibility */
#if !defined(SND_LIB_VERSION)
#define SND_LIB_VERSION 0
#endif
#if SND_LIB_VERSION > 0x00010109
  exit(0);
#else
# error not present
#endif
exit(0);
]])],
  [AC_MSG_RESULT(yes)
   enable_atopology="yes"],
  [AC_MSG_RESULT(no)]
)
AC_LANG_POP([C])
fi

dnl Now that we know that we have the right version, let's see if we have the library and not just the headers.
if test "x$enable_alsatest" = "xyes"; then
AC_CHECK_LIB([asound], [snd_ctl_open],,
	[ifelse([$3], , [AC_MSG_ERROR(No linkable libasound was found.)])
	 alsa_found=no]
)
if test "x$enable_atopology" = "xyes"; then
alsa_topology_found=yes
alsa_save_LIBS2="$LIBS"
AC_CHECK_LIB([atopology], [snd_tplg_new],,
	[ifelse([$3], , [AC_MSG_ERROR(No linkable libatopology was found.)])
	 alsa_topology_found=no,
]
)
LIBS="$alsa_save_LIBS2"
fi
else
if test "x$enable_atopology" = "xyes"; then
  alsa_topology_found=yes
fi
fi

if test "x$alsa_found" = "xyes" ; then
   ifelse([$2], , :, [$2])
   LIBS=`echo $LIBS | sed 's/-lasound//g'`
   LIBS=`echo $LIBS | sed 's/  //'`
   LIBS="-lasound $LIBS"
fi
if test "x$alsa_found" = "xno" ; then
   ifelse([$3], , :, [$3])
   CFLAGS="$alsa_save_CFLAGS"
   LDFLAGS="$alsa_save_LDFLAGS"
   LIBS="$alsa_save_LIBS"
   ALSA_CFLAGS=""
   ALSA_LIBS=""
   ALSA_TOPOLOGY_LIBS=""
fi

dnl add the alsa topology library; must be at the end
AC_MSG_CHECKING(for ALSA topology LDFLAGS)
if test "x$alsa_topology_found" = "xyes"; then
  ALSA_TOPOLOGY_LIBS="$ALSA_TOPOLOGY_LIBS -latopology"
fi
AC_MSG_RESULT($ALSA_TOPOLOGY_LIBS)

dnl That should be it.  Now just export out symbols:
AC_SUBST(ALSA_CFLAGS)
AC_SUBST(ALSA_LIBS)
AC_SUBST(ALSA_TOPOLOGY_LIBS)
])
