# Configure paths for GTK--
# Erik Andersen	30 May 1998
# Modified by Tero Pulkkinen (added the compiler checks... I hope they work..)
# Modified by Thomas Langen 16 Jan 2000 (corrected CXXFLAGS)

dnl Test for GTKMM, and define GTKMM_CFLAGS and GTKMM_LIBS
dnl   to be used as follows:
dnl AM_PATH_GTKMM([MINIMUM-VERSION, [ACTION-IF-FOUND [, ACTION-IF-NOT-FOUND]]])
dnl
AC_DEFUN([AM_PATH_GTKMM],
[dnl 
dnl Get the cflags and libraries from the gtkmm-config script
dnl
AC_ARG_WITH(gtkmm-prefix,[  --with-gtkmm-prefix=PREFIX
                          Prefix where GTK-- is installed (optional)],
            gtkmm_config_prefix="$withval", gtkmm_config_prefix="")
AC_ARG_WITH(gtkmm-exec-prefix,[  --with-gtkmm-exec-prefix=PREFIX
                          Exec prefix where GTK-- is installed (optional)],
            gtkmm_config_exec_prefix="$withval", gtkmm_config_exec_prefix="")
AC_ARG_ENABLE(gtkmmtest, [  --disable-gtkmmtest     Do not try to compile and run a test GTK-- program],
		    , enable_gtkmmtest=yes)

  if test x$gtkmm_config_exec_prefix != x ; then
     gtkmm_config_args="$gtkmm_config_args --exec-prefix=$gtkmm_config_exec_prefix"
     if test x${GTKMM_CONFIG+set} != xset ; then
        GTKMM_CONFIG=$gtkmm_config_exec_prefix/bin/gtkmm-config
     fi
  fi
  if test x$gtkmm_config_prefix != x ; then
     gtkmm_config_args="$gtkmm_config_args --prefix=$gtkmm_config_prefix"
     if test x${GTKMM_CONFIG+set} != xset ; then
        GTKMM_CONFIG=$gtkmm_config_prefix/bin/gtkmm-config
     fi
  fi

  AC_PATH_PROG(GTKMM_CONFIG, gtkmm-config, no)
  min_gtkmm_version=ifelse([$1], ,0.10.0,$1)

  AC_MSG_CHECKING(for GTK-- - version >= $min_gtkmm_version)
  AC_LANG_SAVE
  no_gtkmm=""
  if test "$GTKMM_CONFIG" = "no" ; then
    no_gtkmm=yes
  else
    AC_LANG_CPLUSPLUS

    GTKMM_CFLAGS=`$GTKMM_CONFIG $gtkmm_config_args --cflags`
    GTKMM_LIBS=`$GTKMM_CONFIG $gtkmm_config_args --libs`
    gtkmm_config_major_version=`$GTKMM_CONFIG $gtkmm_config_args --version | \
           sed 's/\([[0-9]]*\).\([[0-9]]*\).\([[0-9]]*\)/\1/'`
    gtkmm_config_minor_version=`$GTKMM_CONFIG $gtkmm_config_args --version | \
           sed 's/\([[0-9]]*\).\([[0-9]]*\).\([[0-9]]*\)/\2/'`
    gtkmm_config_micro_version=`$GTKMM_CONFIG $gtkmm_config_args --version | \
           sed 's/\([[0-9]]*\).\([[0-9]]*\).\([[0-9]]*\)/\3/'`
    if test "x$enable_gtkmmtest" = "xyes" ; then
      ac_save_CXXFLAGS="$CXXFLAGS"
      ac_save_LIBS="$LIBS"
      CXXFLAGS="$CXXFLAGS $GTKMM_CFLAGS"
      LIBS="$LIBS $GTKMM_LIBS"
dnl
dnl Now check if the installed GTK-- is sufficiently new. (Also sanity
dnl checks the results of gtkmm-config to some extent
dnl
      rm -f conf.gtkmmtest
      AC_TRY_RUN([
#include <gtk--.h>
#include <stdio.h>
#include <stdlib.h>

int 
main ()
{
  int major, minor, micro;
  char *tmp_version;

  system ("touch conf.gtkmmtest");

  /* HP/UX 0 (%@#!) writes to sscanf strings */
  tmp_version = g_strdup("$min_gtkmm_version");
  if (sscanf(tmp_version, "%d.%d.%d", &major, &minor, &micro) != 3) {
     printf("%s, bad version string\n", "$min_gtkmm_version");
     exit(1);
   }

  if ((gtkmm_major_version != $gtkmm_config_major_version) ||
      (gtkmm_minor_version != $gtkmm_config_minor_version) ||
      (gtkmm_micro_version != $gtkmm_config_micro_version))
    {
      printf("\n*** 'gtkmm-config --version' returned %d.%d.%d, but GTK-- (%d.%d.%d)\n", 
             $gtkmm_config_major_version, $gtkmm_config_minor_version, $gtkmm_config_micro_version,
             gtkmm_major_version, gtkmm_minor_version, gtkmm_micro_version);
      printf ("*** was found! If gtkmm-config was correct, then it is best\n");
      printf ("*** to remove the old version of GTK--. You may also be able to fix the error\n");
      printf("*** by modifying your LD_LIBRARY_PATH enviroment variable, or by editing\n");
      printf("*** /etc/ld.so.conf. Make sure you have run ldconfig if that is\n");
      printf("*** required on your system.\n");
      printf("*** If gtkmm-config was wrong, set the environment variable GTKMM_CONFIG\n");
      printf("*** to point to the correct copy of gtkmm-config, and remove the file config.cache\n");
      printf("*** before re-running configure\n");
    } 
/* GTK-- does not have the GTKMM_*_VERSION constants */
/* 
  else if ((gtkmm_major_version != GTKMM_MAJOR_VERSION) ||
	   (gtkmm_minor_version != GTKMM_MINOR_VERSION) ||
           (gtkmm_micro_version != GTKMM_MICRO_VERSION))
    {
      printf("*** GTK-- header files (version %d.%d.%d) do not match\n",
	     GTKMM_MAJOR_VERSION, GTKMM_MINOR_VERSION, GTKMM_MICRO_VERSION);
      printf("*** library (version %d.%d.%d)\n",
	     gtkmm_major_version, gtkmm_minor_version, gtkmm_micro_version);
    }
*/
  else
    {
      if ((gtkmm_major_version > major) ||
        ((gtkmm_major_version == major) && (gtkmm_minor_version > minor)) ||
        ((gtkmm_major_version == major) && (gtkmm_minor_version == minor) && (gtkmm_micro_version >= micro)))
      {
        return 0;
       }
     else
      {
        printf("\n*** An old version of GTK-- (%d.%d.%d) was found.\n",
               gtkmm_major_version, gtkmm_minor_version, gtkmm_micro_version);
        printf("*** You need a version of GTK-- newer than %d.%d.%d. The latest version of\n",
	       major, minor, micro);
        printf("*** GTK-- is always available from ftp://ftp.gtk.org.\n");
        printf("***\n");
        printf("*** If you have already installed a sufficiently new version, this error\n");
        printf("*** probably means that the wrong copy of the gtkmm-config shell script is\n");
        printf("*** being found. The easiest way to fix this is to remove the old version\n");
        printf("*** of GTK--, but you can also set the GTKMM_CONFIG environment to point to the\n");
        printf("*** correct copy of gtkmm-config. (In this case, you will have to\n");
        printf("*** modify your LD_LIBRARY_PATH enviroment variable, or edit /etc/ld.so.conf\n");
        printf("*** so that the correct libraries are found at run-time))\n");
      }
    }
  return 1;
}
],, no_gtkmm=yes,[echo $ac_n "cross compiling; assumed OK... $ac_c"])
       CXXFLAGS="$ac_save_CXXFLAGS"
       LIBS="$ac_save_LIBS"
     fi
  fi
  if test "x$no_gtkmm" = x ; then
     AC_MSG_RESULT(yes)
     ifelse([$2], , :, [$2])     
  else
     AC_MSG_RESULT(no)
     if test "$GTKMM_CONFIG" = "no" ; then
       echo "*** The gtkmm-config script installed by GTK-- could not be found"
       echo "*** If GTK-- was installed in PREFIX, make sure PREFIX/bin is in"
       echo "*** your path, or set the GTKMM_CONFIG environment variable to the"
       echo "*** full path to gtkmm-config."
       echo "*** The gtkmm-config script was not available in GTK-- versions"
       echo "*** prior to 0.9.12. Perhaps you need to update your installed"
       echo "*** version to 0.9.12 or later"
     else
       if test -f conf.gtkmmtest ; then
        :
       else
          echo "*** Could not run GTK-- test program, checking why..."
          CXXFLAGS="$CXXFLAGS $GTKMM_CFLAGS"
          LIBS="$LIBS $GTKMM_LIBS"
          AC_TRY_LINK([
#include <gtk--.h>
#include <stdio.h>
],      [ return ((gtkmm_major_version) || (gtkmm_minor_version) || (gtkmm_micro_version)); ],
        [ echo "*** The test program compiled, but did not run. This usually means"
          echo "*** that the run-time linker is not finding GTK-- or finding the wrong"
          echo "*** version of GTK--. If it is not finding GTK--, you'll need to set your"
          echo "*** LD_LIBRARY_PATH environment variable, or edit /etc/ld.so.conf to point"
          echo "*** to the installed location  Also, make sure you have run ldconfig if that"
          echo "*** is required on your system"
	  echo "***"
          echo "*** If you have an old version installed, it is best to remove it, although"
          echo "*** you may also be able to get things to work by modifying LD_LIBRARY_PATH" ],
        [ echo "*** The test program failed to compile or link. See the file config.log for the"
          echo "*** exact error that occured. This usually means GTK-- was incorrectly installed"
          echo "*** or that you have moved GTK-- since it was installed. In the latter case, you"
          echo "*** may want to edit the gtkmm-config script: $GTKMM_CONFIG" ])
          CXXFLAGS="$ac_save_CXXFLAGS"
          LIBS="$ac_save_LIBS"
       fi
     fi
     GTKMM_CFLAGS=""
     GTKMM_LIBS=""
     ifelse([$3], , :, [$3])
  fi
  AC_LANG_RESTORE
  AC_SUBST(GTKMM_CFLAGS)
  AC_SUBST(GTKMM_LIBS)
  rm -f conf.gtkmmtest
])

