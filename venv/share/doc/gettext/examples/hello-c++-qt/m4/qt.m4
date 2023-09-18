# Copyright (C) 2001 David Johnson
# This file is free software; the author gives unlimited permission to  copy
# and/or distribute it, with or without modifications, as long as this notice
# is preserved.

# FUN_TYPE BOOL
# check for a built-in bool type
# HAVE_BOOL will be defined in the config header

AC_DEFUN([FUN_TYPE_BOOL],
[
  AC_REQUIRE([AC_PROG_CXX])
    
  AC_LANG_PUSH(C++)
  AC_CHECK_TYPE(bool, ac_check_bool=yes, ac_check_bool=no) 
  AC_LANG_POP(C++)
  if test "x$ac_check_bool" = "xyes" ; then
    AC_DEFINE(HAVE_BOOL,,[define if bool is a built-in type])
  fi

  AH_BOTTOM([#ifndef HAVE_BOOL])
  AH_BOTTOM([enum booltyp { false, true }; typedef enum booltyp bool;])
  AH_BOTTOM([#endif])
])# FUN_TYPE_BOOL

# FUN_HEADER_STDCXX
# check for standard ISO C++ headers

AC_DEFUN([FUN_HEADER_STDCXX],
[
  AC_REQUIRE([AC_PROG_CXX])

  AC_LANG_PUSH(C++)

  ac_check_headers=no
  AC_CHECK_HEADER(cstdlib,
      ac_check_headers=yes,
      ac_check_headers=no)
  AC_CHECK_HEADER(cstring,
      ac_check_headers=$ac_check_headers,
      ac_check_headers=no)
  AC_CHECK_HEADER(iostream,
      ac_check_headers=$ac_check_headers,
      ac_check_headers=no)

  AC_LANG_POP(C++)

  if test "x$ac_check_headers" = "xno" ; then
    AC_MSG_ERROR(standard ISO C++ headers not found!)
  fi
])#FUN_HEADER_STDCXX

# FUN_CHECK_PTHREAD
# check for posix pthreads
# sets PTHREAD_LIBS and PTHREAD_CFLAGS
# sets HAVE_PTHREADS in the configuration header

AC_DEFUN([FUN_CHECK_PTHREAD],
[
  AC_REQUIRE([AC_CANONICAL_HOST])
  AC_REQUIRE([AC_PROG_CC])

  PTHREAD_LIBS=""
  PTHREAD_CFLAGS=""

  AC_ARG_ENABLE(threads, AC_HELP_STRING([--enable-threads],
                                [enable the use of the threads [[default=no]]]),
    ac_use_threads=$enableval, ac_use_threads=no)

  if test "x$ac_use_threads" = "xyes" ; then

    AC_CHECK_HEADER(pthread.h, ac_posix_threads=yes, ac_posix_threads=no)

    if test "x$ac_posix_threads" = "xyes" ; then

      AC_MSG_CHECKING([whether ${CC} accepts -pthread])
      ac_cflags_save="$CFLAGS"
      CFLAGS="$CFLAGS -pthread"
      AC_TRY_COMPILE([#include <pthread.h>], [pthread_attr_init(0)],
        ac_cc_pthread=yes, ac_cc_pthread=no)
      CFLAGS="$ac_cflags_save"

      if test "x$ac_cc_pthread" = "xyes" ; then
        AC_MSG_RESULT([yes])
        PTHREAD_CFLAGS="-pthread"
      else
        AC_MSG_RESULT([no])
        ac_thread_library=none

        if test "x$ac_thread_library" = "xnone" ; then
          AC_CHECK_LIB(c_r, pthread_self, ac_thread_library=c_r)
        fi
        if test "x$ac_thread_library" = "xnone" ; then
          AC_CHECK_LIB(pthread, pthread_self, ac_thread_library=pthread)
        fi
        if test "x$ac_thread_library" = "xnone" ; then
          AC_CHECK_LIB(pthreads, pthread_self, ac_thread_library=pthreads)
        fi
        if test "x$ac_thread_library" = "xnone" ; then
          AC_CHECK_LIB(thread, pthread_self, ac_thread_library=thread)
        fi
        if test "x$ac_thread_library" = "xnone" ; then
          AC_CHECK_LIB(gthreads, pthread_self, ac_thread_library=gthreads)
        fi
        if test "x$ac_thread_library" = "xnone" ; then
          AC_CHECK_LIB(c, pthread_self, ac_thread_library=c)
        fi
        if test "x$ac_thread_library" = "xnone" ; then
          ac_use_threads=no
        else
          PTHREAD_LIBS="-l$ac_thread_library"
        fi
      fi
    else
      ac_use_threads=no
    fi
  fi

  if test "x$ac_use_threads" = "xyes" ; then
    AC_DEFINE(HAVE_PTHREAD, 1, [Define if you have POSIX threads])
    case $host_os in
      aix* | freebsd*)
        PTHREAD_CFLAGS="$PTHREAD_CFLAGS -D_THREAD_SAFE"
        ;;
      linux* | solaris*)
        PTHREAD_CFLAGS="$PTHREAD_CFLAGS -D_REENTRANT"
        ;;
      *)
        ;;
    esac
  fi
])#FUN_CHECK_PTHREAD

# FUN_CHECK_QT([qt_min_version],[qt_max_version])
# check for qt headers, libs, progs and compilation
# substs QT_CXXFLAGS, QT_LDFLAGS, and QT_LIBS
# substs QTVERSION, MOC and UIC
# LIBQT, MOC and UIC 'precious' variables

AC_DEFUN([FUN_CHECK_QT],
[
  AC_REQUIRE([AC_PROG_CXX])
  AC_REQUIRE([AC_PATH_X])
  AC_REQUIRE([AC_PATH_XTRA])
  AC_REQUIRE([FUN_CHECK_PTHREAD])

  # some 'precious' variables for configure --help
  AC_ARG_VAR(QTMIN, minimum version of Qt to search for  e.g. 220)
  AC_ARG_VAR(QTMAX, maximum version of Qt to search for  e.g. 399)
  AC_ARG_VAR(LIBQT, library flag for the Qt libary  e.g. -lqt)
  AC_ARG_VAR(MOC, QT meta object compiler command)
  AC_ARG_VAR(UIC, Qt UI compiler command)

  AC_CACHE_SAVE

  AC_MSG_NOTICE([checking for Qt])

  # process our args
  if test -z "$1" ; then
    qt_min_version=0
  else
    qt_min_version=$1
  fi
  if test -z "$2" ; then
    qt_max_version=9999
  else
    qt_max_version=$2
  fi
  # adjust for user preferences
  if test "x$QTMIN" != "x" ; then
    if expr $QTMIN '>' $qt_min_version > /dev/null ; then
      qt_min_version=$QTMIN;
    fi
  fi
  if test "x$QTMAX" != "x" ; then
    if expr $QTMAX '<' $qt_max_version > /dev/null ; then
      qt_max_version=$QTMAX;
    fi
  fi

  # set up our configuration options
  qt_dir=""
  qt_includes=""
  qt_libraries=""
  qt_programs=""
  AC_ARG_WITH([qt_dir], AC_HELP_STRING([--with-qt-dir=DIR],
                        [where the Qt package is installed]),
    [ qt_dir="$withval"
      qt_includes="$withval"/include
      qt_libraries="$withval"/lib
      qt_programs="$withval"/bin
    ])
  AC_ARG_WITH([qt_includes], AC_HELP_STRING([--with-qt-includes=DIR],
                             [where the Qt includes are installed]),
    [qt_includes="$withval"])
  AC_ARG_WITH([qt_libraries], AC_HELP_STRING([--with-qt-libraries=DIR],
                              [where the Qt libraries are installed]),
    [qt_libraries="$withval"])
  AC_ARG_WITH([qt_programs], AC_HELP_STRING([--with-qt-programs=DIR],
                             [where the Qt programs are installed]),
    [qt_programs="$withval"])

  QTVERSION="000"

  FUN_QT_HEADERS

  # check for a traditional qt installation tree
  if ls $qt_includes/../lib/libqt* > /dev/null 2> /dev/null; then
    qt_dir="`echo $qt_includes | sed s,'/include',,`"
    qt_libraries="$qt_dir/lib"
    qt_programs="$qt_dir/bin"
  fi

  FUN_QT_LIBRARIES
  FUN_QT_PROGRAMS
  FUN_QT_COMPILE

  AC_MSG_NOTICE([Found Qt version $QTVERSION])

  AC_SUBST(QTVERSION)
  AC_SUBST(MOC)
  AC_SUBST(UIC)
  QT_CXXFLAGS="-I$qt_includes"
  AC_SUBST(QT_CXXFLAGS)
  QT_LDFLAGS="-L$qt_libraries"
  AC_SUBST(QT_LDFLAGS)
  QT_LIBS="$LIBQT"
  AC_SUBST(QT_LIBS)
])#FUN_CHECK_QT

# FUN_QT_HEADERS
# helper function for FUN_CHECK_QT
# check for qt headers in standard locations

AC_DEFUN([FUN_QT_HEADERS],
[
  AC_MSG_CHECKING([for Qt headers])

  if test "x$qt_includes" = "x" ; then
    # look in standard locations
    qt_found_dirs=""
    qt_include_dirs="
      $QTDIR
      /usr/include
      /usr/local/include
      /usr/X11R6/include
      `ls -dr /usr/include/qt* 2>/dev/null`
      `ls -dr /usr/local/include/qt* 2>/dev/null`
      `ls -dr /usr/X11R6/include/qt* 2>/dev/null`
      `ls -dr /usr/lib/qt*/include 2>/dev/null`
      `ls -dr /usr/local/lib/qt*/include 2>/dev/null`
      `ls -dr /usr/X11R6/lib/qt*/include 2>/dev/null`
      `ls -dr /usr/local/qt*/include 2>/dev/null`
      `ls -dr /opt/qt*/include 2>/dev/null` "
    for n in $qt_include_dirs ; do
      if test -r "$n/qglobal.h"; then
        qt_found_dirs="$qt_found_dirs $n"
      fi
    done

    # find the latest version between min_version and max_version
    qt_prev_version=$qt_min_version
    qt_found_version=""
    for n in $qt_found_dirs ; do
      qt_current_version=`grep -w '#define QT_VERSION' $n/qglobal.h |
        sed s/'#define QT_VERSION'//`
      if expr $qt_current_version '>=' $qt_prev_version > /dev/null ; then
        if expr $qt_current_version '<=' $qt_max_version > /dev/null ; then
          qt_includes=$n
          qt_prev_version=$qt_current_version
        fi
      fi
    done
  fi

  if test "x$qt_includes" = "x" ; then
    AC_MSG_RESULT([no])
    AC_MSG_ERROR([cannot find correct Qt headers!])
  else
    dnl TODO need to strip out white space
    QTVERSION=$qt_prev_version;
    AC_MSG_RESULT([$qt_includes])
  fi
])#FUN_QT_HEADERS

# FUN_QT_LIBRARIES
# helper function for FUN_CHECK_QT
# check for qt libs in standard locations

AC_DEFUN([FUN_QT_LIBRARIES],
[
  AC_REQUIRE([FUN_QT_HEADERS])

  AC_MSG_CHECKING([for Qt libraries])

  # check which library to look for
  if test -z "$LIBQT" ; then
    if test "x$ac_use_threads" = "xyes" ; then
      LIBQT="-lqt-mt"
    else
      LIBQT="-lqt"
    fi
  fi

  lib_qt=`echo $LIBQT | sed s/'-l'//`

  if test "x$qt_libraries" = "x" ; then
    # see if it is relative to the includes
    qt_tree="$qt_includes"
    while test "x$qt_tree" != "x" ; do
      # first go around will fail...
      if ls $qt_tree/lib/libqt* > /dev/null 2> /dev/null ; then
        qt_libraries=$qt_tree/lib
        break
      else
        # lop off tail of path
        dnl not as portable as it should be...
        qt_tree="`dirname $qt_tree`"
      fi
    done
  fi  

  if test "x$qt_libraries" = "x" ; then
    AC_MSG_RESULT([no])
    AC_MSG_ERROR([cannot find Qt libraries!])
  else
    # check that we're looking at the right library
    if ls $qt_libraries/lib$lib_qt.* > /dev/null 2> /dev/null ; then
      AC_MSG_RESULT([$qt_libraries])
    else
      AC_MSG_RESULT([no])
      if test "x$ac_use_threads" = "xyes" ; then
        AC_MSG_ERROR([cannot find the threaded Qt library in $qt_libraries!])
      else
        AC_MSG_ERROR([cannot find the non-threaded Qt library in $qt_libraries!])
      fi
    fi
  fi
])#FUN_QT_LIBRARIES

# FUN_QT_PROGRAMS
# helper function for FUN_CHECK_QT
# searches for moc and uic

AC_DEFUN([FUN_QT_PROGRAMS],
[
  AC_REQUIRE([FUN_QT_LIBRARIES])

  AC_MSG_CHECKING([for Qt utilities])

  if test "x$q_programs" = "x" ; then
    # see if it is relative to the libraries
    qt_tree="$qt_libraries"
    while test "x$qt_tree" != "x" ; do
      # first go around will fail
      if ls $qt_tree/bin/moc* > /dev/null 2> /dev/null ; then
        qt_programs=$qt_tree/bin
        break
      else
        # lop off tail of path
        dnl not as portable as it should be...
        qt_tree="`dirname $qt_tree`"
      fi
    done
    # if we haven't found the progs, there's not much more we can do
  fi  

  if test "x$qt_programs" = "x" ; then
    AC_MSG_RESULT([no])
    AC_MSG_ERROR([cannot find Qt utilities!])
  else
    AC_MSG_RESULT([$qt_programs])
    # find the right moc
    if test -z "$MOC" ; then
      AC_CHECK_PROG(MOC, moc, moc)
      if test "x$MOC" = "x" ; then
        # could be renamed to avoid clashes
        if ls $qt_programs/moc > /dev/null 2> /dev/null ; then
          MOC="$qt_programs/moc"
        else
          if expr "$QTVERSION" '>=' "200" > /dev/null ; then
            if ls $qt_programs/moc2 > /dev/null 2> /dev/null ; then
              MOC="$qt_programs/moc2"
            fi
          else
            if expr "$QTVERSION" '>=' "300" > /dev/null ; then
              if $qt_programs/moc3 > /dev/null 2> /dev/null ; then
                MOC="$qt_programs/moc3"
              fi
            fi
          fi
        fi
      fi
      if test "x$MOC" = "x" ; then
        AC_MSG_RESULT([no])
        AC_MSG_ERROR([cannot find Qt meta object compiler!])
      fi
    fi

    # find the right uic
    if expr "$QTVERSION" '>=' "220" > /dev/null ; then
      if test -z "$UIC" ; then
        AC_CHECK_PROG(UIC, uic, uic)
        if test "x$UIC" = "x" ; then
          # could be renamed to avoid clashes
          if ls $qt_programs/uic > /dev/null 2> /dev/null ; then
            UIC="$qt_programs/uic"
          else
            if expr "$QTVERSION" '>=' "300" > /dev/null ; then
              if ls $qt_programs/uic3 > /dev/null 2> /dev/null ; then
                UIC="$qt_programs/uic3"
              fi
            fi
          fi
        fi
      fi
    else
      # if uic is important to the build, change this
      UIC=""
    fi
  fi
])#FUN_QT_PROGRAMS

# FUN_QT_COMPILE
# helper function for FUN_CHECK_QT
# compile a simple qt program

AC_DEFUN([FUN_QT_COMPILE],
[
  AC_REQUIRE([FUN_QT_HEADERS])
  AC_REQUIRE([FUN_QT_LIBRARIES])
  AC_REQUIRE([FUN_QT_PROGRAMS])

  AC_MSG_CHECKING([whether a simple Qt program compiles])

  AC_LANG_PUSH(C++)

  ac_cxxflags_save="$CXXFLAGS"
  ac_ldflags_save="$LDFLAGS"
  ac_libs_save="$LIBS"
  CXXFLAGS="$CXXFLAGS $PTHREAD_CFLAGS -I$qt_includes $X_CFLAGS $all_includes"
  LDFLAGS="$LDFLAGS -L$qt_libraries $X_LIBS "
  LIBS="$LIBS $PTHREAD_LIBS $X_PRE_LIBS $X_EXTRA_LIBS -lXext -lX11 $LIBQT"

  AC_TRY_LINK([
    #include <qglobal.h>
    #include <qmessagebox.h>
    #include <qstring.h>],
    [QString s = "hello world";
    QMessageBox::information(0, s, "no he is not");
    return 0;],
  qt_compile=yes, qt_compile=no)

  CXXFLAGS="$ac_cxxflags_save"
  LDFLAGS="$ac_ldflags_save"
  LIBS="$ac_libs_save"

  AC_LANG_POP(C++)

  if test "x$qt_compile" = "xyes" ; then
    AC_MSG_RESULT([yes])
  else
    AC_MSG_RESULT([no])
    AC_MSG_ERROR([cannot compile a Qt program!])
  fi
])#FUN_QT_COMPILE
