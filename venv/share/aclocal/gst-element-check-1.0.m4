dnl Perform a check for a GStreamer element using gst-inspect-x.y
dnl
dnl GST_ELEMENT_CHECK(ELEMENT-NAME, MIN-VERSION, ACTION-IF-FOUND, ACTION-IF-NOT-FOUND)
dnl
dnl ELEMENT-NAME        : element factory name (mandatory)
dnl MIN-VERSION         : minimum version required, e.g. 1.0 or 1.0.5 (mandatory)
dnl ACTION-IF_FOUND     : action if element exists and is of the desired version
dnl ACTION-IF-NOT-FOUND : action if element does not exist or is too old
dnl
dnl gstapiversion=`echo $2 | tr '.' '\n' | head -n 2 | tr '\n' '.' | sed 's/\.$//'`

AC_DEFUN([GST_ELEMENT_CHECK],
[
  AC_REQUIRE([PKG_PROG_PKG_CONFIG])

  gstapiversion=`echo "$2" | while IFS=. read a b; do echo "$a.0"; done`
  gsttoolsdir=`$PKG_CONFIG --variable=toolsdir gstreamer-$gstapiversion`
  if test "x$gsttoolsdir" != "x"; then
    gstinspect="$gsttoolsdir/gst-inspect-$gstapiversion"
    AC_MSG_CHECKING(GStreamer $gstapiversion element $1 >= $2)
    if [ $gstinspect --exists --atleast-version=$2 $1 ]; then
      AC_MSG_RESULT([found])
      $3
    else
      if [ $gstinspect --exists $1 ]; then
        AC_MSG_RESULT([found, but too old])
      else
        AC_MSG_RESULT([not found])
      fi
      $4
    fi
  fi
])
