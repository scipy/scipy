#! /bin/sh

# Global variables...
AUTOCONF="autoconf"
AUTOHEADER="autoheader"
AUTOM4TE="autom4te"
AUTOMAKE="automake"
ACLOCAL="aclocal"


# We don't use variable here for remembering the type ... strings.
# local variables are not that portable, but we fear namespace issues with
# our includer.  The repeated type calls are not that expensive.
checkAutoconf()
{
  if test -x "`$WHICH autoconf-2.5x`" ; then	
    AUTOCONF="`$WHICH autoconf-2.5x`"
  elif test -x "`$WHICH autoconf-2.54`" ; then
    AUTOCONF="`$WHICH autoconf-2.54`"
  elif test -x "`$WHICH autoconf-2.53`" ; then
    AUTOCONF="`$WHICH autoconf-2.53`"
  elif test -x "`$WHICH autoconf-2.53a`" ; then
    AUTOCONF="`$WHICH autoconf-2.53a`"
  elif test -x "`$WHICH autoconf-2.52`" ; then
    AUTOCONF="`$WHICH autoconf-2.52`"
  elif test -x "`$WHICH autoconf2.50`" ; then
    AUTOCONF="`$WHICH autoconf2.50`"
  fi
}

checkAutoheader()
{
  if test -x "`$WHICH autoheader-2.5x`" ; then
    AUTOHEADER="`$WHICH autoheader-2.5x`"
    AUTOM4TE="`$WHICH autom4te-2.5x`"
  elif test -x "`$WHICH autoheader-2.54`" ; then
    AUTOHEADER="`$WHICH autoheader-2.54`"
    AUTOM4TE="`$WHICH autom4te-2.54`"
  elif test -x "`$WHICH autoheader-2.53`" ; then
    AUTOHEADER="`$WHICH autoheader-2.53`"
    AUTOM4TE="`$WHICH autom4te-2.53`"
  elif test -x "`$WHICH autoheader-2.53a`" ; then
    AUTOHEADER="`$WHICH autoheader-2.53a`"
    AUTOM4TE="`$WHICH autom4te-2.53a`"
  elif test -x "`$WHICH autoheader-2.52`" ; then
    AUTOHEADER="`$WHICH autoheader-2.52`"
  elif test -x "`$WHICH autoheader2.50`" ; then
    AUTOHEADER="`$WHICH autoheader2.50`"
  fi
}

checkAutomakeAclocal ()
{
  if test -z "$UNSERMAKE"; then
   if test -x "`$WHICH automake-1.6`" ; then
      AUTOMAKE="`$WHICH automake-1.6`"
      ACLOCAL="`$WHICH aclocal-1.6`"
    elif test -x "`$WHICH automake-1.7`" ; then
      AUTOMAKE="`$WHICH automake-1.7`"
      ACLOCAL="`$WHICH aclocal-1.7`"
    fi
  else
     AUTOMAKE="$UNSERMAKE"
  fi
}

checkWhich ()
{
  WHICH=""
  for i in "type -p" "which" "type" ; do
    T=`$i sh 2> /dev/null`
    test -x "$T" && WHICH="$i" && break
  done
}

checkWhich
checkAutoconf
checkAutoheader
checkAutomakeAclocal

export WHICH AUTOHEADER AUTOCONF AUTOM4TE AUTOMAKE ACLOCAL
