/*
     HOST.C

     $Id$

     This is lifted from the MIT Xlib source, in case GIST must be
     used with a non-MIT Xlib which does not have _XGetHostname.
     If present, _XGetHostname should be in libX11.a.
 */

/*
 * $XConsortium: XlibInt.c,v 11.119 89/12/12 11:50:18 jim Exp $
 */

/*

Copyright 1985, 1986, 1987, 1988, 1989 by the
Massachusetts Institute of Technology

Permission to use, copy, modify, distribute, and sell this software and its
documentation for any purpose is hereby granted without fee, provided that
the above copyright notice appear in all copies and that both that
copyright notice and this permission notice appear in supporting
documentation, and that the name of M.I.T. not be used in advertising or
publicity pertaining to distribution of the software without specific,
written prior permission.  M.I.T. makes no representations about the
suitability of this software for any purpose.  It is provided "as is"
without express or implied warranty.

*/

#include <X11/Xos.h>

/*
 * _XGetHostname - similar to gethostname but allows special processing.
 */
int _XGetHostname (char*, int);

int _XGetHostname (buf, maxlen)
    char *buf;
    int maxlen;
{
    int len;

#ifdef hpux				/* stupid makedepend [need if] */
#define NEED_UTSNAME
#endif
#ifdef USG
#define NEED_UTSNAME
#endif

#ifdef NEED_UTSNAME
#include <sys/utsname.h>
    /*
     * same host name crock as in server and xinit.
     */
    struct utsname name;

    uname (&name);
    len = strlen (name.nodename);
    if (len >= maxlen) len = maxlen - 1;
    strncpy (buf, name.nodename, len);
    buf[len] = '\0';
#else
    buf[0] = '\0';
    (void) gethostname (buf, maxlen);
    buf [maxlen - 1] = '\0';
    len = strlen(buf);
#endif /* hpux */
    return len;
}
