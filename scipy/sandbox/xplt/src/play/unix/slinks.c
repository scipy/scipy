/*
 * slinks.c -- $Id$
 * u_track_link and u_find_exe functions from playu.h
 *
 * Copyright (c) 2000.  See accompanying LEGAL file for details.
 */

/* need these to pick up readlink prototype (unistd.h) */
#ifndef _XOPEN_SOURCE
#define _XOPEN_SOURCE 1
#endif
#ifndef _XOPEN_SOURCE_EXTENDED
#define _XOPEN_SOURCE_EXTENDED 1
#endif

#include "config.h"

#include "playu.h"
#include "pstdlib.h"
#include "pstdio.h"

#include <unistd.h>

char *
u_track_link(const char *name)
{
  if (name) {
    char link[P_WKSIZ+1];
    int i, len;
    if (name != p_wkspc.c)
      for (i=0 ; (p_wkspc.c[i] = name[i]) && i<P_WKSIZ ; i++);
    for (;;) {
      len = readlink(p_wkspc.c, link, P_WKSIZ);
      if (len < 0) break;
      for (i=0 ; i<len ; i++) p_wkspc.c[i] = link[i];
      p_wkspc.c[i] = '\0';
    }
    return p_wkspc.c;
  } else {
    return 0;
  }
}

char *
u_find_exe(const char *argv0)
{
  char *wkspc = p_wkspc.c;
  int i = 0;
  if (!argv0) return 0;

  while (argv0[i] && argv0[i]!='/') i++;

  if (!argv0[i]) {   /* search for argv0 on PATH environment variable */
    char *path = getenv("PATH");
    char c = path? path[0] : 0;
    int s, j, k=0;
    while (c) {
      while (c && c!=':') c = path[k++];
      if (k > 1) {
        for (j=0 ; j<k-1 && j<P_WKSIZ ; j++) wkspc[j] = path[j];
        if (wkspc[j-1] == '/') s = 0;
        else s = 1, wkspc[j] = '/';
        for (; j<k+i && j<P_WKSIZ ; j++) wkspc[j+s] = argv0[j-k+1];
        if (access(wkspc, X_OK) >= 0) break;
      }
      path += k;
      k = 0;
      c = path[0];
    }
    return k? wkspc : 0;
  }

  if (i) {           /* argv0 refers to subdirectory of cwd */
    wkspc = p_getcwd();
    if (wkspc) {
      int j;
      for (j=0 ; wkspc[j] ; j++);
      if (j && wkspc[j-1] != '/') wkspc[j++] = '/';
      while (argv0[0]=='.' && argv0[1]=='/') argv0 += 2;
      for (i=j ; argv0[i-j] && i<P_WKSIZ ; i++) wkspc[i] = argv0[i-j];
    } else {
      i = 0;
      wkspc = p_wkspc.c;
    }

  } else {           /* argv0 is absolute pathname */
    for (i=0 ; argv0[i] && i<P_WKSIZ ; i++) wkspc[i] = argv0[i];
  }

  wkspc[i] = '\0';
  return (access(wkspc, X_OK) >= 0)? wkspc : 0;
}
