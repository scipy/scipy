/*
 * pathnm.c -- $Id$
 * p_getenv and u_pathname pathname preprocessing for UNIX machines
 *
 * Copyright (c) 1998.  See accompanying LEGAL file for details.
 */

#ifndef _POSIX_SOURCE
#define _POSIX_SOURCE 1
#endif

#include "config.h"
#include "play.h"
#include "playu.h"
#include "pstdlib.h"

#include <string.h>

#ifndef NO_PASSWD
#include <sys/types.h>
#include <unistd.h>
#include <pwd.h>
#endif

char *
p_getenv(const char *name)
{
  return getenv(name);
}

char *
u_pathname(const char *name)
{
  const char *tmp;
  long len = 0;
  long left = P_WKSIZ;

  p_wkspc.c[0] = '\0';

  if (name[0]=='$') {
    int delim = *(++name);
    if (delim=='(')      { delim = ')';  name++; }
    else if (delim=='{') { delim = '}';  name++; }
    else                   delim = '/';
    tmp = strchr(name, delim);
    if (!tmp) tmp = name+strlen(name);
    if (tmp>name+left) return p_wkspc.c;
    if (tmp>name) {
      char *env;
      len = tmp-name;
      strncat(p_wkspc.c, name, len);
      env = getenv(p_wkspc.c);
      if (!env) return &p_wkspc.c[len];
      len = strlen(env);
      left -= len;
      if (left<0) return &p_wkspc.c[len];
      strcpy(p_wkspc.c, env);
      if (delim!='/') tmp++;
      name = tmp;
    }

  } else if (name[0]=='~') {
    char *home = 0;
    tmp = strchr(++name, '/');
    if (!tmp) {
      len = strlen(name);
      if (len>left) return p_wkspc.c;
      strcpy(p_wkspc.c, name);
      name += len;
    } else {
      len = tmp-name;
      if (tmp>name+left) return p_wkspc.c;
      if (len) strncat(p_wkspc.c, name, len);
      name = tmp;
    }
#ifndef NO_PASSWD
    if (!p_wkspc.c[0]) {
      home = getenv("HOME");
      if (!home) {
        struct passwd *pw = getpwuid(getuid());  /* see also usernm.c */
        if (pw) home = pw->pw_dir;
      }
    } else {
      struct passwd *pw = getpwnam(p_wkspc.c);
      if (pw) home = pw->pw_dir;
    }
#else
    home = p_wkspc.c[0]? 0 : getenv("HOME");
#endif
    if (!home) return &p_wkspc.c[len];
    len = strlen(home);
    left -= len;
    if (left<0) return &p_wkspc.c[len];
    strcpy(p_wkspc.c, home);
  }

  if (strlen(name)<=left) strcpy(p_wkspc.c+len, name);
  else p_wkspc.c[0] = '\0';
  return p_wkspc.c;
}
