/*
 * pathnm.c -- $Id$
 * p_getenv and w_pathname pathname preprocessing for MS Windows
 *
 * Copyright (c) 1999.  See accompanying LEGAL file for details.
 */

#include "playw.h"
#include "pstdlib.h"

#include <string.h>

/* GetEnvironmentVariable is in kernel32.lib */

char *
p_getenv(const char *name)
{
  DWORD flag = GetEnvironmentVariable(name, p_wkspc.c, P_WKSIZ);
  return (flag>0 && flag<P_WKSIZ)? p_wkspc.c : 0;
}

char *
w_pathname(const char *name)
{
  const char *tmp;
  long len = 0;
  unsigned long left = P_WKSIZ;

  if (name[0]=='$') {
    int delim = *(++name);
    if (delim=='(')      { delim = ')';  name++; }
    else if (delim=='{') { delim = '}';  name++; }
    else                   delim = '\0';
    for (tmp=name ; tmp[0] ; tmp++)
      if ((delim && tmp[0]==delim) ||
          tmp[0]=='/' || tmp[0]=='\\') break;
    if (tmp>name+1024) {
      p_wkspc.c[0] = '\0';
      return p_wkspc.c;
    }
    if (tmp>name) {
      char env_name[1024];
      len = tmp-name;
      env_name[0] = '\0';
      strncat(env_name, name, len);
      len = GetEnvironmentVariable(env_name, p_wkspc.c, P_WKSIZ);
      if (len>P_WKSIZ) len = left+1;
      left -= len;
      if (left<0) {
        p_wkspc.c[0] = '\0';
        return p_wkspc.c;
      }
      if (!delim) tmp++;
      name = tmp;
    }
  }

  if (strlen(name)<=left) strcpy(p_wkspc.c+len, name);
  else p_wkspc.c[0] = '\0';

  for (left=0 ; p_wkspc.c[left] ; left++)
    if (p_wkspc.c[left]=='/') p_wkspc.c[left] = '\\';
  return p_wkspc.c;
}

char *
w_unixpath(const char *wname)
{
  int i;
  for (i=0 ; i<P_WKSIZ && wname[i] ; i++)
    p_wkspc.c[i] = (wname[i]=='\\')? '/' : wname[i];
  return p_wkspc.c;
}
