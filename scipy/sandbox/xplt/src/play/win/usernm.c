/*
 * usernm.c -- $Id$
 * p_getuser for MS Windows
 *
 * Copyright (c) 1999.  See accompanying LEGAL file for details.
 */

#include "pstdlib.h"
#include "play.h"

/* GetUserName in advapi32.lib */
#include <windows.h>

char *
p_getuser(void)
{
  DWORD sz = P_WKSIZ;
  if (GetUserName(p_wkspc.c, &sz))
    return p_wkspc.c;
  else
    return 0;
}
