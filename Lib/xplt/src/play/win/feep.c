/*
 * feep.c -- $Id$
 * p_feep for MS Windows
 *
 * Copyright (c) 1999.  See accompanying LEGAL file for details.
 */

#include "playw.h"

/* ARGSUSED */
void
p_feep(p_win *w)
{
  MessageBeep(MB_OK);
}
