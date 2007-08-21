/*
 * pmain.c -- $Id$
 * play main function for UNIX with raw X11 (no X toolkit)
 *
 * Copyright (c) 1998.  See accompanying LEGAL file for details.
 */

#include "play.h"

BEGIN_EXTERN_C
extern void p_mminit(void);
extern int u_main_loop(int (*on_launch)(int,char**), int, char **);
END_EXTERN_C

int
main(int argc, char *argv[])
{
  p_mminit();
  return u_main_loop(on_launch, argc, argv);
}
