/*
 * cygmain.c -- $Id$
 * cygwin (or uwin?) main program stub
 *
 * Copyright (c) 2000.  See accompanying LEGAL file for details.
 */

#include <windows.h>

int WINAPI
WinMain(HINSTANCE me, HINSTANCE prev, LPSTR cmd_line, int show0)
{
  extern int on_launch(int, char **);
  extern int cyg_app(int (*on_launch)(int, char **),
                     HINSTANCE me, LPSTR cmd_line, int show0);
  return cyg_app(on_launch, me, cmd_line, show0);
}
