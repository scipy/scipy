/*
 * conterm.c -- $Id$
 * console specific part of wstdio.c
 *
 * Copyright (c) 1999.  See accompanying LEGAL file for details.
 */

#include "playw.h"

static void con_sender(void *context);
static void con_stdout(char *output_line, long len);
static void con_stderr(char *output_line, long len);

static HANDLE w_hin = INVALID_HANDLE_VALUE;
static HANDLE w_hout = INVALID_HANDLE_VALUE;
static HANDLE w_herr = INVALID_HANDLE_VALUE;
static HANDLE stdin_ready = 0;
static HANDLE stdin_accepted = 0;
static DWORD WINAPI stdin_main(LPVOID lp);

int
con_stdinit(void (**wout)(char*,long), void (**werr)(char*,long))
{
  HANDLE h;
  *wout = con_stdout;
  *werr = con_stderr;
  /* documentation for MsgWaitForMultipleObjects says:
   * Windows 95: No handle may be a duplicate of another handle
   * created using DuplicateHandle.
   * - except for true console input, files are not on the list
   *   of things any of the wait functions can wait for
   * - apparently, emacs uses the recommended procedure for starting
   *   a subprocess with redirected stdin, which involves DuplicateHandle
   * - in any event, just using w_hin to w_add_input fails in
   *   Win95 for use with emacs shell mode */
  w_hin = GetStdHandle(STD_INPUT_HANDLE);
  w_hout = GetStdHandle(STD_OUTPUT_HANDLE);
  if (w_hin==INVALID_HANDLE_VALUE || w_hout==INVALID_HANDLE_VALUE)
    return 1;
  w_herr = GetStdHandle(STD_ERROR_HANDLE);
  if (w_herr == INVALID_HANDLE_VALUE) w_herr = w_hout;
  stdin_ready = CreateEvent(0,0,0,0);
  stdin_accepted = CreateEvent(0,0,0,0);
  if (stdin_ready && stdin_accepted &&
      w_hin!=INVALID_HANDLE_VALUE && w_hout!=INVALID_HANDLE_VALUE) {
    UINT id;
    h = CreateThread(0,0, stdin_main, 0, 0, &id);
    if (h) {
      w_add_input(stdin_ready, con_sender, 0);
      Sleep(0);
      CloseHandle(h);
    } else {
      CloseHandle(stdin_ready);
      CloseHandle(stdin_accepted);
    }
  }
  return 0;
}

/* separate thread blocks waiting for stdin, signals on arrival */
static DWORD WINAPI
stdin_main(LPVOID lp)
{
  int i = 0;
  DWORD n = 0;
  DWORD sz = 256;
  char *buf = w_sendbuf(sz);
  char tmp[256];
  int ii, no_wait = 1;
  for (;; i=0,sz=256) {
    for (;;) {
      if (!ReadFile(w_hin, tmp, 256, &n, 0))
        return 1;  /* serious error -- stdin unreadable */
      if (tmp[0] == '\03') {
        /* stdin_ready was set -- reset may avoid calling con_sender */
        if (!no_wait) ResetEvent(stdin_ready);
        w_sigint(1);
        /* be sure stdin_accepted is set whether or not con_sender called */
        if (!no_wait) SetEvent(stdin_accepted);
        i = 0;
        break;
      }
      if (no_wait || WaitForSingleObject(stdin_accepted, 0)==WAIT_OBJECT_0) {
        no_wait = 1;
        for (ii=0 ; (unsigned int)ii<n ; ii++) buf[i+ii] = tmp[ii];
        i += ii;
        if (!i || buf[i-1]=='\n') break;
      } else {
        MessageBeep(MB_OK);
      }
      if ((unsigned int)i+256 >= sz) {
        sz += sz;
        buf = w_sendbuf(sz);
      }
    }
    buf[i] = '\0';
    if (i) {
      SetEvent(stdin_ready);
      no_wait = 0;
    }
    if (i>=256) w_sendbuf(256);
  }
  return 0;
}

/* ARGSUSED */
static void
con_sender(void *context)
{
  SetEvent(stdin_accepted);
  w_deliver(w_sendbuf(-1));
}

void con_stdout(char *line, long len)
{
  DWORD n;
  WriteFile(w_hout, line, len, &n, 0);
  FlushFileBuffers(w_hout);
}

void con_stderr(char *line, long len)
{
  DWORD n;
  WriteFile(w_herr, line, len, &n, 0);
  FlushFileBuffers(w_hout);
}
