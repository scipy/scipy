/*
 * wstdio.c -- $Id$
 * p_stdinit, p_stdout, p_stdin for MS Windows
 *
 * Copyright (c) 1999.  See accompanying LEGAL file for details.
 */

#include "playw.h"
#include "pstdlib.h"

int (*w_stdinit)(void(**)(char*,long), void(**)(char*,long))= 0;
int w_no_mdi = 0;

static void (*w_on_stdin)(char *input_line)= 0;
static void w_formout(char *line, void (*func)(char *, long));

static void (*w_stdout)(char *line, long len)= 0;
static void (*w_stderr)(char *line, long len)= 0;

int console_mode = 0;

void
p_stdinit(void (*on_stdin)(char *input_line))
{
  char *mdi = p_getenv("NO_MDI");
  if (mdi && mdi[0] && (mdi[0]!='0' || mdi[1])) w_no_mdi |= 1;
  if (!w_no_mdi && !w_stdinit) w_no_mdi = 1, AllocConsole();
  if (!w_no_mdi || con_stdinit(&w_stdout, &w_stderr)) {
    w_stdinit(&w_stdout, &w_stderr);
  } else if (w_main_window) {
    /* without this, actual first window created does not show properly */
    ShowWindow(w_main_window, SW_SHOWNORMAL);
    ShowWindow(w_main_window, SW_HIDE);
  }
  w_on_stdin = on_stdin;
}

void p_stdout(char *output_line)
{
  w_formout(output_line, w_stdout);
}

void p_stderr(char *output_line)
{
  w_formout(output_line, w_stderr);
}

char *w_sendbuf(long len)
{
  /* console app: called by worker to get buffer to hold stdin
   * gui app: called by boss to get buffer to hold input line
   *          called by worker w_sendbuf(-1) to retrieve boss buffer
   * therefore must use raw Windows memory manager routines */
  static char *buf = 0;
  if (len >= 0) {
    HANDLE heap = GetProcessHeap();
    if (len <= 0x100) len = 0x100;
    else len = ((len-1)&0x3ff) + 0x100;
    buf = buf? HeapReAlloc(heap, HEAP_GENERATE_EXCEPTIONS, buf, len+1) :
              HeapAlloc(heap, HEAP_GENERATE_EXCEPTIONS, len+1);
  }
  return buf;
}

void
w_deliver(char *buf)
{
  int cr, i, j;
  for (i=j=0 ; buf[j] ; j=i) {
    for (; buf[i] ; i++) if (buf[i]=='\r' || buf[i]=='\n') break;
    cr = (buf[i]=='\r');
    buf[i++] = '\0';
    if (cr && buf[i]=='\n') i++;
    /* deliver one line at a time, not including newline */
    if (w_on_stdin) w_on_stdin(buf+j);
  }
}

static void
w_formout(char *line, void (*func)(char *, long))
{
  if (!p_signalling && line) {
    static char *buf = 0;
    HANDLE heap = GetProcessHeap();
    long len = 256, j = 0, i = 0;
    char c = line[i];
    if (!buf)
      buf = HeapAlloc(heap, HEAP_GENERATE_EXCEPTIONS, 258);
    do {
      if (j >= len) {
        buf = HeapReAlloc(heap, HEAP_GENERATE_EXCEPTIONS, buf, 2*len+2);
        len *= 2;
      }
      if (line[i]!='\n' || c=='\r' || console_mode) {
        if (line[i] == '\a') {
          i++;
          p_feep(0);
          continue;
        } else if (console_mode && line[i]=='\r' && line[i+1]=='\n') {
          i++;
          continue;
        }
      } else {
        buf[j++] = '\r';
      }
      buf[j++] = c = line[i++];
    } while (c);
    func(buf, j-1);
    if (len > 256)
      buf = HeapReAlloc(heap, HEAP_GENERATE_EXCEPTIONS, buf, 258);
  }
}
