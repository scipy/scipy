/*
 * cygapp.c -- $Id$
 * cygwin (or uwin?) main program stub
 *
 * Copyright (c) 2000.  See accompanying LEGAL file for details.
 */

#include "playw.h"
#include "pstdlib.h"

extern int cyg_app(int (*on_launch)(int, char **),
                   HINSTANCE me, LPSTR cmd_line, int show0);

static int (*w_on_launch)(int ,char **)= 0;
static void w_get_cmds(LPSTR cmd_line);
static int w_argc = 0;
static char **w_argv = 0;
static char *w_str_cpy(HANDLE heap, const char *text, long len);
static int w_result = 0;

static int w_loop(void);
static BOOL WINAPI w_ctrlc_handler(DWORD type);

static void cyg_quit(void);

int
cyg_app(int (*on_launch)(int, char **),
        HINSTANCE me, LPSTR cmd_line, int show0)
{
  HMODULE hm = me;
  HWND hw = CreateWindow(w_win_class, "dummy main wnd", 0,
                        CW_USEDEFAULT, 0, 1, 1, 0, 0, hm, 0);

  w_initialize(hm, hw, cyg_quit, 0, 0);

  /* this should not be necessary for windows subsystem */
  SetConsoleCtrlHandler(w_ctrlc_handler, 1);

  w_on_launch = on_launch;
  w_get_cmds(cmd_line);
  w_protect(w_loop);
  if (w_result) return w_result;

  return w_on_quit();
}

static void
cyg_quit(void)
{
  PostQuitMessage(0);
}

static int
w_loop(void)
{
  MSG msg;
  if (w_on_launch) {
    int (*on_launch)(int, char **)= w_on_launch;
    w_on_launch = 0;
    p_mminit();
    w_pollinit();
    w_result = on_launch(w_argc, w_argv);
    if (w_result) return w_result;
  }
  /* event loop closely matches mfc default CWinApp::Run */
  for (;;) {
    while (!PeekMessage(&msg, 0, 0,0, PM_REMOVE))
      if (!w_work_idle()) {
        GetMessage(&msg, 0, 0,0);
        break;
      }
    if (msg.message == WM_QUIT) break;
    if (!w_app_msg(&msg)) {
      TranslateMessage(&msg);
      DispatchMessage(&msg);
    }
  }
  return msg.wParam;
}

static void
w_get_cmds(LPSTR cmd_line)
{
  HANDLE heap = GetProcessHeap();
  char module_name[1028];
  DWORD len = GetModuleFileName(w_app_instance, module_name, 1024);
  module_name[len] = '\0';
  w_argc = 0;
  w_argv = (char**)HeapAlloc(heap, HEAP_GENERATE_EXCEPTIONS, sizeof(char *)*9);
  w_argv[w_argc++] = w_str_cpy(heap, w_unixpath(module_name), len);
  if (cmd_line) {
    char *c = cmd_line;
    char delim;
    for (;;) {
      while (c[0]==' ' || c[0]=='\t' || c[0]=='\r' || c[0]=='\n') c++;
      delim = c[0];
      if (!delim) break;
      cmd_line = c;
      if (delim=='"' || delim=='\'') {
        cmd_line = ++c;
        while (c[0] && c[0]!=delim) c++;
      } else {
        while (c[0] && c[0]!=' ' && c[0]!='\t' &&
               c[0]!='\r' && c[0]!='\n') c++;
        delim = 'x';
      }
      if (w_argc>1 || cmd_line[0]!='-'||cmd_line[1]!='n'||cmd_line[2]!='o'||
          cmd_line[3]!='m'||cmd_line[4]!='d'||cmd_line[5]!='i') {
        if (!(w_argc&7))
          w_argv = (char **)HeapReAlloc(heap, HEAP_GENERATE_EXCEPTIONS,
                                       w_argv, sizeof(char *)*(2*w_argc+1));
        w_argv[w_argc++] = w_str_cpy(heap, cmd_line, c - cmd_line);
      } else {
        w_no_mdi = 1;
      }
      if (c[0] == delim) c++;
      cmd_line = c;
    }
  }
  w_argv[w_argc] = 0;
}

static char *
w_str_cpy(HANDLE heap, const char *text, long len)
{
  if (!len) while (text[len]) len++;
  {
    char *buf = (char *)HeapAlloc(heap, HEAP_GENERATE_EXCEPTIONS, len+1);
    char *buffer = buf;
    while (len--) *buf++= *text++;
    buf[0] = '\0';
    return buffer;
  }
}

static BOOL WINAPI
w_ctrlc_handler(DWORD type)
{ /* this function runs in its own thread */
  return (type==CTRL_C_EVENT)? w_sigint(1) : FALSE;
}
