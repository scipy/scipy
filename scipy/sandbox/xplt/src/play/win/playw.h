/*
 * playw.h -- $Id$
 * MS Windows-private portability layer declarations
 *
 * Copyright (c) 1999.  See accompanying LEGAL file for details.
 */

#include "play.h"

#ifndef __AFX_H__
#include <windows.h>
#endif

#if defined(__cplusplus) || defined(c_plusplus)
extern "C" {
#endif

extern void w_initialize(HINSTANCE i, HWND w,  void (*wquit)(void),
                         int (*wstdinit)(void(**)(char*,long),
                                         void(**)(char*,long)),
                         HWND (*wparent)(int, int, char *, int));
extern int w_on_quit(void);

extern void w_caught(void);
extern int w_sigint(int delay);
extern void w_siginit(void);
extern void w_fpu_setup(void);
extern int w_protect(int (*run)(void));

extern DWORD w_id_worker;  /* required for PostThreadMessage */
extern HANDLE w_worker;
extern int w_work_idle(void);
extern void w_pollinit(void);
extern UINT w_add_msg(void (*on_msg)(MSG *));
extern int w_app_msg(MSG *msg);
extern int (*w_msg_hook)(MSG *msg);
extern int w_add_input(HANDLE wait_obj, void (*on_input)(void *),
                       void *context);

extern int w_no_mdi;
extern int con_stdinit(void(**)(char*,long), void(**)(char*,long));
extern int (*w_stdinit)(void (**wout)(char*,long),
                        void (**werr)(char*,long));
extern void w_deliver(char *buf);   /* calls on_stdin */
extern char *w_sendbuf(long len);

extern int w_nwins;  /* count of graphics windows */
extern int w_nputs;  /* count of w_add_input input sources */

extern HINSTANCE w_app_instance;
extern HWND w_main_window;
extern HWND (*w_parent)(int width, int height, char *title, int hints);
extern LRESULT CALLBACK w_winproc(HWND hwnd, UINT uMsg,
                                  WPARAM wParam, LPARAM lParam);

extern char *w_pathname(const char *name);
extern char *w_unixpath(const char *wname);

/* if set non-zero, p_abort will call this */
extern void (*w_abort_hook)(void);

/* ------------------------------------------------------------------------ */

#define W_FONTS_CACHED 6

struct p_scr {
  int width, height, depth;
  int x0, y0;                  /* usable area may not start at (0,0) */
  int does_linetypes;          /* Win9x can't do dashed lines */
  int does_rotext;             /* may not be able to rotate text */

  COLORREF sys_colors[15], sys_index[15];
  PALETTEENTRY *sys_pal;
  int sys_offset;

  HFONT gui_font;
  int font_order[W_FONTS_CACHED];  /* indices into font_cache */
  struct {
    HFONT hfont;
    int font, pixsize, orient;
  } font_cache[W_FONTS_CACHED];
  p_win *font_win;             /* most recent window that called p_font */

  HCURSOR cursors[P_NONE];

  p_win *first_menu;
  p_win *active;               /* last p_win which set foreground palette */
};

struct p_win {
  void *ctx;
  p_scr *s;

  HDC dc;
  HWND w;
  HBITMAP bm;
  int menu;

  p_win *parent;
  HCURSOR cursor;
  HPALETTE palette;

  /* DC keeps three separate colors, only set when actually used */
  COLORREF *pixels;
  int n_pixels, rgb_mode;
  unsigned long color;       /* pending color */

  int pbflag;                /* 1 set if null brush installed in dc
                              * 2 set if null pen installed in dc */
  HBRUSH brush;              /* brush to be installed in dc */
  unsigned long brush_color; /* actual color of brush */
  HPEN pen;                  /* pen to be installed in dc */
  unsigned long pen_color;   /* actual color of pen */
  int pen_width, pen_type;   /* pending pen width and type */

  unsigned long font_color;  /* actual color of current font */
  int font, pixsize, orient;

  unsigned long bg;          /* for p_clear */

  unsigned long keydown;     /* required to interpret WM_CHAR */
  HWND ancestor;             /* p_destroy or p_resize communication */
};

extern p_scr w_screen;
extern HCURSOR w_cursor(int cursor);
extern LPCTSTR w_win_class;
extern LPCTSTR w_menu_class;

extern COLORREF w_color(p_win *w, unsigned long color);
extern HDC w_getdc(p_win *w, int flags);

/* all bitmap functions seem to get colors backwards?? */
#define W_SWAPRGB(c) (((c)&0xff00)|(((c)>>16)&0xff)|(((c)&0xff)<<16))

extern POINT w_pt_list[2050];
extern int w_pt_count;

#if defined(__cplusplus) || defined(c_plusplus)
}
#endif
