/*
 * pwin.c -- $Id$
 * routines to create graphics devices for MS Windows
 * 
 * CHANGES:
 * 12/06/04 mdh Correct bug in p_destroy: 
 *              The call to DestroyWindow effectively calls the w_winproc 
 *              window routine with the WM_DESTROY message. It does not simply 
 *              post this message, but calls w_winproc immediately. In the 
 *              w_winproc function, pw is free'd:  p_free(pw);
 *              (line 400 in src/play/win/pscr.c). So, by the time we get to 
 *              the end of the p_destroy routine, pw has been freed and we 
 *              cannot access pw->pbflag, pw->pen, and pw->brush. Moving the 
 *              last two lines to the top of the p_destroy routine solves the 
 *              problem.
 *
 * Copyright (c) 2000.  See accompanying LEGAL file for details.
 */

#include "playw.h"
#include "pstdlib.h"

HWND (*w_parent)(int width, int height, char *title, int hints)= 0;

static p_win *w_pwin(void *ctx, p_scr *s, HDC dc, HBITMAP bm,
                     unsigned long bg);

void
p_cursor(p_win *w, int cursor)
{
  if (cursor==P_NONE) {
    w->cursor = 0;
    SetCursor(0);
    return;
  }
  if (cursor<0 || cursor>P_NONE) cursor = P_SELECT;
  w->cursor = w->s->cursors[cursor];
  if (!w->cursor)
    w->cursor = w->s->cursors[cursor] = w_cursor(cursor);
  SetCursor(w->cursor);
}

p_win *
p_window(p_scr *s, int width, int height, char *title,
         unsigned long bg, int hints, void *ctx)
{
  p_win *pw = w_pwin(ctx, s, 0, 0, bg);
  HWND hw;
  HWND parent = 0;
  DWORD xstyle = 0;
  DWORD style = 0;

  if (!pw) return 0;

  if (w_parent && !(hints & P_DIALOG))
    parent = w_parent(width, height, title, hints);
  pw->ancestor = parent;
  if (parent) {
    style |= WS_CHILD;
  } else {
    RECT rect;
    style |= WS_OVERLAPPEDWINDOW;
    if (hints & P_NORESIZE)
      style &= ~(WS_THICKFRAME | WS_MAXIMIZEBOX);
    rect.left = 0;
    rect.top = 0;
    rect.right = width;
    rect.bottom = height;
    if (AdjustWindowRectEx(&rect, style, 0, xstyle)) {
      width = rect.right - rect.left;
      height = rect.bottom - rect.top;
    }
  }

  hw = CreateWindowEx(xstyle, w_win_class, parent? 0 : title, style,
                     CW_USEDEFAULT, 0, width, height, parent, 0,
                     w_app_instance, pw);
  /* CreateWindow already causes several calls to w_winproc
   * (including WM_CREATE) which will not have GWL_USERDATA set
   * -- WM_CREATE handler fills in and initializes pw */
  if (hw) {
    if (hints & P_RGBMODEL) {
      p_palette(pw, p_595, 225);
      pw->rgb_mode = 1;
    }
    if (!(hints&P_NOKEY)) SetActiveWindow(hw);
    ShowWindow(hw, (hints&P_NOKEY)? SW_SHOWNA : SW_SHOW);
  } else {
    p_destroy(pw);
    pw = 0;
  }

  return pw;
}

static p_win *
w_pwin(void *ctx, p_scr *s, HDC dc, HBITMAP bm, unsigned long bg)
{
  p_win *pw = p_malloc(sizeof(p_win));
  if (pw) {
    pw->ctx = ctx;
    pw->s = s;
    pw->dc = dc;
    pw->w = 0;
    pw->bm = bm;
    pw->menu = 0;

    pw->parent = 0;
    pw->cursor = 0;
    pw->palette = 0;

    if (dc || bm) {
      pw->pixels = 0;
    } else {
      COLORREF *pixels = p_malloc(sizeof(COLORREF)*256);
      if (pixels) {
        int i;
        for (i=0 ; i<=P_EXTRA ; i++) pixels[i] = w_screen.sys_colors[255-P_FG];
        for (; i<256 ; i++) pixels[i] = w_screen.sys_colors[255-i];
        pw->pixels = pixels;
      } else {
        p_free(pw);
        return 0;
      }
      pw->cursor = s->cursors[P_SELECT];
    }
    pw->n_pixels = 0;
    pw->rgb_mode = 0;

    pw->color = P_FG;

    pw->pbflag = 0;
    pw->brush = 0;
    pw->brush_color = P_BG;
    pw->pen = 0;
    pw->pen_color = P_BG;
    pw->pen_width = 1;
    pw->pen_type = P_SOLID;

    pw->font_color = P_BG;
    pw->font = pw->pixsize = pw->orient = 0;

    pw->bg = bg;

    pw->keydown = 0;
    pw->ancestor = 0;

    w_getdc(pw, 24);
  }
  return pw;
}

void
p_destroy(p_win *pw)
{
  if (pw->s && pw->s->font_win==pw) pw->s->font_win = 0;
  if ((pw->pbflag&1) && pw->pen)
    DeleteObject(pw->pen), pw->pen = 0;
  if ((pw->pbflag&2) && pw->brush)
    DeleteObject(pw->brush), pw->brush = 0;
  if (pw->w) {
    pw->ctx = 0;           /* do not call on_destroy handler */
    DestroyWindow(pw->w);  /* WM_DESTROY handler does most of work */
  } else if (pw->bm) {
    /* which order is correct here? */
    DeleteDC(pw->dc), pw->dc = 0;
    DeleteObject(pw->bm), pw->bm = 0;
    pw->pixels = 0;  /* really belongs to parent window */
  } else {
    HENHMETAFILE meta = CloseEnhMetaFile(pw->dc);
    if (meta) DeleteEnhMetaFile(meta);
    pw->dc = 0;
    pw->pixels = 0;  /* really belongs to parent window */
  }
}

p_win *
p_offscreen(p_win *parent, int width, int height)
{
  p_win *pw = 0;
  HDC dc = w_getdc(parent, 0);
  dc = (parent->parent || !dc)? 0 : CreateCompatibleDC(dc);
  if (dc) {
    HBITMAP bm = CreateCompatibleBitmap(parent->dc, width, height);
    if (bm) {
      pw = w_pwin(0, parent->s, dc, bm, parent->bg);
      if (pw) {
        pw->parent = parent;
        SetBkColor(dc, w_color(pw, pw->bg));
        SelectObject(dc, bm);
      } else {
        DeleteObject(bm);
        bm = 0;
      }
    }
    if (!bm) DeleteDC(dc);
  }
  return pw;
}

p_win *
p_menu(p_scr *s, int width, int height, int x, int y,
       unsigned long bg, void *ctx)
{
  p_win *pw = w_pwin(ctx, s, 0, 0, bg);
  HWND hw;
  HWND parent = 0;
  DWORD xstyle = WS_EX_TOPMOST;
  DWORD style = WS_POPUP | WS_VISIBLE;

  if (!pw) return 0;
  pw->menu = 1;

  hw = CreateWindowEx(xstyle, w_menu_class, 0, style,
                     x+s->x0, y+s->y0, width, height, w_main_window, 0,
                     w_app_instance, pw);
  /* see p_window comment */
  if (hw) {
    SetActiveWindow(hw);
    ShowWindow(hw, SW_SHOW);
  } else {
    p_destroy(pw);
    pw = 0;
  }

  return pw;
}

p_win *
p_metafile(p_win *parent, char *filename,
           int x0, int y0, int width, int height, int hints)
{
  /* filename should end in ".emf" or ".EMF"
   * add missing extension here if necessary
   * -- note that this feature can be overridden by terminating
   *    filename with "." (Windows ignores trailing . in filenames) */
  p_win *pw = 0;
  HDC mf = 0;
  HDC dc = parent->parent? 0 : w_getdc(parent, 0);

  if (dc) {
    long dxmm = GetDeviceCaps(dc, HORZSIZE);
    long dymm = GetDeviceCaps(dc, VERTSIZE);
    long dx = GetDeviceCaps(dc, HORZRES);
    long dy = GetDeviceCaps(dc, VERTRES);
    RECT r;  /* in 0.01 mm units */
    int i, n;
    r.left = r.top = 0;
    r.right = (width * dxmm * 100) / dx;
    r.bottom = (height * dymm * 100) / dy;

#ifdef CYGWIN
    filename = u_pathname(filename);  /* result always in p_wkspc.c */
#else
    filename = w_pathname(filename);  /* result always in p_wkspc.c */
#endif
    for (n=0 ; filename[n] ; n++);
    for (i=1 ; i<=4 && i<n ; i++) if (filename[n-i]=='.') break;
    if (i>4) filename[n] = '.', filename[n+1] = 'e', filename[n+2] = 'm',
               filename[n+3] = 'f', filename[n+4] = '\0';
    mf = CreateEnhMetaFile(dc, filename, &r, 0);
    if (mf) {
      pw = w_pwin(0, parent->s, mf, 0, parent->bg);
      if (pw) {
        pw->parent = parent;
        pw->rgb_mode = 1;  /* inhibit p_palette */
        pw->pixels = parent->pixels;
        pw->n_pixels = parent->n_pixels;
      } else {
        DeleteEnhMetaFile((HENHMETAFILE)mf);
      }
    }
  }

  return pw;
}
