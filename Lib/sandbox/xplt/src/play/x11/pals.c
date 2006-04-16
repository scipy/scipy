/*
 * pals.c -- $Id$
 * palette handling for X11
 *
 * Copyright (c) 1998.  See accompanying LEGAL file for details.
 */

#include "config.h"
#include "playx.h"
#include "pstdlib.h"

static int x_use_shared(p_scr *s, p_col_t color, p_col_t *pixel);
static void x_lose_shared(p_scr *s, p_col_t *pixels, int n);
static void x_all_shared(p_scr *s, XColor *map /*[256]*/);
static p_col_t x_best_shared(p_scr *s, p_col_t color,
                             XColor *map, int n);
static void x_cavailable(void *p, p_hashkey key, void *ctx);
static void x_list_dead(void *p, p_hashkey key, void *ctx);
static void x_find_dead(void *p, p_hashkey key, void *ctx);
static void x_mark_shared(void *p);

void
p_palette(p_win *w, p_col_t *colors, int n)
{
  p_scr *s = w->s;
  p_col_t r, g, b;
  int i;
  if (n>240) n = 240;

  if (w->parent) w = w->parent;  /* pixmap uses parent palette */
  s->gc_color = -1;

  if (s->vclass==TrueColor || s->vclass==DirectColor) {
    for (i=0 ; i<n ; i++) {
      r = P_R(colors[i]);
      g = P_G(colors[i]);
      b = P_B(colors[i]);
      w->pixels[i] = (s->pixels[r]&s->rmask) |
        (s->pixels[g]&s->gmask) | (s->pixels[b]&s->bmask);
    }
    w->n_palette = n;

  } else if (s->vclass!=PseudoColor) {
    for (i=0 ; i<n ; i++) {
      r = P_R(colors[i]);
      g = P_G(colors[i]);
      b = P_B(colors[i]);
      w->pixels[i] = s->pixels[(r+g+b)/3];
    }
    w->n_palette = n;

  } else if (w->rgb_pixels) {
    for (i=0 ; i<n ; i++) {
      r = P_R(colors[i]);
      g = P_G(colors[i]);
      b = P_B(colors[i]);
      r = (r+32)>>6;
      g = (g+16)>>5;
      b = (b+32)>>6;
      g += b+(b<<3);
      w->pixels[i] = w->rgb_pixels[r+g+(g<<2)];  /* r + 5*g * 45*b */
    }
    w->n_palette = n;

  } else {
    XColor c;
    int p;
    p_col_t fg_pixel = s->colors[1].pixel;
    Display *dpy = s->xdpy->dpy;
    Visual *visual = DefaultVisual(dpy,s->scr_num);
    int map_size = visual->map_entries;
    if (map_size>256) map_size = 256;

    if (w->cmap==None) {
      /* this window uses default colormap (read-only shared colors) */
      int n_old = w->n_palette;

      /* (1) free our use of existing colors */
      w->n_palette = 0;
      x_lose_shared(s, w->pixels, n_old);
      for (i=0 ; i<n_old ; i++) w->pixels[i] = fg_pixel;
      if (n<=0) return;

      /* (2) go for it, but prudently ask in bit reversed order */
      for (i=0 ; i<256 ; i++) {
        p = p_bit_rev[i];
        if (p>=n) continue;
        if (!x_use_shared(s, colors[p], &w->pixels[p])) break;
        w->n_palette++;
      }

      /* (3) fallback if didn't get full request */
      if (w->n_palette<n) {
        XColor map[256];
        int nsh;

        /* query all colors, reserving a use of all sharable colors */
        x_all_shared(s, map);
        for (nsh=0 ; nsh<256 ; nsh++) {
          if (map[nsh].flags) break;
          map[nsh].red =   (map[nsh].red  >>8)&0xff;
          map[nsh].green = (map[nsh].green>>8)&0xff;
          map[nsh].blue =  (map[nsh].blue >>8)&0xff;
        }

        /* take closest sharable pixel to those requested
         * - allocate extra uses for colors used more than once */
        for (; i<256 ; i++) {
          p = p_bit_rev[i];
          if (p>=n) continue;
          w->pixels[p] = x_best_shared(s, colors[p], map, nsh);
          w->n_palette++;
        }

        /* release shared colors we didn't ever use */
        x_lose_shared(s, (p_col_t *)0, 0);
      }

    } else {
      /* this window has a private colormap (read-write private colors) */
      char used[256];

      /* private colormaps will flash display when they are installed
       * take two steps to minimize the annoyance:
       * (1) be sure to allocate the 16 standard colors
       *     s->colors[0:15] the same as default cmap on screen
       *     -- x_color algorithm fails without this,
       *        done when cmap created
       * (2) allocate palette colors from top down, since most
       *     X servers seem to pass out colors from bottom up and
       *     important permanent apps (window manager) started first
       *     -- allocate colors at bottom same as in default cmap */
      for (i=0 ; i<map_size ; i++) used[i] = 0;
      for (i=0 ; i<16 ; i++)
        if (s->colors[i].pixel<map_size) used[s->colors[i].pixel] = 1;
      for (p=map_size-1,i=0 ; p>=0 && i<n ; p--) {
        if (used[p]) continue;
        c.pixel = w->pixels[i] = p;
        c.red = P_R(colors[i]) << 8;
        c.green = P_G(colors[i]) << 8;
        c.blue = P_B(colors[i]) << 8;
        i++;
        c.flags = DoRed | DoGreen | DoBlue;
        XStoreColor(dpy, w->cmap, &c);
      }
      for (; i<240 ; i++) w->pixels[i] = fg_pixel;
      if (p>=0) {
        /* restore as much of the default colormap as possible */
        Colormap cmap = DefaultColormap(dpy,s->scr_num);
        XColor map[256];
        for (i=0 ; i<=p ; i++) map[i].pixel = i;
        XQueryColors(dpy, cmap, map, p+1);
        for (i=0 ; i<=p ; i++) {
          if (used[i]) continue;
          map[i].flags = DoRed | DoGreen | DoBlue;
          XStoreColor(dpy, w->cmap, &map[i]);
        }
      }
      w->n_palette = n;
    }
  }
  if (p_signalling) p_abort();
}

int
x_rgb_palette(p_win *w)
{
  if (w->parent) w = w->parent;
  if (!w->rgb_pixels) {
    p_scr *s = w->s;
    p_col_t *pixels;
    int i;
    if (s->vclass!=PseudoColor) return 0;

    /* would be friendlier to translate preexisting
     * palette, but nuking it is far simpler */
    p_palette(w, p_595, 225);
    x_tmpzap(&s->tmp);
    pixels = s->tmp = p_malloc(sizeof(p_col_t)*256);
    if (!pixels) return 0;
    for (i=0 ; i<256 ; i++) pixels[i] = w->pixels[i];
    s->tmp = 0;
    w->rgb_pixels = pixels;
    p_palette(w, (p_col_t *)0, 0);
  }
  return 1;
}

/*------------------------------------------------------------------------*/

struct x_cshared {
  unsigned long *usepxl;   /* [uses,pixel] for each pixel */
  unsigned long nextpxl;   /* index into next free usepxl pair */
  p_hashtab *bypixel;      /* returns usepxl index given pixel value */
  p_hashtab *bycolor;      /* returns usepxl index given color value */
  /* note: bypixel hash table is unnecessary if we were guaranteed
   *    that the pixel values are <256 for PseudoColor displays,
   *    as I assume in a few places (e.g.- x_all_shared below)
   *    - nevertheless, it ensures that x_use_shared, x_lose_shared
   *      and x_nuke_shared work for all possible X servers
   */
};

struct x_deadpix {
  unsigned long list[256], keys[256];
  int n, m;
  x_cshared *shared;
};
struct x_c2pix {
  XColor *map;
  int n;
};

static int
x_use_shared(p_scr *s, p_col_t color, p_col_t *pixel)
{
  x_cshared *shared = s->shared;
  p_hashkey colkey = P_IHASH(color);
  unsigned long *usepxl;

  usepxl = p_hfind(shared->bycolor, colkey);

  if (!usepxl) {
    Display *dpy = s->xdpy->dpy;
    Colormap cmap = DefaultColormap(dpy, s->scr_num);
    p_hashkey pixkey;
    XColor c;
    unsigned long nextpxl = shared->nextpxl;
    if (nextpxl>=512) return 0;      /* only provide for 256 shared colors */
    c.red = P_R(color) << 8;
    c.green = P_G(color) << 8;
    c.blue = P_B(color) << 8;
    if (!XAllocColor(dpy, cmap, &c)) return 0; /* default colormap is full */
    pixkey = P_IHASH(c.pixel);
    usepxl = p_hfind(shared->bypixel, pixkey);
    if (usepxl) {                  /* different colors may give same pixel */
      XFreeColors(dpy, cmap, &c.pixel, 1, 0UL);
      p_hinsert(shared->bycolor, colkey, usepxl);
    } else {
      usepxl = shared->usepxl + nextpxl;
      shared->nextpxl = usepxl[0];
      usepxl[0] = 0;     /* this will be first use */
      usepxl[1] = c.pixel;
      p_hinsert(shared->bypixel, pixkey, usepxl);
      p_hinsert(shared->bycolor, colkey, usepxl);
    }
  }

  usepxl[0]++;
  *pixel = usepxl[1];
  return 1;
}

static void
x_lose_shared(p_scr *s, p_col_t *pixels, int n)
{
  x_cshared *shared = s->shared;
  unsigned long *usepxl;
  struct x_deadpix deadpix;
  int i;

  if (!shared) {
    shared = p_malloc(sizeof(x_cshared));
    if (!shared) return;
    shared->bycolor = p_halloc(256);
    shared->bypixel = p_halloc(256);
    shared->usepxl = p_malloc(sizeof(unsigned long)*512);
    if (!shared->bycolor || !shared->usepxl) return;
    shared->nextpxl = 0;
    for (i=0 ; i<512 ; i+=2) shared->usepxl[i] = i+2;
    s->shared = shared;
  }

  while ((n--)>0) {
    usepxl = p_hfind(shared->bypixel, P_IHASH(pixels[n]));
    if (usepxl && usepxl[0]) usepxl[0]--;
  }
  deadpix.shared = shared;
  deadpix.n = deadpix.m = 0;
  p_hiter(shared->bycolor, &x_find_dead, &deadpix);
  for (i=0 ; i<deadpix.m ; i++)
    p_hinsert(shared->bycolor, deadpix.keys[i], (void*)0);
  p_hiter(shared->bypixel, &x_list_dead, &deadpix); /* also unlinks */
  for (i=0 ; i<deadpix.n ; i++)
    p_hinsert(shared->bypixel, P_IHASH(deadpix.list[i]), (void*)0);
  if (deadpix.n) {
    Display *dpy = s->xdpy->dpy;
    XFreeColors(dpy, DefaultColormap(dpy,s->scr_num),
                deadpix.list, deadpix.n, 0UL);
  }
}

static void
x_all_shared(p_scr *s, XColor *map /*[256]*/)
{
  x_cshared *shared = s->shared;
  if (shared) {
    int i, n;
    unsigned long *usepxl;
    struct x_c2pix available;
    Display *dpy = s->xdpy->dpy;
    Colormap cmap = DefaultColormap(dpy,s->scr_num);
    Visual *visual = DefaultVisual(dpy,s->scr_num);
    int map_size = visual->map_entries;
    if (map_size>256) map_size = 256;

    for (i=n=0 ; i<map_size ; i++)
      if (!p_hfind(shared->bypixel, P_IHASH(i)))
        map[n++].pixel = i;
    if (!n) return;       /* we already own all shared colors */

    XQueryColors(dpy, cmap, map, n);
    for (i=0 ; i<n ; i++)
      if (XAllocColor(dpy, cmap, &map[i])) {
        /* this color is sharable, record in bypixel table
         * -includes the standard colors */
        usepxl = shared->usepxl + shared->nextpxl;
        shared->nextpxl = usepxl[0];
        usepxl[0] = 0;     /* we haven't actually used it yet */
        usepxl[1] = map[i].pixel;
        p_hinsert(shared->bypixel, P_IHASH(map[i].pixel), usepxl);
        if (shared->nextpxl>=512) break;
      }
    available.map = map;
    available.n = 0;
    p_hiter(shared->bypixel, &x_cavailable, &available);
    XQueryColors(dpy, cmap, map, available.n);
    for (i=0 ; i<available.n ; i++) map[i].flags = 0;
    if (i<256) map[i].flags = 1;
  }
}

static p_col_t
x_best_shared(p_scr *s, p_col_t color, XColor *map, int n)
{
  unsigned long *usepxl;
  long d1, d0, tmp;
  int j, k;
  int r = P_R(color);
  int g = P_G(color);
  int b = P_B(color);
  d0 = 3*256*256;
  for (k=j=0 ; k<n ; k++) {
    d1 =  ((tmp = map[k].red   - r), tmp*tmp);
    d1 += ((tmp = map[k].green - g), tmp*tmp);
    d1 += ((tmp = map[k].blue  - b), tmp*tmp);
    if (d1<d0) j = k, d0 = d1;
  }
  usepxl = p_hfind(s->shared->bypixel, P_IHASH(map[j].pixel));
  usepxl[0]++;
  return usepxl[1];
}

void
x_nuke_shared(p_scr *s)
{
  x_cshared *shared = s->shared;
  if (shared) {
    unsigned long *usepxl = shared->usepxl;
    void (*noaction)(void *)= 0;
    Display *dpy = s->xdpy->dpy;
    int i, n;
    s->shared = 0;
    p_hfree(shared->bypixel, noaction);
    p_hfree(shared->bycolor, &x_mark_shared);
    for (i=n=0 ; i<512 ; i+=2)
      if (usepxl[i]==1) usepxl[n++] = usepxl[i+1];
    if (n)
      XFreeColors(dpy, DefaultColormap(dpy,s->scr_num), usepxl, n, 0UL);
    p_free(usepxl);
    p_free(shared);
  }
  if (p_signalling) p_abort();
}

/* ARGSUSED */
static void
x_cavailable(void *p, p_hashkey key, void *ctx)
{
  unsigned long *usepxl = p;
  struct x_c2pix *available = ctx;
  available->map[available->n++].pixel = usepxl[1];
}

/* ARGSUSED */
static void
x_list_dead(void *p, p_hashkey key, void *ctx)
{
  unsigned long *usepxl = p;
  if (!usepxl[0]) {
    struct x_deadpix *deadpix = ctx;
    deadpix->list[deadpix->n++] = usepxl[1];
    usepxl[0] = deadpix->shared->nextpxl;
    deadpix->shared->nextpxl = usepxl - deadpix->shared->usepxl;
  }
}

static void
x_find_dead(void *p, p_hashkey key, void *ctx)
{
  unsigned long *usepxl = p;
  if (!usepxl[0]) {
    struct x_deadpix *deadpix = ctx;
    deadpix->keys[deadpix->m++] = key;
  }
}

static void
x_mark_shared(void *p)
{
  unsigned long *usepxl = p;
  usepxl[0] = 1;  /* 1 cannot be in nextpxl free list (always even) */
}
