/*
 * test2d.c -- $Id$
 * test code to exercise play features
 *
 * Copyright (c) 1998.  See accompanying LEGAL file for details.
 */

/* set these to 1/0 during development to turn on/off
 * exercising the corresponding play features */
#define TRY_STDINIT   1
#define TRY_HANDLER   1
#define TRY_GUI       1

#define TRY_PSTDIO    1

/* TRY_GRAPHICS prerequisites:
 *      p_connect, p_window, p_text, p_lines, and p_fill
 * -implement those functions first
 * -some features (e.g.- rotated text) need not be implemented initially */
#define TRY_GRAPHICS  1
#define TRY_RECT      1
#define TRY_ELLIPSE   1
#define TRY_DOTS      1
#define TRY_SEGMENTS  1
#define TRY_NDXCELL   1
#define TRY_RGBCELL   1

#define TRY_PALETTE   1
#define TRY_CLIP      1
#define TRY_CURSOR    1
#define TRY_RESIZE    1

#define TRY_CLIPBOARD 1

#define TRY_OFFSCREEN 1
#define TRY_MENU      1
#define TRY_METAFILE  1

#define TRY_RGBREAD   1

/* #include "config.h" -- really must work without this */

#include "play.h"
#include "pstdlib.h"
#include "pstdio.h"
#include <stdio.h>
#include <string.h>

#if TRY_HANDLER
#include <signal.h>
#ifdef SIGBUS
# define MY_SIGBUS SIGBUS
#else
# define MY_SIGBUS 0
#endif
#ifdef SIGPIPE
# define MY_SIGPIPE SIGPIPE
#else
# define MY_SIGPIPE 0
#endif
static int sig_table[] = {
  0, 0, SIGINT, SIGFPE, SIGSEGV, SIGILL, MY_SIGBUS, MY_SIGPIPE };
static char *sig_name[] = {
  "PSIG_NONE", "PSIG_SOFT", "PSIG_INT", "PSIG_FPE", "PSIG_SEGV",
  "PSIG_ILL", "PSIG_BUS", "PSIG_IO", "PSIG_OTHER" };
#endif

extern int on_idle(void);
extern void on_alarm(void *context);
static void idle_act(char *args);
static void set_act(char *args);
static void clr_act(char *args);

extern int on_quit(void);

#if TRY_STDINIT
extern void on_stdin(char *input_line);
static void help_act(char *args);
#endif

#if TRY_PSTDIO
static void ls_act(char *args);
static void cd_act(char *args);
static void pwd_act(char *args);
static void mv_act(char *args);
static void rm_act(char *args);
static void mkdir_act(char *args);
static void rmdir_act(char *args);
#endif

#if TRY_HANDLER
extern void on_exception(int signal, char *errmsg);
static void raise_act(char *args);
static void abort_act(char *args);
#endif

#if TRY_GRAPHICS
extern int phints;
extern p_scr *scr, *scr2, *scremote;
p_scr *scr, *scr2, *scremote;
extern p_win *win1, *win2, *win3, *win4;
p_win *win1, *win2, *win3, *win4, *offscr;
extern void dis_graphics(void);
extern void simple_create(p_win **winp);
extern void simple_draw(p_scr *scr, p_win *w, int x0, int y0);
extern void advanced_draw(p_win *w);
extern void seg_draw(p_win *w, int x0, int y0, int x1, int y1);
extern void box_draw(p_win *w, int x0, int y0, int x1, int y1, int fillit);
extern void tri_draw(p_win *w, int x0,int y0, int x1,int y1, int x2,int y2,
                     int fillit);
extern void curve_draw(p_win *w, int x0, int y0, int width, int height);
static void feep_act(char *args);
static void redraw_act(char *args);
static void private_act(char *args);
static void rgb_act(char *args);
extern void hilite(p_win *w, int x0, int y0, int x1, int y1, int onoff);
static void make_rainbow(int pn, p_col_t *colors);
# if TRY_RESIZE
static void resize_act(char *args);
# endif
# if TRY_CLIPBOARD
static void paste_act(char *args);
# endif
static int hcour, bcour, wcour;
# if TRY_GUI
extern void on_focus(void *c,int in);
extern void on_key(void *c,int k,int md);
extern void on_click(void *c,int b,int md,int x,int y,
                     unsigned long ms);
extern void on_motion(void *c,int md,int x,int y);
extern void on_expose(void *c, int *xy);
extern void on_deselect(void *c);
extern void on_destroy(void *c);
extern void on_resize(void *c,int w,int h);
static int xfoc, xkey, xclick, xmotion, ygui, xhi0, yhi0, xhi1, yhi1, hion;
extern void on_panic(p_scr *screen);
# endif
# if TRY_METAFILE
static void meta_act(char *args);
# endif
# if TRY_PALETTE
static void pal_act(char *args);
static int pal_number = 1;
static void reset_palette(void);
# endif
# if TRY_RGBCELL
static void rgbcell_act(char *args);
# endif
# if TRY_OFFSCREEN
extern void on_animate(void *context);
# endif
#endif
static int rgbcell_on = 0;

static int app_argc = 0;
static char **app_argv = 0;
static int app_quit = 0;

static int prompt_issued = 0;

int
on_quit(void)
{
#if TRY_GRAPHICS
  dis_graphics();
#endif
#if TRY_STDINIT
  sprintf(p_wkspc.c, "\non_quit called, returning %d\n", app_quit);
  p_stdout(p_wkspc.c);
#endif
  return app_quit;
}

#if TRY_STDINIT
static void
quit_act(char *args)
{
  app_quit = (int)strtol(args, (char **)0, 0);
  p_quit();
}

struct command {
  char *name;
  void (*action)(char *args);
  char *desc;
} commands[] = {
  { "help", &help_act, "help\n" },
  { "quit", &quit_act, "quit [value]\n" },
  { "idle", &idle_act, "idle [ncalls]\n" },
  { "set", &set_act, "set alarm_secs\n" },
  { "clr", &clr_act, "clr [alarm_number]\n" },
#if TRY_HANDLER
  { "raise", &raise_act, "raise [signal_number]\n" },
  { "abort", &abort_act, "abort\n" },
#endif
#if TRY_PSTDIO
  { "ls", &ls_act, "ls [path]\n" },
  { "cd", &cd_act, "cd [path]\n" },
  { "pwd", &pwd_act, "pwd\n" },
  { "mv", &mv_act, "mv old_path new_path\n" },
  { "rm", &rm_act, "rm path\n" },
  { "mkdir", &mkdir_act, "mkdir path\n" },
  { "rmdir", &rmdir_act, "rmdir path\n" },
#endif
#if TRY_GRAPHICS
  { "feep", &feep_act, "feep\n" },
  { "redraw", &redraw_act, "redraw\n" },
  { "private", &private_act, "private\n" },
  { "rgb", &rgb_act, "rgb\n" },
# if TRY_RESIZE
  { "resize", &resize_act, "resize width height\n" },
# endif
# if TRY_CLIPBOARD
  { "paste", &paste_act, "paste\n" },
# endif
# if TRY_METAFILE
  { "meta", &meta_act, "meta filename\n" },
# endif
# if TRY_PALETTE
  { "palette", &pal_act, "palette\n" },
# endif
# if TRY_RGBCELL
  { "rgbcell", &rgbcell_act, "rgbcell\n" },
# endif
#endif
  { (char*)0, (void (*)(char*))0 }
};

/* ARGSUSED */
static void
help_act(char *args)
{
  int i;
  for (i=0 ; commands[i].desc ; i++)
    p_stdout(commands[i].desc);
}

void
on_stdin(char *input_line)
{
  struct command *cmnd;
  int n = strspn(input_line, " \t\n\r");
  input_line += n;
  n = strcspn(input_line, " \t\n\r");
  prompt_issued = 0;
  if (n) {
    for (cmnd=commands ; cmnd->name ; cmnd++) {
      if (strncmp(cmnd->name, input_line, n)) continue;
      if (!cmnd->action) break;
      strcpy(p_wkspc.c, "doing: ");
      strncat(p_wkspc.c, input_line, n);
      strcat(p_wkspc.c, "\n");
      p_stdout(p_wkspc.c);
      input_line += n;
      input_line += strspn(input_line, " \t\n\r");
      cmnd->action(input_line);
      return;
    }

    strcpy(p_wkspc.c, "\ntest2d command not recognized: ");
    strncat(p_wkspc.c, input_line, n);
    strcat(p_wkspc.c, "\n");
    p_stderr(p_wkspc.c);
  }
}
#endif

#if TRY_HANDLER

void
on_exception(int signal, char *errmsg)
{
  p_qclear();

#if TRY_STDINIT
  if (signal<0) signal = 0;
  else if (signal>PSIG_OTHER) signal = PSIG_OTHER;
  sprintf(p_wkspc.c, "test2d received signal %s\n", sig_name[signal]);
  p_stdout(p_wkspc.c);
  if (errmsg) {
    sprintf(p_wkspc.c, "  with errmsg = %s\n", errmsg);
    p_stdout(p_wkspc.c);
  }
  prompt_issued = 0;
#endif
}

static void
raise_act(char *args)
{
  char *after = args;
  int n = args[0]? (int)strtol(args, &after, 0) : 0;
#if TRY_STDINIT
  if (after==args) {
    while (args[-1]=='\n' || args[-1]=='\r' ||
           args[-1]==' ' || args[-1]=='\t') args--;
    n = (args[-1]-'e');  /* this should be non-optimizable n=0 */
    sprintf(p_wkspc.c, "SIGFPE handling broken: 1.0/%d = %g\n",
            n, 1.0/(double)n);
    p_stdout(p_wkspc.c);
    return;
  }
#endif
  if (n<0) n = 0;
  else if (n>=PSIG_OTHER) n = 1001;
  else n = sig_table[n];
  raise(n);
}

static void
abort_act(char *args)
{
  p_abort();
}
#endif

static unsigned int app_ncalls = 0;
static int panic_count = 0;

int
on_launch(int argc, char *argv[])
{
  app_argc = argc;
  app_argv = argv;
  p_quitter(&on_quit);
  p_idler(&on_idle);
#if TRY_STDINIT
  p_stdinit(&on_stdin);
  p_stdout("test2d initialized stdio\n");
#endif
#if TRY_HANDLER
  p_handler(&on_exception);
#endif
#if TRY_GUI
  p_gui(&on_expose, &on_destroy, &on_resize, &on_focus,
        &on_key, &on_click, &on_motion, &on_deselect, &on_panic);
#endif
  return 0;
}

int
on_idle(void)
{
  if (app_ncalls) {
#if TRY_STDINIT
    p_stdout("on_idle called after idle command\n");
#endif
    app_ncalls--;
    prompt_issued = (app_ncalls!=0);
  }
#if TRY_GRAPHICS
  if (!scr && (panic_count<3)) {
    int sw, sh;
    scr = p_connect(0);
    simple_create(&win1);
    win2 = p_window(scr, 400, 400, "test2d advanced", P_BG,
                   phints | P_NORESIZE | P_NOKEY | P_NOMOTION, &win2);
    p_sshape(scr, &sw, &sh);
# if TRY_STDINIT
    sprintf(p_wkspc.c, "test2d: screen shape %d X %d pixels\n", sw, sh);
    p_stdout(p_wkspc.c);
    prompt_issued = 0;
# endif
  }
#endif
#if TRY_STDINIT
  if (!prompt_issued) {
    p_stdout("idle> ");
    prompt_issued = 1;
  }
#endif
  return (app_ncalls!=0);
}

static struct alarm_data {
  double timeout, set_at, rang_at;
  int state;
} alarm_list[5];

void
on_alarm(void *context)
{
  struct alarm_data *alarm = context;
  int n = alarm-alarm_list;
  if (alarm->state==-1) {
    alarm->state = -2;
#if TRY_STDINIT
    sprintf(p_wkspc.c, "test2d: ERROR cancelled alarm #%d rang\n", n);
    p_stderr(p_wkspc.c);
#endif
  } else {
    alarm->rang_at = p_wall_secs();
    alarm->state = 0;
#if TRY_STDINIT
    sprintf(p_wkspc.c,
            "test2d: alarm #%d rang after %g secs (requested %g)\n", n,
            alarm->rang_at-alarm->set_at, alarm->timeout);
    p_stdout(p_wkspc.c);
#endif
  }
  prompt_issued = 0;
}

static void
idle_act(char *args)
{
  app_ncalls = (unsigned int)strtol(args, (char **)0, 0);
  if (app_ncalls<1) app_ncalls = 1;
  if (app_ncalls>5) app_ncalls = 5;
}

static void
set_act(char *args)
{
  double t = strtod(args, (char **)0);
  int n;
  for (n=0 ; n<5 ; n++) if (alarm_list[n].state<=0) break;
  if (n>=5) {
#if TRY_STDINIT
    p_stdout("test2d: (no free alarms, set ignored)\n");
#endif
    return;
  }
  alarm_list[n].timeout = t>0.? t : 0.;
  alarm_list[n].state = 1;
  alarm_list[n].set_at = p_wall_secs();
  alarm_list[n].rang_at = 0.0;
  p_set_alarm(t, &on_alarm, &alarm_list[n]);
#if TRY_STDINIT
  sprintf(p_wkspc.c,
          "test2d: set alarm #%d for %g secs\n", n, alarm_list[n].timeout);
  p_stdout(p_wkspc.c);
#endif
}

static void
clr_act(char *args)
{
  int n = (int)strtol(args, (char **)0, 0);
  if (!args[0])
    for (n=0 ; n<5 ; n++) if (alarm_list[n].state==1) break;
  if (n<0 || n>=5) {
    p_stdout("test2d: (no such alarm, clr ignored)\n");
    return;
  }
  p_clr_alarm((n%2)? &on_alarm : 0, &alarm_list[n]);
  alarm_list[n].state = -1;
#if TRY_STDINIT
  sprintf(p_wkspc.c, "test2d: clr alarm #%d\n", n);
  p_stdout(p_wkspc.c);
#endif
}

#if TRY_PSTDIO
static void
ls_act(char *args)
{
  p_dir *dir;
  char *dot = ".";
  char *path = strtok(args, " \t\r\n");
  if (!path) path = dot;
  dir = p_dopen(path);
  if (dir) {
    p_file *f;
    int is_dir;
    char *name = p_dnext(dir, &is_dir);
    long len = strlen(path);
    char *dirname = p_strncat(path, (len<1 || path[len-1]!='/')? "/":"", 0);
    for ( ; name ; name=p_dnext(dir, &is_dir)) {
      if (!is_dir) {
        strcpy(p_wkspc.c, dirname);
        strncat(p_wkspc.c, name, strlen(name));
        f = p_fopen(p_wkspc.c, "rb");
        if (f) {
          len = p_fsize(f);
          p_fclose(f);
        } else {
          len = -1;
        }
        sprintf(p_wkspc.c, "%s %-20s %ld\n", "FIL", name, len);
      } else {
        sprintf(p_wkspc.c, "%s %s\n", "DIR", name);
      }
      p_stdout(p_wkspc.c);
    }
    p_free(dirname);
    p_dclose(dir);
  } else {
    p_stdout("test2d: p_dopen cant find directory\n");
  }
}

static void
cd_act(char *args)
{
  char *home = p_getenv("HOME");
  char *path = strtok(args, " \t\r\n");
  if (!path) path = home;
  if (!path) p_stdout("test2d: p_getenv cant find HOME\n");
  else if (p_chdir(path)) p_stdout("test2d: p_chdir failed\n");
}

static void
pwd_act(char *args)
{
  char *path = p_getcwd();
  if (!path) p_stdout("test2d: p_getcwd failed\n");
  else {
    p_stdout(path);
    p_stdout("\n");
  }
}

static void
mv_act(char *args)
{
  char *path1 = strtok(args, " \t\r\n");
  char *path2 = strtok((char *)0, " \t\r\n");
  if (!path1 || !path2)
    p_stdout("test2d: syntax: mv old_path new_path\n");
  else if (p_rename(path1, path2))
    p_stdout("test2d: p_rename failed\n");
}

static void
rm_act(char *args)
{
  char *path = strtok(args, " \t\r\n");
  if (!path) p_stdout("test2d: syntax: rm path\n");
  else if (p_remove(path)) p_stdout("test2d: p_remove failed\n");
}

static void
mkdir_act(char *args)
{
  char *path = strtok(args, " \t\r\n");
  if (!path) p_stdout("test2d: syntax: mkdir path\n");
  else if (p_mkdir(path)) p_stdout("test2d: p_mkdir failed\n");
}

static void
rmdir_act(char *args)
{
  char *path = strtok(args, " \t\r\n");
  if (!path) p_stdout("test2d: syntax: rmdir path\n");
  else if (p_rmdir(path)) p_stdout("test2d: p_rmdir failed\n");
}
#endif

#if TRY_GRAPHICS
void
dis_graphics(void)
{
  if (win4) p_destroy(win4);
  if (win3) p_destroy(win3);
  if (win2) p_destroy(win2);
  if (win1) p_destroy(win1);
  if (offscr) p_destroy(offscr);
  if (scremote) p_disconnect(scremote);
  if (scr2) p_disconnect(scr2);
  if (scr) p_disconnect(scr);
}

int phints = 0;

#define SIMDX 450
#define SIMDY 450

static int w0 = SIMDX;
static int h0 = SIMDY;
static int yoffscr, xsub, ysub, doffscr, animating, xpal, ypal;
static int win3_hi, win3_flag, win4_hi;

void
simple_create(p_win **winp)
{
  *winp = p_window(scr, SIMDX, SIMDY, "test2d window", P_BG, phints, winp);
# if !TRY_GUI
  simple_draw(scr, *winp, w0, h0);
# endif
}

# if TRY_OFFSCREEN
void
on_animate(void *context)
{
  p_win *w = context;
  p_clear(offscr);
  p_color(offscr, P_FG);
  p_pen(offscr, 1, P_SOLID);
  box_draw(offscr, xsub, yoffscr+ysub, xsub+15, yoffscr+ysub+15, 0);
  p_bitblt(w, 0, 16, offscr, xsub, ysub+16, xsub+16, ysub+h0-16);
  if (animating) p_set_alarm(0.05, &on_animate, context);
  xsub += 1;
  if (xsub>=16) xsub = 0;
  ysub += 1;
  if (ysub>=16) ysub = 0;
  if (!doffscr) {
    yoffscr += 4;
    if (yoffscr>h0-32) {
      yoffscr = h0-32;
      doffscr = 1;
    }
  } else {
    yoffscr -= 4;
    if (yoffscr<16) {
      yoffscr = 16;
      doffscr = 0;
    }
  }
}
# endif

void
simple_draw(p_scr *scr, p_win *w, int x0, int y0)
{
  int fontb;
  int fonth = p_txheight(scr, P_GUI_FONT, 12, &fontb);
  int ytxt = 2*fonth;
  int dx, dx1, dx2;

  p_clear(w);
  p_color(w, P_FG);
  p_pen(w, 1, P_SOLID);
  p_font(w, P_GUI_FONT, 12, 0);

  /* draw 16x16 boxes in each corner (tests origin and size) */
  box_draw(w,     0,     0,   15,   15, 0);
  box_draw(w,     0, y0-16,   15, y0-1, 0);
  box_draw(w, x0-16,     0, x0-1,   15, 0);
  box_draw(w, x0-16, y0-16, x0-1, y0-1, 0);
  p_text(w, 16, ytxt, "Corner boxes show four sides?", 29);
  ytxt += 2*fonth;

# if TRY_OFFSCREEN
  if (offscr) {
    xsub = ysub = 0;
    yoffscr = 16;
    doffscr = 0;
    p_clr_alarm(&on_animate, (void *)0);
    animating = 0;
    on_animate(w);
    p_color(w, P_FG);
    p_pen(w, 1, P_SOLID);
  }
# endif

  /* string alignment test */
  dx = p_txwidth(scr, "String aligned within boxes?", 28, P_GUI_FONT, 12);
  box_draw(w, 16, ytxt-fontb+fonth, 16+dx, ytxt-fontb, 0);
  dx1 = p_txwidth(scr, "String ", 7, P_GUI_FONT, 12);
  dx2 = p_txwidth(scr, "aligned withinGARBAGE", 14, P_GUI_FONT, 12);
  box_draw(w, 16+dx1, ytxt-fontb, 16+dx1+dx2, ytxt, 0);
  p_text(w, 16, ytxt, "String aligned within boxes?GARBAGE", 28);
  ytxt += 2*fonth;

  /* fill alignment test */
  /* alignment ticks */
  seg_draw(w, x0-40,  9, x0-40, 19);
  seg_draw(w, x0-40, 61, x0-40, 71);
  seg_draw(w, x0- 9, 40, x0-19, 40);
  seg_draw(w, x0-61, 40, x0-71, 40);
  /* upper right corner */
  box_draw(w, x0-19, 30, x0-29, 40, 1);
  box_draw(w, x0-39, 30, x0-29, 20, 1);
  tri_draw(w, x0-29, 20, x0-29, 30, x0-19, 30, 1);
  tri_draw(w, x0-29, 40, x0-39, 30, x0-39, 40, 1);
  /* lower right, y --> 81-y */
  box_draw(w, x0-19, 51, x0-29, 41, 1);
  box_draw(w, x0-39, 51, x0-29, 61, 1);
  tri_draw(w, x0-29, 61, x0-29, 51, x0-19, 51, 1);
  tri_draw(w, x0-29, 41, x0-39, 51, x0-39, 41, 1);
  /* lower left, x0-x --> 81-(x0-x) */
  box_draw(w, x0-60, 51, x0-50, 41, 1);
  box_draw(w, x0-40, 51, x0-50, 61, 1);
  tri_draw(w, x0-50, 61, x0-50, 51, x0-60, 51, 1);
  tri_draw(w, x0-50, 41, x0-40, 51, x0-40, 41, 1);
  /* upper left, y --> 81-y */
  box_draw(w, x0-60, 30, x0-50, 40, 1);
  box_draw(w, x0-40, 30, x0-50, 20, 1);
  tri_draw(w, x0-50, 20, x0-50, 30, x0-60, 30, 1);
  tri_draw(w, x0-50, 40, x0-40, 30, x0-40, 40, 1);
  p_text(w, 16, ytxt, "Fills aligned (UR corner)?", 26);
  ytxt += 2*fonth;
  p_flush(w);

# if TRY_DOTS
  {
    int x[4], y[4];
    x[0] = x0-40; x[1] = x0-40; x[2] = x0-61; x[3] = x0-21;
    y[0] = 80; y[1] = 120; y[2] = 100; y[3] = 100;
    p_i_pnts(w, x, y, 4);
    p_dots(w);
  }
  seg_draw(w, x0-40,  82, x0-40, 118);
  seg_draw(w, x0-59, 100, x0-23, 100);
  p_text(w, 16, ytxt, "Dots aligned (.+.)?", 19);
  ytxt += 2*fonth;
# endif

# if TRY_SEGMENTS
  {
    int x[4], y[4];
    x[0] = x0-40; x[1] = x0-40; x[2] = x0-45; x[3] = x0-35;
    y[0] = 140; y[1] = 150; y[2] = 145; y[3] = 145;
    p_i_pnts(w, x, y, 4);
    p_segments(w);
  }
  seg_draw(w, x0-40, 130, x0-40, 138);
  seg_draw(w, x0-40, 152, x0-40, 160);
  seg_draw(w, x0-33, 145, x0-25, 145);
  seg_draw(w, x0-47, 145, x0-55, 145);
  p_text(w, 16, ytxt, "Segments aligned (-+-)?", 23);
  ytxt += 2*fonth;
# endif

# if TRY_ELLIPSE
  box_draw(w, x0-61, 170, x0-42, 189, 0);
  box_draw(w, x0-40, 170, x0-21, 189, 0);
  p_ellipse(w, x0-59, 172, x0-44, 187, 0);
  p_ellipse(w, x0-38, 172, x0-23, 187, 1);
  p_text(w, 16, ytxt, "Circs/rects OK (1 pixel gaps)?", 30);
  ytxt += 2*fonth;
# endif

# if TRY_RECT
  box_draw(w, x0-61, 200, x0-21, 230, 0);
  p_pen(w, 4, P_SOLID | P_SQUARE);
  p_rect(w, x0-57, 204, x0-24, 227, 1);
  p_pen(w, 3, P_SOLID | P_SQUARE);
  p_rect(w, x0-53, 208, x0-29, 222, 1);
  p_pen(w, 1, P_SOLID);
  p_rect(w, x0-50, 211, x0-32, 219, 1);
  p_rect(w, x0-48, 213, x0-33, 218, 0);
# endif

  ytxt -= fonth;
  if (ytxt<190) ytxt = 190;
  /* top set of curves tests thin linestyles */
  curve_draw(w, 20, ytxt, -50, -50);
  p_pen(w, 1, P_DASH);
  curve_draw(w, 40, ytxt, 50, 50);
  p_pen(w, 1, P_DOT);
  curve_draw(w, 110, ytxt, -50, 50);
  p_pen(w, 1, P_DASHDOT);
  curve_draw(w, 130, ytxt, 50, -50);
  p_pen(w, 1, P_DASHDOTDOT);
  curve_draw(w, 200, ytxt, -50, -50);
  p_pen(w, 1, P_SOLID);
  curve_draw(w, 220, ytxt, 50, 50);

# if TRY_CLIP
  p_clip(w, 315,ytxt+15, 345,ytxt+45);
  box_draw(w, 300,ytxt, 330,ytxt+30, 1);
  box_draw(w, 330,ytxt+30, 360,ytxt+60, 1);
  p_clip(w, 0,0, 0,0);
  box_draw(w, 313,ytxt+13, 346,ytxt+46, 0);
# endif
  ytxt += 70;

  /* telescope segments test line widths, should be centered on y=ytxt */
  seg_draw(w, 20, ytxt, 60, ytxt);
  p_pen(w, 2, P_SOLID);
  seg_draw(w, 60, ytxt, 100, ytxt);
  p_pen(w, 3, P_SOLID);
  seg_draw(w, 100, ytxt, 140, ytxt);
  p_pen(w, 4, P_SOLID);
  seg_draw(w, 140, ytxt, 180, ytxt);
  p_pen(w, 5, P_SOLID);
  seg_draw(w, 180, ytxt, 220, ytxt);
  p_pen(w, 6, P_SOLID);
  seg_draw(w, 220, ytxt, 260, ytxt);
  p_pen(w, 7, P_SOLID);
  seg_draw(w, 260, ytxt, 300, ytxt);
  p_pen(w, 8, P_SOLID);
  seg_draw(w, 300, ytxt, 340, ytxt);
  p_pen(w, 9, P_SOLID);
  seg_draw(w, 340, ytxt, 380, ytxt);
  p_pen(w, 5, P_SOLID);
  ytxt += 20;

  /* lower set of curves test thick linestyles */
  curve_draw(w, 20, ytxt, -50, -50);
  p_pen(w, 5, P_DASH);
  curve_draw(w, 40, ytxt, 50, 50);
  p_pen(w, 5, P_DOT);
  curve_draw(w, 110, ytxt, -50, 50);
  p_pen(w, 5, P_DASHDOT);
  curve_draw(w, 130, ytxt, 50, -50);
  p_pen(w, 5, P_DASHDOTDOT);
  curve_draw(w, 200, ytxt, -50, -50);
  p_pen(w, 5, P_SOLID);
  curve_draw(w, 220, ytxt, 50, 50);

  /* triangles test square joins --
   * upper one should have sharp points, no notch at closure,
   * lowere one should have rounded corners */
  p_pen(w, 5, P_SOLID | P_SQUARE);
  tri_draw(w, 300,ytxt, 340,ytxt, 300,ytxt+40, 0);
  p_pen(w, 5, P_SOLID);
  tri_draw(w, 350,ytxt+50, 350,ytxt+10, 310,ytxt+50, 0);
  ytxt += 60+fonth;

  hcour = p_txheight(scr, P_COURIER | P_BOLD, 14, &bcour);
  wcour = p_txwidth(scr, "M", 1, P_COURIER | P_BOLD, 14);
# if TRY_GUI
  ygui = ytxt;
  ytxt = ygui+fonth;
  p_font(w, P_COURIER | P_BOLD, 14, 0);
  p_text(w, 20, ytxt, "on_focus", 8);
  xfoc = 20+3*wcour;
  p_text(w, xfoc+7*wcour, ytxt, "on_key", 6);
  xkey = xfoc+6*wcour;
  p_text(w, xkey+10*wcour, ytxt, "on_click", 8);
  xclick = xkey+11*wcour;
  yhi0 = ytxt-bcour;
  yhi1 = yhi0+hcour;
  xhi0 = xclick-wcour;
  xhi1 = xhi0 + p_txwidth(scr, "on_click", 8, P_COURIER | P_BOLD, 14);
  hion = 0;
  p_text(w, xclick+9*wcour, ytxt, "on_motion", 9);
  xmotion = xclick+9*wcour;
# endif
}

void
seg_draw(p_win *w, int x0, int y0, int x1, int y1)
{
  int x[2], y[2];
  x[0] = x0;  x[1] = x1;  y[0] = y0;  y[1] = y1;
  p_i_pnts(w, x, y, 2);
  p_lines(w);
}

void
box_draw(p_win *w, int x0, int y0, int x1, int y1, int fillit)
{
  int x[5], y[5];
  x[0] = x0;  x[1] = x1;  x[2] = x1;  x[3] = x0;  x[4] = x0;
  y[0] = y0;  y[1] = y0;  y[2] = y1;  y[3] = y1;  y[4] = y0;
  p_i_pnts(w, x, y, 5);
  if (!fillit) p_lines(w);
  else         p_fill(w, 2);
}

void
tri_draw(p_win *w, int x0,int y0, int x1,int y1, int x2,int y2, int fillit)
{
  int x[4], y[4];
  x[0] = x0;  x[1] = x1;  x[2] = x2;  x[3] = x0;
  y[0] = y0;  y[1] = y1;  y[2] = y2;  y[3] = y0;
  p_i_pnts(w, x, y, 4);
  if (!fillit) p_lines(w);
  else         p_fill(w, 2);
}

extern double cv[6];
double cv[6] = { 1., 1., 0.896575, 0.517638, 0.267949, 0. };

void
curve_draw(p_win *w, int x0, int y0, int width, int height)
{
  int x[6], y[6], i;
  if (width<0) x0 -= width;
  if (height<0) y0 -= height;
  for (i=0 ; i<6 ; i++) {
    x[i] = x0 + (int)(cv[i]*width);
    y[i] = y0 + (int)(cv[5-i]*height);
  }
  p_i_pnts(w, x, y, 6);
  p_lines(w);
}

void
advanced_draw(p_win *w)
{
  int x0, y0;
  int xtxt = 20;
  int ytxt = 10;
  int base, height, width;
  int i, j, ii, jj;

  p_clear(win2);
# if TRY_PALETTE
  reset_palette();
# endif
  p_color(w, P_FG);

  /* cant tell what order simple_draw and advanced_draw run... */
  hcour = p_txheight(scr, P_COURIER | P_BOLD, 14, &bcour);
  wcour = p_txwidth(scr, "M", 1, P_COURIER | P_BOLD, 14);

  ytxt += p_txheight(scr, P_COURIER, 14, &base);
  p_font(w, P_COURIER, 14, 0);
  p_text(w, xtxt,ytxt, "Courier 14", 10);
  xtxt += p_txwidth(scr, "Courier 14", 10, P_COURIER, 14);
  p_font(w, P_COURIER | P_BOLD, 14, 0);
  p_text(w, xtxt,ytxt, " Bold", 5);
  xtxt += p_txwidth(scr, " Bold", 5, P_COURIER | P_BOLD, 14);
  p_font(w, P_COURIER | P_ITALIC, 14, 0);
  p_text(w, xtxt,ytxt, " Italic", 7);
  xtxt += p_txwidth(scr, " Italic", 7, P_COURIER | P_ITALIC, 14);
  p_font(w, P_COURIER | P_BOLD | P_ITALIC, 14, 0);
  p_text(w, xtxt,ytxt, " Bold-Italic", 12);

  ytxt += 6 + p_txheight(scr, P_HELVETICA, 10, &base);
  xtxt = 20;
  p_font(w, P_HELVETICA, 10, 0);
  p_text(w, xtxt,ytxt, "Helvetica 10", 12);
  xtxt += p_txwidth(scr, "Helvetica 10", 12, P_HELVETICA, 10);
  p_font(w, P_HELVETICA | P_BOLD, 10, 0);
  p_text(w, xtxt,ytxt, " Bold", 5);
  xtxt += p_txwidth(scr, " Bold", 5, P_HELVETICA | P_BOLD, 10);
  p_font(w, P_HELVETICA | P_ITALIC, 10, 0);
  p_text(w, xtxt,ytxt, " Italic", 7);
  xtxt += p_txwidth(scr, " Italic", 7, P_HELVETICA | P_ITALIC, 10);
  p_font(w, P_HELVETICA | P_BOLD | P_ITALIC, 10, 0);
  p_text(w, xtxt,ytxt, " Bold-Italic", 12);

  ytxt += 6 + p_txheight(scr, P_TIMES, 12, &base);
  xtxt = 20;
  p_font(w, P_TIMES, 12, 0);
  p_text(w, xtxt,ytxt, "Times 12", 8);
  xtxt += p_txwidth(scr, "Times 12", 8, P_TIMES, 12);
  p_font(w, P_TIMES | P_BOLD, 12, 0);
  p_text(w, xtxt,ytxt, " Bold", 5);
  xtxt += p_txwidth(scr, " Bold", 5, P_TIMES | P_BOLD, 12);
  p_font(w, P_TIMES | P_ITALIC, 12, 0);
  p_text(w, xtxt,ytxt, " Italic", 7);
  xtxt += p_txwidth(scr, " Italic", 7, P_TIMES | P_ITALIC, 12);
  p_font(w, P_TIMES | P_BOLD | P_ITALIC, 12, 0);
  p_text(w, xtxt,ytxt, " Bold-Italic", 12);

  ytxt += 6 + p_txheight(scr, P_NEWCENTURY, 18, &base);
  xtxt = 20;
  p_font(w, P_NEWCENTURY, 18, 0);
  p_text(w, xtxt,ytxt, "Newcentury 18", 13);
  xtxt += p_txwidth(scr, "Newcentury 18", 13, P_NEWCENTURY, 18);
  p_font(w, P_NEWCENTURY | P_BOLD, 18, 0);
  p_text(w, xtxt,ytxt, " Bold", 5);
  xtxt += p_txwidth(scr, " Bold", 5, P_NEWCENTURY | P_BOLD, 18);
  p_font(w, P_NEWCENTURY | P_ITALIC, 18, 0);
  p_text(w, xtxt,ytxt, " Italic", 7);
  xtxt += p_txwidth(scr, " Italic", 7, P_NEWCENTURY | P_ITALIC, 18);
  p_font(w, P_NEWCENTURY | P_BOLD | P_ITALIC, 18, 0);
  p_text(w, xtxt,ytxt, " Bold-Italic", 12);

  ytxt += 6 + p_txheight(scr, P_SYMBOL, 14, &base);
  xtxt = 20;
  p_font(w, P_SYMBOL, 14, 0);
  p_text(w, xtxt,ytxt, "Symbol 14", 9);
  xtxt += p_txwidth(scr, "Symbol 14", 9, P_SYMBOL, 14);

  xtxt += 20;
  p_pen(w, 1, P_SOLID);
  p_font(w, P_COURIER | P_BOLD, 14, 0);
  p_text(w, xtxt,ytxt, "Square", 6);
  height = p_txheight(scr, P_COURIER | P_BOLD, 14, &base);
  width = p_txwidth(scr, "Square", 6, P_COURIER | P_BOLD, 14);
  x0 = xtxt;
  y0 = ytxt -= base;
  base = height-base;  /* descent */
  box_draw(w, xtxt,ytxt, xtxt+width,ytxt+height, 0);
  p_font(w, P_COURIER | P_BOLD, 14, 3);
  p_text(w, xtxt+width+base,ytxt, "Square", 6);
  box_draw(w, xtxt+width,ytxt, xtxt+width+height,ytxt+width, 0);
  p_font(w, P_COURIER | P_BOLD, 14, 2);
  p_text(w, xtxt+width+height,ytxt+width+base, "Square", 6);
  box_draw(w, xtxt+width+height,ytxt+width, xtxt+height,ytxt+width+height, 0);
  p_font(w, P_COURIER | P_BOLD, 14, 1);
  p_text(w, xtxt+height-base,ytxt+width+height, "Square", 6);
  box_draw(w, xtxt+height,ytxt+width+height, xtxt,ytxt+height, 0);

  p_font(w, P_COURIER | P_BOLD, 14, 0);
  xtxt += width+3*height;
  width /= 6;
  p_color(w, P_WHITE);
  box_draw(w, xtxt, ytxt, xtxt+16*width, ytxt+6*height, 1);
  p_color(w, P_BLACK);
  p_text(w, xtxt+width,ytxt+height-base, "black", 5);
  box_draw(w, xtxt, ytxt+height, xtxt+7*width, ytxt+2*height, 1);
  p_color(w, P_WHITE);
  p_text(w, xtxt+width,ytxt+2*height-base, "white", 5);
  p_color(w, P_RED);
  p_text(w, xtxt+8*width,ytxt+height-base, "red", 3);
  p_color(w, P_GREEN);
  p_text(w, xtxt+8*width,ytxt+2*height-base, "green", 5);
  p_color(w, P_BLUE);
  p_text(w, xtxt+8*width,ytxt+3*height-base, "blue", 4);
  p_color(w, P_CYAN);
  p_text(w, xtxt+8*width,ytxt+4*height-base, "cyan", 4);
  p_color(w, P_MAGENTA);
  p_text(w, xtxt+8*width,ytxt+5*height-base, "magenta", 7);
  p_color(w, P_YELLOW);
  p_text(w, xtxt+8*width,ytxt+6*height-base, "yellow", 6);
  p_color(w, P_BG);
  box_draw(w, xtxt, ytxt+2*height, xtxt+7*width, ytxt+4*height, 1);
  p_color(w, P_FG);
  p_text(w, xtxt+width,ytxt+3*height-base, "fg", 2);
  p_text(w, xtxt+width,ytxt+4*height-base, "xorfg", 5);
  p_color(w, P_XOR);
  box_draw(w, xtxt, ytxt+3*height, xtxt+7*width, ytxt+4*height, 1);
  ytxt += 6*height + 4;

# if TRY_NDXCELL
  {
    unsigned char cells[1200];
    p_col_t fgpix, bgpix, pix;
    for (j=0 ; j<400 ; j+=80)
      for (i=0 ; i<20 ; i+=4)
        for (jj=0 ; jj<80 ; jj+=20)
          for (ii=0 ; ii<4 ; ii++)
            cells[j+i+jj+ii] = (unsigned char)(P_BG - (i/4+j/16)%10);
    p_color(w, P_FG);
    xtxt = 20;
    box_draw(w, xtxt-1, ytxt-1, xtxt+20, ytxt+20, 0);
    p_ndx_cell(w, cells, 20,20, xtxt,ytxt,xtxt+20,ytxt+20);
    xtxt += 30;
    p_color(w, P_FG);
    box_draw(w, xtxt-1, ytxt-1, xtxt+10, ytxt+10, 0);
    p_ndx_cell(w, cells, 20,20, xtxt,ytxt,xtxt+10,ytxt+10);
    xtxt += 20;
    p_color(w, P_FG);
    box_draw(w, xtxt-1, ytxt-1, xtxt+40, ytxt+40, 0);
    p_ndx_cell(w, cells, 20,20, xtxt,ytxt,xtxt+40,ytxt+40);
    xtxt += 50;
    p_color(w, P_FG);
    box_draw(w, xtxt-1, ytxt-1, xtxt+10, ytxt+40, 0);
    p_ndx_cell(w, cells, 20,20, xtxt,ytxt,xtxt+10,ytxt+40);
    xtxt += 20;
    p_color(w, P_FG);
    box_draw(w, xtxt-1, ytxt-1, xtxt+40, ytxt+10, 0);
    p_ndx_cell(w, cells, 20,20, xtxt,ytxt,xtxt+40,ytxt+10);

    xtxt += 50;
    for (jj=0 ; jj<25 ; jj+=5)
      for (ii=0 ; ii<5 ; ii++)
        cells[jj+ii] = (unsigned char)(P_BG - (ii+jj)%10);
    box_draw(w, xtxt-1, ytxt-1, xtxt+40, ytxt+40, 0);
    p_ndx_cell(w, cells, 5,5, xtxt,ytxt,xtxt+40,ytxt+40);

#  if TRY_RGBREAD
    p_rgb_read(w, cells, x0,y0, x0+20,y0+20);
    fgpix = P_RGB(cells[0],cells[1],cells[2]);
    bgpix = P_RGB(cells[1197],cells[1198],cells[1199]);
    for (j=0 ; j<400 ; j+=20)
      for (i=0 ; i<20 ; i++) {
        pix = P_RGB(cells[3*(i+j)],cells[3*(i+j)+1],cells[3*(i+j)+2]);
        if (pix==fgpix) cells[i+j] = P_FG;
        else if (pix==bgpix) cells[i+j] = P_BG;
        else cells[i+j] = P_RED;
      }
    xtxt += 50;
    p_ndx_cell(w, cells, 20,20, xtxt,ytxt,xtxt+20,ytxt+20);
#  endif

    ytxt += 50;
  }
# endif

# if TRY_PALETTE
  p_color(w, P_FG);
  xtxt = 20;
  p_text(w, xtxt,ytxt+height-base, "palette:", 8);
  xpal = xtxt + 8*width + 10;
  ypal = ytxt;
  p_color(win2, P_BG);
  box_draw(win2, xpal-3*wcour,ypal+2*hcour-bcour,
           xpal-2*wcour,ypal+3*hcour-bcour, 1);
  p_color(win2, P_FG);
  p_font(win2, P_COURIER | P_BOLD, 14, 0);
  {
    char palno[8];
    palno[0] = pal_number+'0';
    palno[1] = '\0';
    p_text(win2, xpal-3*wcour,ypal+2*hcour, palno, 1);
    for (i=0 ; i<200 ; i++) {
      p_color(win2, i);
      box_draw(win2, xpal+i,ypal, xpal+i+1,ypal+40, 1);
    }
  }
# endif

#if TRY_RGBCELL
  rgbcell_on = !rgbcell_on;
  rgbcell_act(0);
#endif
}

# if TRY_PALETTE
static void
reset_palette(void)
{
  int i;
  p_col_t colors[200];
  if (!win2) return;
  if (pal_number==0) {
    p_palette(win2, colors, 0);
  } else {
    make_rainbow(pal_number, colors);
    p_palette(win2, colors, 200);
    for (i=0 ; i<200 ; i++) colors[i] = 0;  /* check for no effect */
  }
}
# endif

static void
make_rainbow(int pn, p_col_t *colors)
{
  long i, j, k;
  for (i=0 ; i<200 ; i++) {
    if (pn==1) {
      j = (255*i)/199;
      colors[i] = P_RGB(j,j,j);
    } else if (pn==2) {
      if (i<=33) {
        j = (251*(33-i))/66;
        if (j>127) k=(255*(255-j))/j, j=255;
        else k=255, j=(255*j)/(255-j);
        colors[i] = P_RGB(k,0,j);
      } else if (i<=100) {
        j = (251*(100-i))/66;
        if (j>127) k=(255*(255-j))/j, j=255;
        else k=255, j=(255*j)/(255-j);
        colors[i] = P_RGB(j,k,0);
      } else if (i<=167) {
        j = (251*(167-i))/66;
        if (j>127) k=(255*(255-j))/j, j=255;
        else k=255, j=(255*j)/(255-j);
        colors[i] = P_RGB(0,j,k);
      } else {
        j = (251*(233-i))/66;
        if (j>127) k=(255*(255-j))/j, j=255;
        else k=255, j=(255*j)/(255-j);
        colors[i] = P_RGB(k,0,j);
      }
    }
  }
}

static void
feep_act(char *args)
{
  if (win1) p_feep(win1);
}

static void
redraw_act(char *args)
{
  if (win1) {
    p_clear(win1);
    simple_draw(scr, win1, w0, h0);
  } else {
    win1 = p_window(scr, w0, h0, "test2d window", P_BG, phints, &win1);
  }
  if (!win2) {
    int hints = phints | P_NORESIZE | P_NOKEY | P_NOMOTION;
    if (rgbcell_on) hints |= P_RGBMODEL;
    win2 = p_window(scr, 400, 400, "test2d advanced", P_BG, hints, &win2);
  } else {
    advanced_draw(win2);
  }
}

static void
private_act(char *args)
{
  phints ^= P_PRIVMAP;
  if (win2) {
    int hints = phints | P_NORESIZE | P_NOKEY | P_NOMOTION;
    p_win *w = win2;
    win2 = 0;
    p_destroy(w);
    if (rgbcell_on) hints |= P_RGBMODEL;
    win2 = p_window(scr, 400, 400, "test2d advanced", P_BG, hints, &win2);
  }
# if TRY_STDINIT
  sprintf(p_wkspc.c, "test2d: using %s colormap\n",
          (phints&P_PRIVMAP)? "private" : "shared");
  p_stdout(p_wkspc.c);
# endif
}

static void
rgb_act(char *args)
{
  phints ^= P_RGBMODEL;
  if (win2) {
    p_win *w = win2;
    win2 = 0;
    p_destroy(w);
    win2 = p_window(scr, 400, 400, "test2d advanced", P_BG,
                   phints | P_NORESIZE | P_NOKEY | P_NOMOTION, &win2);
  }
# if TRY_STDINIT
  sprintf(p_wkspc.c, "test2d: using %s colormap\n",
          (phints&P_RGBMODEL)? "5-9-5 true" : "pseudo");
  p_stdout(p_wkspc.c);
# endif
}

# if TRY_RESIZE
static void
resize_act(char *args)
{
  char *arg2 = args;
  int w = (int)strtol(args, &arg2, 0);
  int h = (int)strtol(arg2, (char **)0, 0);
  if (win3) return;
  if (win1) p_raise(win1);
  if (win1 && w>50 && w<1000 && h>50 && w<1000) {
    p_resize(win1, w, h);
    if (w!=w0 || h!=h0) {
#  if TRY_OFFSCREEN
      if (offscr) {
        p_destroy(offscr);
        offscr = p_offscreen(win1, 32, h);
      }
#  endif
      simple_draw(scr, win1, w, h);
      w0 = w;
      h0 = h;
    }
  }
#  if TRY_STDINIT
  sprintf(p_wkspc.c, "test2d: %d pixels X %d pixels\n", w, h);
  p_stdout(p_wkspc.c);
#  endif
}
# endif

# if TRY_CLIPBOARD
static void
paste_act(char *args)
{
  if (!win1) return;
#  if TRY_STDINIT
  p_stdout("paste got:");
  p_stdout(p_spaste(win1));
  p_stdout("\n");
#  endif
}
# endif

# if TRY_METAFILE
static void
meta_act(char *args)
{
  char *path = strtok(args, " \t\r\n");
  if (path && win1) {
    p_win *meta = p_metafile(win1, path, 0, 0, w0, h0, 0);
    if (meta) {
      simple_draw(scr, meta, w0, h0);
      p_destroy(meta);
#  if TRY_STDINIT
      sprintf(p_wkspc.c, "test2d: wrote %s\n", path);
      p_stdout(p_wkspc.c);
#  endif
    } else {
#  if TRY_STDINIT
      p_stdout("test2d: p_metafile() not implemented\n");
#  endif
    }
  } else {
#  if TRY_STDINIT
    p_stdout("test2d: meta needs pathname and basic window\n");
#  endif
  }
}
# endif

# if TRY_PALETTE
static void
pal_act(char *args)
{
  pal_number++;
  if (pal_number>2) pal_number = 0;
  advanced_draw(win2);
}
# endif

# if TRY_RGBCELL
static void
rgbcell_act(char *args)
{
  rgbcell_on = !rgbcell_on;
  if (win2 && rgbcell_on) {
    int i, j, k;
    unsigned char *rgbs = p_malloc(3*200*50);
    p_col_t colors[200];
    make_rainbow(2, colors);
    for (i=j=0 ; i<600 ; i+=3,j++) {
      rgbs[i  ] = (unsigned char)P_R(colors[j]);
      rgbs[i+1] = (unsigned char)P_G(colors[j]);
      rgbs[i+2] = (unsigned char)P_B(colors[j]);
    }
    for (j=600,k=49 ; j<30000 ; j+=600,k--)
      for (i=0 ; i<600 ; i+=3) {
        rgbs[j+i  ] = (k*(int)rgbs[i  ])/50;
        rgbs[j+i+1] = (k*(int)rgbs[i+1])/50;
        rgbs[j+i+2] = (k*(int)rgbs[i+2])/50;
      }
    p_rgb_cell(win2, rgbs, 200, 50, xpal,ypal+50, xpal+200,ypal+100);
    p_free(rgbs);
  } else if (win2) {
    p_color(win2, P_BG);
    box_draw(win2, xpal,ypal+50, xpal+200,ypal+100, 1);
  }
}
# endif
#endif

#if TRY_GUI
void
on_focus(void *c,int in)
{
  p_win *w = *(p_win**)c;
  if (w!=win1 || win3 || !w) return;
  p_color(w, P_BG);
  box_draw(w, xfoc, ygui-bcour+hcour, xfoc+3*wcour, ygui-bcour, 1);
  p_color(w, P_FG);
  p_font(w, P_COURIER | P_BOLD, 14, 0);
  if (in) p_text(w, xfoc, ygui, "IN", 2);
  else p_text(w, xfoc, ygui, "OUT", 3);
}

void
on_key(void *c,int k,int md)
{
  char txt[8];
  p_win *w = *(p_win**)c;
  if (w!=win1 || win3) return;
  p_color(w, P_BG);
  box_draw(w, xkey, ygui-bcour+hcour, xkey+8*wcour, ygui-bcour, 1);
  p_color(w, P_FG);
  p_font(w, P_COURIER | P_BOLD, 14, 0);
  txt[0] = (md&P_SHIFT)? 'S' : ' ';
  txt[1] = (md&P_CONTROL)? 'C' : ' ';
  txt[2] = (md&P_META)? 'M' : ' ';
  txt[3] = (md&P_ALT)? 'A' : ' ';
  txt[4] = (md&P_COMPOSE)? '+' : ' ';
  txt[5] = (md&P_KEYPAD)? '#' : ' ';
  if (k=='\177') {
    txt[6] = '^';
    txt[7] = '?';
  } else if (k<P_LEFT) {
    txt[6] = (k<' ')? '^' : ' ';
    txt[7] = (k<' ')? k+'@' : k;
  } else if (k<=P_F12) {
    txt[6] = 'F';
    if (k==P_LEFT) txt[7] = '<';
    else if (k==P_RIGHT) txt[7] = '>';
    else if (k==P_UP) txt[7] = '^';
    else if (k==P_DOWN) txt[7] = 'V';
    else if (k==P_PGUP) txt[7] = 'U';
    else if (k==P_PGDN) txt[7] = 'D';
    else if (k==P_HOME) txt[7] = 'H';
    else if (k==P_END) txt[7] = 'E';
    else if (k==P_INSERT) txt[7] = 'I';
    else if (k<P_F10) txt[7] = (k-P_F0)+'0';
    else txt[7] = (k-P_F10)+'a';
  } else {
    txt[6] = '?';
    txt[7] = '?';
  }
  p_text(w, xkey, ygui, txt, 8);
}

static unsigned long multi_ms = 0;

void
on_click(void *c,int b,int md,int x,int y, unsigned long ms)
{
  char txt[64];
  p_win *w = *(p_win**)c;
  int down = (md&(1<<(b+2)))==0;
  if (w==win1 && !win3) {
    unsigned long dms = multi_ms? ms-multi_ms : 0;
    if (down) multi_ms = ms;
    p_color(w, P_BG);
    box_draw(w, xkey, ygui-bcour+hcour, xkey+8*wcour, ygui-bcour, 1);
    box_draw(w, xclick, ygui-bcour+hcour, xclick+7*wcour, ygui-bcour, 1);
    p_color(w, P_FG);
    p_font(w, P_COURIER | P_BOLD, 14, 0);
    txt[0] = (md&P_SHIFT)? 'S' : ' ';
    txt[1] = (md&P_CONTROL)? 'C' : ' ';
    txt[2] = (md&P_META)? 'M' : ' ';
    txt[3] = (md&P_ALT)? 'A' : ' ';
    txt[4] = (md&P_COMPOSE)? '+' : ' ';
    txt[5] = (md&P_KEYPAD)? '#' : ' ';
    txt[6] = '-';
    txt[7] = '>';
    p_text(w, xkey, ygui, txt, 8);
    sprintf(txt, "%03o %03o", md&0370, 1<<(b+2));
    p_text(w, xclick, ygui, txt, 7);
    if (down && dms && dms<250 && x>=xhi0 && x<xhi1 && y>=yhi0 && y<yhi1) {
      p_color(w, P_XOR);
      box_draw(w, xhi0, yhi0, xhi1, yhi1, 1);
      hion = !hion;
# if TRY_CLIPBOARD
      if (hion) p_scopy(w, "on_click", 8);
      else p_scopy(w, (char *)0, 0);
# endif
    }
# if TRY_OFFSCREEN
    if (y<=ygui && x<16 && !down) {
      animating = !animating;
      if (animating) p_set_alarm(0.05, &on_animate, w);
      else p_clr_alarm(&on_animate, (void *)0);
    }
# endif
# if TRY_MENU
    if (!win3 && down && dms && dms<250 && y<ygui) {
      int x0, y0;
      p_winloc(w, &x0, &y0);
      win3 = p_menu(scr, 12*wcour, 6*hcour, x0+x, y0+y, P_GRAYA, &win3);
      win3_hi = win3_flag = win4_hi = 0;
    }
# endif
  } else if (win3) {
    if (down || win3_flag) {
      p_win *w3 = win3;
      if (win3_hi!=2) {  /* win4 already destroyed */
        win3 = 0;
        p_destroy(w3);
# if TRY_STDINIT
        sprintf(p_wkspc.c, "test2d: menu item %d\n", win3_hi);
        p_stdout(p_wkspc.c);
        prompt_issued = 0;
# endif
      } else {
        p_win *w4 = win4;
        int skip = 0;
        if (!win4_hi) {
          if (w!=win3) {  /* get x,y relative to win3 */
            int x0, y0;
            p_winloc(w, &x0, &y0);
            x += x0;
            y += y0;
            p_winloc(win3, &x0, &y0);
            x -= x0;
            y -= y0;
          }
          skip = (x>=0 && x<12*wcour && y>=0 && y<6*hcour);
        }
        if (!skip) {
          win4 = win3 = 0;
          p_destroy(w4);
          p_destroy(w3);
# if TRY_STDINIT
          sprintf(p_wkspc.c, "test2d: submenu item %d\n", win4_hi);
          p_stdout(p_wkspc.c);
          prompt_issued = 0;
# endif
        }
      }
    }
  }
}

static int cur_cursor = 0;

void
on_motion(void *c,int md,int x,int y)
{
  char txt[64];
  p_win *w = *(p_win**)c;
  if (w==win1 && !win3) {
    p_color(w, P_BG);
    box_draw(w, xmotion, ygui-bcour+hcour, xmotion+12*wcour, ygui-bcour, 1);
    p_color(w, P_FG);
    p_font(w, P_COURIER | P_BOLD, 14, 0);
    sprintf(txt, "% 4d,% 4d  ", x, y);
    p_text(w, xmotion, ygui, txt, strlen(txt));
# if TRY_CURSOR
    {
      int i;
      if (y > yhi1) {
        i = ((P_NONE+1)*x)/w0;
        if (i<0) i = 0;
        if (i>P_NONE) i = P_NONE;
      } else {
        i = P_SELECT;
      }
      if (i!=cur_cursor) {
        cur_cursor = i;
        p_cursor(w, i);
      }
    }
# endif
  } else if (win3) {
    int in_menu = (x>=4 && x<12*wcour-4 && y>=4 && y<6*hcour-4);
    int old3hi = win3_hi;
    int old4hi = win4_hi;
    in_menu = in_menu && (w==win3 || w==win4);
    if (!in_menu) {
      int x0, y0;
      p_winloc(w, &x0, &y0);
      x0 += x;
      y0 += y;
      if (w!=win3) {
        p_winloc(win3, &x, &y);
        x = x0 - x;
        y = y0 - y;
        in_menu = (x>=4 && x<12*wcour-4 && y>=4 && y<6*hcour-4);
        if (in_menu) w = win3;
      }
      if (win4 && w!=win4) {
        p_winloc(win4, &x, &y);
        x = x0 - x;
        y = y0 - y;
        in_menu = (x>=4 && x<12*wcour-4 && y>=4 && y<6*hcour-4);
        if (in_menu) w = win4;
      }
    }
    if (in_menu) {
      int *pwin_hi = (w==win3)? &win3_hi : &win4_hi;
      if (y<2*hcour)      *pwin_hi = 1;
      else if (y<4*hcour) *pwin_hi = 2;
      else                *pwin_hi = 3;
    } else {
      if (!win4) win3_hi = 0;
      win4_hi = 0;
    }
    if (win3_hi!=old3hi) {
      if (old3hi) hilite(win3, 0, 2*(old3hi-1)*hcour,
                         12*wcour, 2*old3hi*hcour, 0);
      if (win3_hi) hilite(win3, 0, 2*(win3_hi-1)*hcour,
                          12*wcour, 2*win3_hi*hcour, 1);
    }
    if (win4 && win4_hi!=old4hi) {
      if (old4hi) hilite(win4, 0, 2*(old4hi-1)*hcour,
                         12*wcour, 2*old4hi*hcour, 0);
      if (win4_hi) hilite(win4, 0, 2*(win4_hi-1)*hcour,
                          12*wcour, 2*win4_hi*hcour, 1);
    }
    if (!win3_flag && win3_hi) win3_flag = 1;
    if (win3_hi==2) {
# if TRY_MENU
      if (!win4) {
        int x0, y0;
        p_winloc(win3, &x0, &y0);
        x0 += 12*wcour;
        y0 += 2*hcour;
        win4 = p_menu(scr, 12*wcour, 6*hcour, x0, y0, P_GRAYA, &win4);
        win4_hi = 0;
      }
# endif
    } else if (win4) {
      p_win *w4 = win4;
      win4 = 0;
      p_destroy(w4);
    }
  }
}

void
hilite(p_win *w, int x0, int y0, int x1, int y1, int onoff)
{
  p_color(w, onoff? P_WHITE : P_GRAYA);
  box_draw(w, x0, y0+4, x0+4, y1, 1);
  box_draw(w, x0, y0, x1, y0+4, 1);
  if (onoff) p_color(w, P_BLACK);
  box_draw(w, x0+4, y1-4, x1, y1, 1);
  box_draw(w, x1-4, y0+4, x1, y1-4, 1);
  if (onoff) {
    tri_draw(w, x0, y1, x0+4, y1, x0+4, y1-4, 1);
    tri_draw(w, x1-4, y0+4, x1, y0+4, x1, y0, 1);
  }
}

void
on_expose(void *c, int *xy)
{
  p_win **w = c;
  if (w == &win1) {
# if TRY_OFFSCREEN
    if (!offscr) offscr = p_offscreen(*w, 32, h0);
# endif
    simple_draw(scr, *w, w0, h0);
  } else if (w == &win2) {
    advanced_draw(*w);
  } else if (w == &win3) {
    p_win *ww = *w;
    p_color(ww, P_BLACK);
    p_font(ww, P_COURIER | P_BOLD, 14, 0);
    p_text(ww, wcour, 2*hcour-bcour, "item no. 1", 10);
    p_text(ww, wcour, 4*hcour-bcour, "sub menu >", 10);
    p_text(ww, wcour, 6*hcour-bcour, "item no. 3", 10);
    win3_hi = 0;
  } else if (w == &win4) {
    p_win *ww = *w;
    p_color(ww, P_BLACK);
    p_font(ww, P_COURIER | P_BOLD, 14, 0);
    p_text(ww, wcour, 2*hcour-bcour, "subitem #1", 10);
    p_text(ww, wcour, 4*hcour-bcour, "subitem #2", 10);
    p_text(ww, wcour, 6*hcour-bcour, "subitem #3", 10);
    win4_hi = 0;
  }
}

void
on_deselect(void *c)
{
  if (hion) {
    p_win *w = *(p_win **)c;
    if (w != win1) return;
    p_color(w, P_XOR);
    box_draw(w, xhi0, yhi0, xhi1, yhi1, 1);
    hion = 0;
# if TRY_STDINIT
    p_stdout("test2d: on_deselect called");
    p_stdout("\n");
    prompt_issued = 0;
# endif
  }
}

void
on_destroy(void *c)
{
  p_win *w = *(p_win **)c;
  int in_menu = (w==win4 || w==win3);
  if (w==win4) win4 = 0;
  else if (w==win3) win3 = 0;
  else if (w==win2) win2 = 0;
  else if (w==win1) win1 = 0;
# if TRY_STDINIT
  if (!in_menu) {
    p_stdout("test2d: on_destroy called\n");
    prompt_issued = 0;
  }
# endif
}

void
on_resize(void *c,int w,int h)
{
  if (w!=w0 || h!=h0) {
    p_win **win = c;
    if (win != &win1) return;
# if TRY_OFFSCREEN
    if (offscr) {
      p_destroy(offscr);
      offscr = p_offscreen(*win, 32, h);
    }
# endif
    if (w>0 && h>0) {
      simple_draw(scr, *win, w, h);
      w0 = w;
      h0 = h;
    } else {
      w0 = h0 = 1;
    }
  }
}

static char *scr_name[] = {
  "<nil bug>", "<unrecognized bug>", "scr", "scr2", "scremote" };

void
on_panic(p_scr *screen)
{
  int scrn = 1;
  if (screen) {
    if (screen==scr) {
      scrn = 2;
      scr = 0;
      panic_count++;
    } else if (screen==scr2) {
      scrn = 3;
      scr2 = 0;
    } else if (screen==scremote) {
      scrn = 4;
      scremote = 0;
    }
    scrn = 1;
  }
# if TRY_STDINIT
  sprintf(p_wkspc.c, "test2d: on_panic called on screen %s\n",
          scr_name[scrn]);
  p_stdout(p_wkspc.c);
  prompt_issued = 0;
# endif
}
#endif
