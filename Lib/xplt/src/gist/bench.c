/*
 * BENCH.C
 *
 * $Id$
 *
 * Exercise GIST graphics
 *
 */
/*    Copyright (c) 1994.  The Regents of the University of California.
                    All rights reserved.  */

#include "pstdlib.h"
#include "pstdio.h"
#include "hlevel.h"

#include <stdio.h>
#include <string.h>
#include <math.h>

#ifndef NO_XLIB
#include "xfancy.h"
#else
#include "play.h"
/* Animated versions of low level tests use these directly--
   ordinarily they are used indirectly via hlevel.h, as in the animate
   command in on_stdin.  */
static int GxAnimate(Engine *engine, GpBox *viewport);
static int GxStrobe(Engine *engine, int clear);
static int GxDirect(Engine *engine);
static int GxAnimate(Engine *engine, GpBox *viewport)
{
  return 0;
}
static int GxStrobe(Engine *engine, int clear)
{
  return 0;
}
static int GxDirect(Engine *engine)
{
  return 0;
}
#endif

#define PRINTF p_stdout
#define PRINTF1(f,x) sprintf(p_wkspc.c,f,x);p_stdout(p_wkspc.c)
#define PRINTF2(f,x,y) sprintf(p_wkspc.c,f,x,y);p_stdout(p_wkspc.c)
#define PRINTF3(f,x,y,z) sprintf(p_wkspc.c,f,x,y,z);p_stdout(p_wkspc.c)

#define DPI 75

#define PI 3.14159265359

/* Mesh for benchmarks 0-4 */
#define NCOLORS 240
#define NFRAMES 50
#define AMPLITUDE 0.5
#define NMID 25
#define N (2*NMID+1)
#define DPHI (2.0*PI/8.3)
GpReal cosine[N], amplitude[NFRAMES];
GaQuadMesh mesh;
GpReal *xmesh=0, *ymesh=0;
int *ireg=0;
GpColor *cmesh=0;
int meshReg= 0;
int meshBnd= 0;

#define PALETTE_FILE "earth.gp"
int actualColors;
GpColorCell *palette;

/* Level function and vector function */
GpReal *z=0, *u=0, *v=0;
GpReal zLevels[10]= { -.9, -.7, -.5, -.3, -.1, .1, .3, .5, .7, .9 };

/* Lissajous test */
#define LISSPTS 400
GpReal xliss[LISSPTS], yliss[LISSPTS];
int animate= 0, markers= 0, markType= M_POINT;
#define LISSFRAMES 50
#define LISSSIZE 40.0
#define NA1 1
#define NB1 5
#define NA2 2
#define NB2 7
#define RC1 40.0
#define RC2 160.0
#define DPHASE (2*PI/(LISSFRAMES-1))
#define DTHETA (PI/(LISSFRAMES-1))
#define DELPHI (2*PI/(LISSPTS-1))

GpReal xSmall[7]= { 1.1, 1.2, 9.1, 9.5, 11.5, 12.5, 13.6 };
GpReal ySmall[7]= { 1.0, 10., 5.0, 9.0,  8.0,  5.0,  4.0 };
GpReal xqSmall[7];
GpReal yqSmall[7];

static void GenerateLiss(int i, GpReal rc, GpReal size, int na, int nb);

static void GenerateLiss(int i, GpReal rc, GpReal size, int na, int nb)
{
  double theta= i*DTHETA;
  double phase= i*DPHASE;
  double xoff= rc*cos(theta),  yoff= rc*sin(theta);
  int j;
  for (j=0 ; j<LISSPTS ; j++) {
    xliss[j]= xoff+size*cos(na*DELPHI*j);
    yliss[j]= yoff+size*sin(nb*DELPHI*j+phase);
  }
}

static void GenerateMesh(int firstPass, GpReal amplFrame);

static void GenerateMesh(int firstPass, GpReal amplFrame)
{
  int i, j, k;
  GpReal xsplit;

  if (firstPass) {
    for (i=k=0 ; i<N ; i++, k=0)
      for (j=0 ; j<N*N ; j+=N) ymesh[j+i]= k++/(N-1.0);
  }

  k= 0;
  for (j=0 ; j<N*N ; j+=N) {
    xsplit= 0.5*(1.0+amplFrame*cosine[k++]);
    for (i=0 ; i<NMID ; i++) xmesh[j+i]= i*xsplit/NMID;
    for (i=NMID ; i<N ; i++) xmesh[j+i]= xsplit+(i-NMID)*(1.0-xsplit)/NMID;
  }
}

static void GenerateColors(void);

static void GenerateColors(void)
{
  GpReal dxmin= (1.0-AMPLITUDE)/(2.0*NMID);
  GpReal dxmax= (1.0+AMPLITUDE)/(2.0*NMID);
  GpReal range= dxmax-dxmin;
  int i, j;

  for (j=N ; j<N*N ; j+=N) {
    cmesh[j]= 0;
    for (i=1 ; i<N ; i++) cmesh[j+i]=
      (GpColor)((actualColors-1)*(xmesh[j+i]-xmesh[j+i-1]-dxmin)/range);
  }
}

static void SetOddReg(int r);

static void SetOddReg(int r)
{
  int i, j;
  for (j=N ; j<N*N ; j+=N) {
    ireg[j]= 0;
    for (i=1 ; i<N ; i++) ireg[j+i]= 1;
    if (j<N*N/4)
      for (i=3*N/8 ; i<5*N/8 ; i++) ireg[j+i]= r;
    else if (j<N*N/2)
      for (i=N/4 ; i<3*N/4 ; i++) ireg[j+i]= r;
    else if (j<5*N*N/8)
      for (i=1 ; i<N/4 ; i++) ireg[j+i]= r;
    else if (j<7*N*N/8)
      for (i=N/2 ; i<N ; i++) ireg[j+i]= r;
    else
      for (i=3*N/8 ; i<5*N/8 ; i++) ireg[j+i]= r;
  }
}

static void TenseZ(void);

static void TenseZ(void)
{
  /* Set diagonal from (10,40) to (20,30) to -1, and diagonal from
     (25,15) to (35,25) to +1.  */
  int i;
  for (i=0 ; i<11 ; i++) {
    z[10+40*N - (N-1)*i]= -1.0;  /* 10+i+N*(40-i) */
    z[25+15*N + (N+1)*i]= +1.0;  /* 25+i+N*(15+i) */
  }
}

static GpReal zs[N];
#define WGT_C 0.50
#define WGT_E (1./3.)
#define WGT_I 0.25

static void RelaxZ(void);

static void RelaxZ(void)
{
  /* Average each point with its 4 nearest neighbors.  */
  int i, j;
  register GpReal zss;

  /* Do 1st row */
  zss= z[0];
  z[0]= (1.-2.*WGT_C)*zss + WGT_C*(z[1]+z[N]);
  zs[0]= zss;
  for (i=1 ; i<N-1 ; i++) {
    zss= z[i];
    z[i]= (1.-3.*WGT_E)*zss + WGT_E*(zs[i-1]+z[i+1]+z[i+N]);
    zs[i]= zss;
  }
  zss= z[i];
  z[i]= 0.5*zss + 0.25*(zs[i-1]+z[i+N]);
  zs[i]= zss;

  /* Do interior rows */
  for (j=N ; j<N*(N-1) ; j+=N) {
    zss= z[j];
    z[j]= (1.-3.*WGT_E)*zss + WGT_E*(zs[0]+z[j+1]+z[j+N]);
    zs[0]= zss;
    for (i=1 ; i<N-1 ; i++) {
      zss= z[j+i];
      z[j+i]= (1.-4.*WGT_I)*zss + WGT_I*(zs[i-1]+zs[i]+z[j+i+1]+z[j+i+N]);
      zs[i]= zss;
    }
    zss= z[j+i];
    z[j+i]= (1.-3.*WGT_E)*zss + WGT_E*(zs[i-1]+zs[i]+z[j+i+N]);
    zs[i]= zss;
  }

  /* Do last row */
  zss= z[j];
  z[j]= (1.-2.*WGT_C)*zss + WGT_C*(zs[0]+z[j+1]);
  zs[0]= zss;
  for (i=1 ; i<N-1 ; i++) {
    zss= z[j+i];
    z[j+i]= (1.-3.*WGT_E)*zss + WGT_E*(zs[i-1]+zs[i]+z[j+i+1]);
    zs[i]= zss;
  }
  zss= z[j+i];
  z[j+i]= (1.-2.*WGT_C)*zss + WGT_C*(zs[i-1]+zs[i]);
  zs[i]= zss;
}

GpReal onePixel= 1.00001/512.5;

extern void GenerateImage(int now, GpReal cx, GpReal cy,
                          GpReal *px, GpReal *py, GpReal *qx, GpReal *qy);

void GenerateImage(int now, GpReal cx, GpReal cy,
                   GpReal *px, GpReal *py, GpReal *qx, GpReal *qy)
{
  GpReal halfwid= 0.5*onePixel*(now+1)*N;
  if (now==0) {
    int i, j;
    for (j=N ; j<N*N ; j+=N)
      for (i=1 ; i<N ; i++)
        cmesh[j+i]= (GpColor)((actualColors-1)*(z[j+i]+1.0)*0.5);
  }
  *px= cx - halfwid;
  *py= cy - halfwid;
  *qx= cx + halfwid;
  *qy= cy + halfwid;
}

static void get_time(double *cpu, double *sys, double *wall);
static double initWall= -1.0e31;

static void
get_time(double *cpu, double *sys, double *wall)
{
  *cpu = p_cpu_secs(sys);
  *wall = p_wall_secs();
  if (initWall<-1.0e30) initWall = *wall;
  *wall -= initWall;
}

static GpBox unitBox= {0.0, 1.0, 0.0, 1.0};
static GpBox zoomBox= {0.25, 0.75, 0.25, 0.75};
static GpTransform meshTrans;

static int calibFont= T_TIMES;

Engine *engine= 0;

int clear= 1;        /* interframe clear flag */
int colorMode= 0;    /* colorMode for hcp device */
int noCopy= 0;       /* noCopy mode flags */

int noKeyboard= 0;
char *noKeyTest[]= { "plg1", "fma", "plm1", "fma", "plf1", "fma",
                     "plc1", "fma", "pli", "fma", "plt", "fma", "pls", "fma",
                     "pldj", "fma", "txin", "txout", "fma", "calib", "fma",
                     "cfont", "calib", "fma", "quit", (char *)0 };

extern int on_idle(void);
extern int on_quit(void);
extern void on_stdin(char *line);
extern void on_exception(int signal, char *errmsg);

static int prompt_issued = 0;

void
on_stdin(char *line)
{
  double user, system, wall, user0, system0, wall0;
  int i;
  char *token;
  prompt_issued = 0;

  GdRevertLimits(1);  /* undo mouse zooms, if any */

  token= strtok(line, ", \t\n");
  if (!token) return;
  if (strcmp(token, "q")==0 || strcmp(token, "quit")==0) {
    p_quit();
    return;
  }

/* ------- A and P level tests -------- */

  if (strcmp(token, "0")==0) {
    PRINTF("   Begin benchmark 0 -- GaMesh (direct).\n");
    get_time(&user0, &system0, &wall0);
    gistA.l.color= FG_COLOR;
    GpSetTrans(&meshTrans);
    for (i=0 ; i<NFRAMES ; i++) {
      GenerateMesh(i==0, amplitude[i]);
      if (clear) GpClear(0, ALWAYS);
      GaMesh(&mesh, meshReg, meshBnd, 0);
      GpFlush(engine);
    }
    get_time(&user, &system, &wall);
    PRINTF3("elapsed seconds: %f user, %f system, %f wall\n",
           user-user0, system-system0, wall-wall0);
    PRINTF1("Plots per wall second= %f\n", NFRAMES/(wall-wall0));

  } else if (strcmp(token, "1")==0) {
    PRINTF("   Begin benchmark 1 -- GaMesh (animated).\n");
    get_time(&user0, &system0, &wall0);
    gistA.l.color= FG_COLOR;
    GpSetTrans(&meshTrans);
    GxAnimate(engine, &meshTrans.viewport);
    for (i=0 ; i<NFRAMES ; i++) {
      GenerateMesh(i==0, amplitude[i]);
      GaMesh(&mesh, meshReg, meshBnd, 0);
      GxStrobe(engine, 1);
      GpFlush(engine);
    }
    GxDirect(engine);  /* turn off animation mode */
    get_time(&user, &system, &wall);
    PRINTF3("elapsed seconds: %f user, %f system, %f wall\n",
           user-user0, system-system0, wall-wall0);
    PRINTF1("Plots per wall second= %f\n", NFRAMES/(wall-wall0));

  } else if (strcmp(token, "2")==0) {
    PRINTF("   Begin benchmark 2 -- GaFillMesh (direct).\n");
    get_time(&user0, &system0, &wall0);
    GpSetTrans(&meshTrans);
    for (i=0 ; i<NFRAMES ; i++) {
      GenerateMesh(i==0, amplitude[i]);
      GenerateColors();
      if (clear) GpClear(0, ALWAYS);
      GaFillMesh(&mesh, meshReg, cmesh+N+1, N);
      GpFlush(engine);
    }
    get_time(&user, &system, &wall);
    PRINTF3("elapsed seconds: %f user, %f system, %f wall\n",
           user-user0, system-system0, wall-wall0);
    PRINTF1("Plots per wall second= %f\n", NFRAMES/(wall-wall0));

  } else if (strcmp(token, "3")==0) {
    PRINTF("   Begin benchmark 3 -- GaFillMesh (animated).\n");
    get_time(&user0, &system0, &wall0);
    GpSetTrans(&meshTrans);
    GxAnimate(engine, &meshTrans.viewport);
    for (i=0 ; i<NFRAMES ; i++) {
      GenerateMesh(i==0, amplitude[i]);
      GenerateColors();
      GaFillMesh(&mesh, meshReg, cmesh+N+1, N);
      GxStrobe(engine, 1);
      GpFlush(engine);
    }
    GxDirect(engine);  /* turn off animation mode */
    get_time(&user, &system, &wall);
    PRINTF3("elapsed seconds: %f user, %f system, %f wall\n",
           user-user0, system-system0, wall-wall0);
    PRINTF1("Plots per wall second= %f\n", NFRAMES/(wall-wall0));

  } else if (strcmp(token, "im")==0) {
    GpReal cx= 0.5*(meshTrans.window.xmin+meshTrans.window.xmax);
    GpReal cy= 0.5*(meshTrans.window.ymin+meshTrans.window.ymax);
    GpReal px, qx, py, qy;
    PRINTF("   Begin image benchmark.\n");
    get_time(&user0, &system0, &wall0);
    GpSetTrans(&meshTrans);
    GpClear(0, ALWAYS);
    for (i=0 ; i<10 ; i++) {
      GenerateImage(i, cx, cy, &px, &py, &qx, &qy);
      GpCells(px, py, qx, qy, N-1, N-1, N, cmesh+N+1);
      GpFlush(engine);
    }
    get_time(&user, &system, &wall);
    PRINTF3("elapsed seconds: %f user, %f system, %f wall\n",
           user-user0, system-system0, wall-wall0);
    PRINTF1("Plots per wall second= %f\n", 10./(wall-wall0));

  } else if (strcmp(token, "clr")==0) {
    if (clear) {
      clear= 0;
      PRINTF("   Toggle interframe clear (now off)\n");
    } else {
      clear= 1;
      PRINTF("   Toggle interframe clear (now on)\n");
    }

  } else if (strcmp(token, "clip")==0) {
    if (gistClip) {
      gistClip= 0;
      PRINTF("   Toggle floating point clip (now off)\n");
    } else {
      gistClip= 1;
      PRINTF("   Toggle floating point clip (now on)\n");
    }

  } else if (strcmp(token, "zoom")==0) {
    if (meshTrans.window.xmax<1.0) {
      meshTrans.window= unitBox;
      PRINTF("   Toggle zoom (now off)\n");
    } else {
      meshTrans.window= zoomBox;
      PRINTF("   Toggle zoom (now on)\n");
    }

  } else if (strcmp(token, "c")==0) {
    PRINTF("   Clear\n");
    GpClear(0, CONDITIONALLY);

/* ------- attribute and property toggles -------- */

  } else if (strcmp(token, "bnd")==0) {
    if (meshBnd) {
      meshBnd= 0;
      PRINTF("   Toggle mesh boundary (now off)\n");
    } else {
      meshBnd= 1;
      PRINTF("   Toggle mesh boundary (now on)\n");
    }

  } else if (strcmp(token, "odd0")==0) {
    PRINTF("   Set mesh odd region to 0\n");
    SetOddReg(0);
    if (meshReg==2) meshReg= 0;

  } else if (strcmp(token, "odd2")==0) {
    PRINTF("   Set mesh odd region to 2\n");
    SetOddReg(2);

  } else if (strcmp(token, "reg0")==0) {
    PRINTF("   Set mesh plot to region 0\n");
    meshReg= 0;

  } else if (strcmp(token, "reg1")==0) {
    PRINTF("   Set mesh plot to region 1\n");
    meshReg= 1;

  } else if (strcmp(token, "reg2")==0) {
    PRINTF("   Set mesh plot to region 2\n");
    meshReg= 2;

  } else if (strcmp(token, "cmode")==0) {
    if (colorMode) {
      colorMode= 0;
      PRINTF("   Toggle hcp color mode (now no dump)\n");
    } else {
      colorMode= 1;
      PRINTF("   Toggle hcp color mode (now dump)\n");
    }
    GpDumpColors(hcpDefault, colorMode);
    GpDumpColors(engine, colorMode);

  } else if (strcmp(token, "nocopy")==0) {
    if (noCopy) {
      noCopy= 0;
      PRINTF("   Toggle noCopy mode (now copies mesh)\n");
    } else {
      noCopy= NOCOPY_MESH|NOCOPY_COLORS|NOCOPY_UV|NOCOPY_Z;
      PRINTF("   Toggle noCopy mode (now no mesh copies)\n");
    }

  } else if (strcmp(token, "earth")==0) {
    PRINTF("OK\n");
    if (palette) p_free(palette);
    actualColors= GpReadPalette(engine, "earth.gp", &palette, 240);
    GpSetPalette(hcpDefault, palette, actualColors);

  } else if (strcmp(token, "stern")==0) {
    PRINTF("OK\n");
    if (palette) p_free(palette);
    actualColors= GpReadPalette(engine, "stern.gp", &palette, 240);
    GpSetPalette(hcpDefault, palette, actualColors);

  } else if (strcmp(token, "rainbow")==0) {
    PRINTF("OK\n");
    if (palette) p_free(palette);
    actualColors= GpReadPalette(engine, "rainbow.gp", &palette, 240);
    GpSetPalette(hcpDefault, palette, actualColors);

  } else if (strcmp(token, "heat")==0) {
    PRINTF("OK\n");
    if (palette) p_free(palette);
    actualColors= GpReadPalette(engine, "heat.gp", &palette, 240);
    GpSetPalette(hcpDefault, palette, actualColors);

  } else if (strcmp(token, "gray")==0) {
    PRINTF("OK\n");
    if (palette) p_free(palette);
    actualColors= GpReadPalette(engine, "gray.gp", &palette, 240);
    GpSetPalette(hcpDefault, palette, actualColors);

  } else if (strcmp(token, "yarg")==0) {
    PRINTF("OK\n");
    if (palette) p_free(palette);
    actualColors= GpReadPalette(engine, "yarg.gp", &palette, 240);
    GpSetPalette(hcpDefault, palette, actualColors);

  } else if (strcmp(token, "work")==0) {
    PRINTF("OK\n");
    GdReadStyle(ghDevices[0].drawing, "work.gs");

  } else if (strcmp(token, "boxed")==0) {
    PRINTF("OK\n");
    GdReadStyle(ghDevices[0].drawing, "boxed.gs");

  } else if (strcmp(token, "axes")==0) {
    PRINTF("OK\n");
    GdReadStyle(ghDevices[0].drawing, "axes.gs");

  } else if (strcmp(token, "wide")==0) {
    GhGetLines();
    if (gistA.l.width>1.0) {
      gistA.l.width= 1.0;
      PRINTF("   Toggle wide lines (now off)\n");
    } else {
      gistA.l.width= 8.0;
      PRINTF("   Toggle wide lines (now on)\n");
    }
    GhSetLines();

  } else if (strcmp(token, "linetype")==0) {
    GhGetLines();
    gistA.l.type++;
    if (gistA.l.type>L_DASHDOTDOT) gistA.l.type= L_SOLID;
    if (gistA.l.type==L_SOLID)
      PRINTF("   Toggle line type (now solid)\n");
    else if (gistA.l.type==L_DASH)
      PRINTF("   Toggle line type (now dashed)\n");
    else if (gistA.l.type==L_DOT)
      PRINTF("   Toggle line type (now dotted)\n");
    else if (gistA.l.type==L_DASHDOT)
      PRINTF("   Toggle line type (now dash-dot)\n");
    else if (gistA.l.type==L_DASHDOTDOT)
      PRINTF("   Toggle line type (now dash-dot-dot)\n");
    else
      PRINTF("   Toggle line type (now bad type?)\n");
    GhSetLines();

  } else if (strcmp(token, "closed")==0) {
    GhGetLines();
    if (gistA.dl.closed) {
      gistA.dl.closed= 0;
      PRINTF("   Toggle closed lines (now off)\n");
    } else {
      gistA.dl.closed= 1;
      PRINTF("   Toggle closed lines (now on)\n");
    }
    GhSetLines();

  } else if (strcmp(token, "smooth")==0) {
    GhGetLines();
    if (gistA.dl.smooth==4) {
      gistA.dl.smooth= 0;
      PRINTF("   Toggle smooth lines (now off)\n");
    } else if (gistA.dl.smooth==1) {
      gistA.dl.smooth= 4;
      PRINTF("   Toggle smooth lines (now max)\n");
    } else {
      gistA.dl.smooth= 1;
      PRINTF("   Toggle smooth lines (now min)\n");
    }
    GhSetLines();

  } else if (strcmp(token, "rays")==0) {
    GhGetLines();
    if (gistA.dl.rays) {
      gistA.dl.rays= 0;
      PRINTF("   Toggle rays (now off)\n");
    } else {
      gistA.dl.rays= 1;
      PRINTF("   Toggle rays (now on)\n");
    }
    GhSetLines();

  } else if (strcmp(token, "edges")==0) {
    GhGetFill();
    if (gistA.e.type!=L_NONE) {
      gistA.e.type= L_NONE;
      PRINTF("   Toggle edges on plf (now off)\n");
    } else {
      gistA.e.type= L_SOLID;
      PRINTF("   Toggle edges on plf (now on)\n");
    }
    GhSetFill();

  } else if (strcmp(token, "limits")==0) {
    char *suffix;
    GpReal limits[4], value;
    int i= 0;
    if (GdGetLimits()) {
      PRINTF("****** No drawing yet\n");
      return;
    }
    PRINTF("   Setting limits\n");
    limits[0]= gistD.limits.xmin;
    limits[1]= gistD.limits.xmax;
    limits[2]= gistD.limits.ymin;
    limits[3]= gistD.limits.ymax;
    while ((token= strtok(0, "=, \t\n"))) {
      if (strcmp(token, "nice")==0) {
        token= strtok(0, "= \t\n");
        if (!token) break;
        if (strcmp(token, "0")==0) gistD.flags&= ~D_NICE;
        else gistD.flags|= D_NICE;
      } else if (strcmp(token, "square")==0) {
        token= strtok(0, "= \t\n");
        if (!token) break;
        if (strcmp(token, "0")==0) gistD.flags&= ~D_SQUARE;
        else gistD.flags|= D_SQUARE;
      } else if (strcmp(token, "restrict")==0) {
        token= strtok(0, "= \t\n");
        if (!token) break;
        if (strcmp(token, "0")==0) gistD.flags&= ~D_RESTRICT;
        else gistD.flags|= D_RESTRICT;
      } else if (i<4) {
        if (strcmp(token, "e")==0) {
          gistD.flags|= (D_XMIN<<i);
        } else {
          value= strtod(token, &suffix);
          if (suffix!=token) limits[i]= value;
          else break;
          gistD.flags&= ~(D_XMIN<<i);
        }
        i++;
      }
    }
    gistD.limits.xmin= limits[0];
    gistD.limits.xmax= limits[1];
    gistD.limits.ymin= limits[2];
    gistD.limits.ymax= limits[3];
    GdSetLimits();
    GhBeforeWait();

  } else if (strcmp(token, "logxy")==0) {
    if (GdGetLimits()) {
      PRINTF("****** No drawing yet\n");
      return;
    }
    PRINTF("   Setting logxy\n");
    token= strtok(0, ", \t\n");
    if (token) {
      if (strcmp(token, "0")==0) gistD.flags&= ~D_LOGX;
      else gistD.flags|= D_LOGX;
      token= strtok(0, ", \t\n");
      if (token) {
        if (strcmp(token, "0")==0) gistD.flags&= ~D_LOGY;
        else gistD.flags|= D_LOGY;
      }
    }
    GdSetLimits();
    GhBeforeWait();

/* ------- D and H level tests -------- */

  } else if (strcmp(token, "plg")==0) {
    PRINTF("   Lissajous test\n");
    get_time(&user0, &system0, &wall0);
    for (i=0 ; i<LISSFRAMES ; i++) {
      GenerateLiss(i, RC1, LISSSIZE, NA1, NB1);
      GhGetLines();
      gistD.legend= "#1 Lissajous pattern";
      gistA.dl.marks= markers;
      GdLines(LISSPTS, xliss, yliss);
      GenerateLiss(i, RC2, LISSSIZE, NA2, NB2);
      gistD.legend= "#2 Lissajous pattern";
      gistA.m.type= 0;  /* default marker */
      GdLines(LISSPTS, xliss, yliss);
      GhFMA();
    }
    get_time(&user, &system, &wall);
    PRINTF3("elapsed seconds: %f user, %f system, %f wall\n",
           user-user0, system-system0, wall-wall0);
    PRINTF1("Plots per wall second= %f\n", LISSFRAMES/(wall-wall0));

  } else if (strcmp(token, "plm")==0) {
    PRINTF("   Begin mesh benchmark.\n");
    get_time(&user0, &system0, &wall0);
    GhGetLines();
    for (i=0 ; i<NFRAMES ; i++) {
      GenerateMesh(i==0, amplitude[i]);
      gistD.legend= "Mesh lines";
      GdMesh(noCopy, &mesh, meshReg, meshBnd, 0);
      GhFMA();
    }
    get_time(&user, &system, &wall);
    PRINTF3("elapsed seconds: %f user, %f system, %f wall\n",
           user-user0, system-system0, wall-wall0);
    PRINTF1("Plots per wall second= %f\n", NFRAMES/(wall-wall0));

  } else if (strcmp(token, "plf")==0) {
    PRINTF("   Begin filled mesh benchmark.\n");
    get_time(&user0, &system0, &wall0);
    for (i=0 ; i<NFRAMES ; i++) {
      GenerateMesh(i==0, amplitude[i]);
      GenerateColors();
      gistD.legend= "Filled mesh";
      GdFillMesh(noCopy, &mesh, meshReg, cmesh+N+1, N);
      GhFMA();
    }
    get_time(&user, &system, &wall);
    PRINTF3("elapsed seconds: %f user, %f system, %f wall\n",
           user-user0, system-system0, wall-wall0);
    PRINTF1("Plots per wall second= %f\n", NFRAMES/(wall-wall0));

  } else if (strcmp(token, "plv")==0) {
    GpReal scale= 0.02;
    PRINTF("   Begin vector benchmark.\n");
    get_time(&user0, &system0, &wall0);
    GhGetVectors();
    for (i=0 ; i<NFRAMES ; i++) {
      GenerateMesh(i==0, amplitude[i]);
      gistD.legend= "Vectors";
      GdVectors(noCopy, &mesh, meshReg, u, v, scale);
      GhFMA();
    }
    get_time(&user, &system, &wall);
    PRINTF3("elapsed seconds: %f user, %f system, %f wall\n",
           user-user0, system-system0, wall-wall0);
    PRINTF1("Plots per wall second= %f\n", NFRAMES/(wall-wall0));

  } else if (strcmp(token, "plc")==0) {
    PRINTF("   Begin contour benchmark.\n");
    get_time(&user0, &system0, &wall0);
    for (i=0 ; i<NFRAMES ; i++) {
      GenerateMesh(i==0, amplitude[i]);
      GhGetLines();
      gistD.legend= "Contours";
      gistA.dl.marks= markers;
      GdContours(noCopy, &mesh, meshReg, z, zLevels, 10);
      GhFMA();
    }
    get_time(&user, &system, &wall);
    PRINTF3("elapsed seconds: %f user, %f system, %f wall\n",
           user-user0, system-system0, wall-wall0);
    PRINTF1("Plots per wall second= %f\n", NFRAMES/(wall-wall0));

  } else if (strcmp(token, "pli")==0) {
    GpReal cx= 0.5*(meshTrans.window.xmin+meshTrans.window.xmax);
    GpReal cy= 0.5*(meshTrans.window.ymin+meshTrans.window.ymax);
    GpReal px, qx, py, qy;
    PRINTF("   Plot image.\n");
    GenerateImage(noCopy, cx, cy, &px, &py, &qx, &qy);
    gistD.legend= "Image";
    GdCells(px, py, qx, qy, N-1, N-1, N, cmesh+N+1);
    GhBeforeWait();

  } else if (strcmp(token, "plg1")==0) {
    int i;
    PRINTF("   Plot small graph.\n");
    GhGetLines();
    gistD.legend= "7 point graph";
    gistA.dl.marks= markers;
    GdLines(7, xSmall, ySmall);
    /* change curve for next pass */
    for (i=0 ; i<7 ; i++) ySmall[i]+= 1.0;
    GhBeforeWait();

  } else if (strcmp(token, "pls")==0) {
    int i;
    PRINTF("   Plot small scatter.\n");
    GhGetLines();
    gistD.legend= "7 point scatter";
    gistA.dl.marks= 1;
    gistA.l.type= L_NONE;
    gistA.m.type= markType++;
    if (markType>M_CROSS) {
      if (markType>='C') markType= M_POINT;
      else if (markType<'A') markType= 'A';
    }
    GdLines(7, xSmall, ySmall);
    /* change curve for next pass */
    for (i=0 ; i<7 ; i++) ySmall[i]+= 1.0;
    GhBeforeWait();

  } else if (strcmp(token, "pldj")==0) {
    int i;
    PRINTF("   Plot disjoint lines.\n");
    GhGetLines();
    gistD.legend= "7 point disjoint";
    for (i=0 ; i<7 ; i++) {
      xqSmall[i]= xSmall[i]+1.5;
      yqSmall[i]= ySmall[i]+0.3;
    }
    GdDisjoint(7, xSmall, ySmall, xqSmall, yqSmall);
    /* change curve for next pass */
    for (i=0 ; i<7 ; i++) ySmall[i]+= 1.0;
    GhBeforeWait();

  } else if (strcmp(token, "plc1")==0) {
    PRINTF("   Plot contours.\n");
    GenerateMesh(1, amplitude[NFRAMES-2]);
    GhGetLines();
    gistD.legend= "Contours";
    gistA.dl.marks= markers;
    GdContours(noCopy, &mesh, meshReg, z, zLevels, 10);
    GhBeforeWait();

  } else if (strcmp(token, "animate")==0) {
    if (animate) {
      animate= 0;
      GhFMAMode(2, 0);
      PRINTF("   Toggle animation mode (now off)\n");
    } else {
      animate= 1;
      GhFMAMode(2, 1);
      PRINTF("   Toggle animation mode (now on)\n");
    }
    GhBeforeWait();

  } else if (strcmp(token, "marks")==0) {
    if (markers) {
      markers= 0;
      PRINTF("   Toggle markers (now off)\n");
    } else {
      markers= 1;
      PRINTF("   Toggle markers (now on)\n");
    }

  } else if (strcmp(token, "legends")==0) {
    if (ghDevices[0].doLegends) {
      ghDevices[0].doLegends= 0;
      PRINTF("   Toggle hcp legends (now off)\n");
    } else {
      ghDevices[0].doLegends= 1;
      PRINTF("   Toggle hcp legends (now on)\n");
    }

  } else if (strcmp(token, "plm1")==0) {
    PRINTF("   GdMesh test\n");
    gistA.l.color= FG_COLOR;
    GhGetLines();
    GenerateMesh(1, amplitude[NFRAMES-2]);
    gistD.legend= "Test mesh";
    GdMesh(noCopy, &mesh, meshReg, meshBnd, 0);
    GhBeforeWait();

  } else if (strcmp(token, "plf1")==0) {
    PRINTF("   GdFillMesh test\n");
    GenerateMesh(1, amplitude[NFRAMES-2]);
    GenerateColors();
    gistD.legend= "Test filled mesh";
    GdFillMesh(noCopy, &mesh, meshReg, cmesh+N+1, N);
    GhBeforeWait();

  } else if (strcmp(token, "txout")==0) {
    PRINTF("   Exterior text\n");
    GhGetText();
    gistA.t.color= MAGENTA_COLOR;
    gistA.t.font= T_NEWCENTURY;
    gistA.t.height= 24.*ONE_POINT;
    gistA.t.opaque= 1;
    gistD.legend= 0;
    GdText(0.2, 0.5, "Outside", 0);
    GhBeforeWait();

  } else if (strcmp(token, "txin")==0) {
    PRINTF("   Interior text\n");
    GhGetText();
    gistA.t.color= GREEN_COLOR;
    gistA.t.font= T_NEWCENTURY;
    gistA.t.height= 10.*ONE_POINT;
    gistA.t.opaque= 1;
    gistD.legend= 0;
    GdText(0.5, 0.5, "Inside", 1);
    GhBeforeWait();

  } else if (strcmp(token, "calib")==0) {
    GpReal y0= 0.85;
    GpReal x[2], y[2], x2[2];
    PRINTF("   Text calibration frame\n");
    GhGetText();
    gistA.t.font= calibFont;
    gistA.t.height= 14.*ONE_POINT;
    gistA.t.alignH= TH_NORMAL;
    gistA.t.alignV= TV_NORMAL;
    gistA.t.opaque= 0;
    gistD.legend= 0;
    GdText(0.15, y0-0.0*gistA.t.height, "Times 14 pt gjpqy", 0);
    GdText(0.15, y0-1.0*gistA.t.height, "Times 14 pt EFTZB", 0);
    GdText(0.15, y0-2.0*gistA.t.height, "Third line, 14 pt", 0);
    GdSetSystem(-1);
    x[0]= x[1]= .30;
    y[0]= y0= 0.75;
    y[1]= y[0]+gistA.t.height;
    GhGetLines();
    GdLines(2, x, y);
    gistA.t.alignH= TH_LEFT;
    GdText(0.30, y0-1.0*gistA.t.height, "Horizontal (lf)", 0);
    gistA.t.alignH= TH_CENTER;
    GdText(0.30, y0-2.5*gistA.t.height, "Horizontal (cn)", 0);
    gistA.t.alignH= TH_RIGHT;
    GdText(0.30, y0-4.0*gistA.t.height, "Horizontal (rt)", 0);
    y[0]= y0-4.5*gistA.t.height;
    y[1]= y0-5.5*gistA.t.height;
    GdLines(2, x, y);
    gistA.t.alignH= TH_NORMAL;
    x[0]= 0.50;
    x[1]= 0.50-gistA.t.height;
    x2[0]= 0.60;
    x2[1]= 0.60+gistA.t.height;
    y[0]= y[1]= y0= .85;
    GdLines(2, x, y);
    GdLines(2, x2, y);
    gistA.t.alignV= TV_TOP;
    GdText(0.50, y0-0.0*gistA.t.height, "Vertical (tp)", 0);
    y[0]= y[1]= y0-2.0*gistA.t.height;
    GdLines(2, x, y);
    GdLines(2, x2, y);
    gistA.t.alignV= TV_CAP;
    GdText(0.50, y0-2.0*gistA.t.height, "Vertical (cp)", 0);
    y[0]= y[1]= y0-4.0*gistA.t.height;
    GdLines(2, x, y);
    GdLines(2, x2, y);
    gistA.t.alignV= TV_HALF;
    GdText(0.50, y0-4.0*gistA.t.height, "Vertical (hf)", 0);
    y[0]= y[1]= y0-6.0*gistA.t.height;
    GdLines(2, x, y);
    GdLines(2, x2, y);
    gistA.t.alignV= TV_BASE;
    GdText(0.50, y0-6.0*gistA.t.height, "Vertical (ba)", 0);
    y[0]= y[1]= y0-8.0*gistA.t.height;
    GdLines(2, x, y);
    GdLines(2, x2, y);
    gistA.t.alignV= TV_BOTTOM;
    GdText(0.50, y0-8.0*gistA.t.height, "Vertical (bt)", 0);
    GdSetSystem(0);
    GhBeforeWait();

  } else if (strcmp(token, "cfont")==0) {
    if (calibFont!=T_TIMES) {
      calibFont= T_TIMES;
      PRINTF("   Toggle calib font (now Times)\n");
    } else {
      calibFont= T_HELVETICA;
      PRINTF("   Toggle calib font (now Helvetica)\n");
    }

  } else if (strcmp(token, "fma")==0) {
    PRINTF("   Frame advance\n");
    GhFMA();

  } else if (strcmp(token, "hcp")==0) {
    PRINTF("   Sent frame to hcp\n");
    GhHCP();

  } else if (strcmp(token, "hcpon")==0) {
    PRINTF("   hcp on fma\n");
    GhFMAMode(1, 2);

  } else if (strcmp(token, "hcpoff")==0) {
    PRINTF("   NO hcp on fma\n");
    GhFMAMode(0, 2);

  } else if (strcmp(token, "redraw")==0) {
    PRINTF("   Redraw\n");
    GhRedraw();

  } else if (strcmp(token, "plt")==0) {
    /* Test text primitive */
    GpReal ynow= 0.90;
    PRINTF("   Text test\n");
    GhGetText();
    gistA.t.color= FG_COLOR;
    gistA.t.font= T_HELVETICA;
    gistA.t.height= ONE_POINT*12.0;
    GdText(0.25, ynow, "Helvetica 12 pt?", 0);
    ynow-= gistA.t.height;
    GdText(0.25, ynow, "Helvetica 12 pt, second line", 0);
    gistA.t.height= ONE_POINT*18.0;
    ynow-= gistA.t.height;
    GdText(0.25, ynow, "Helvetica 18 pt?", 0);
    gistA.t.font= T_NEWCENTURY;
    gistA.t.height= ONE_POINT*12.0;
    ynow-= gistA.t.height;
    GdText(0.25, ynow, "New Century 12 pt, first line", 0);
    ynow-= gistA.t.height;
    gistA.t.font|= T_BOLD | T_ITALIC;
    GdText(0.25, ynow, "New Century 12 pt, bold italic", 0);
    gistA.t.font= T_NEWCENTURY;
    gistA.t.height= ONE_POINT*14.0;
    gistA.t.alignH= TH_RIGHT;
    ynow-= gistA.t.height;
    GdText(0.65, ynow, "New Century 14 pt,\nthree lines\nright justified", 0);
    ynow-= 3*gistA.t.height;
    gistA.t.alignH= TH_NORMAL;
    gistA.t.orient= TX_DOWN;
    GdText(0.25, ynow, "Path down", 0);
    GhBeforeWait();

  } else if (strcmp(token, "help")==0 || strcmp(token, "?")==0) {
    PRINTF("\n   This is the GIST benchmark/test program, commands are:\n");
    PRINTF("\n    Movies to test low level Gist functions--\n");
    PRINTF("0   - raw performance test GaMesh (direct)\n");
    PRINTF("1   - raw performance test GaMesh (animated)\n");
    PRINTF("2   - raw performance test GaFillMesh (direct)\n");
    PRINTF("3   - raw performance test GaFillMesh (animated)\n");
    PRINTF("clr - toggle interframe clear for tests 0-3\n");
    PRINTF("c   - clear screen\n");
    PRINTF("im  - raw performance test GpCells\n");
    PRINTF("clip  - toggle floating point clip\n");
    PRINTF("zoom  - toggle zoom for 0-3, im, pli\n");
    PRINTF("\n    Property toggles for high level tests--\n");
    PRINTF("bnd  - mesh boundary          odd0  - mesh region to 0\n");
    PRINTF("odd2  - mesh region to 2      reg0  - plot region to 0\n");
    PRINTF("reg1  - plot region to 1      reg2  - plot region to 2\n");
    PRINTF("nocopy - toggle noCopy mode for mesh-based tests\n");
    PRINTF("cmode - toggle dumping of colormap in hcp file\n");
    PRINTF("earth - use earthtones palette (default).  Others are:\n");
    PRINTF("     stern   rainbow   heat   gray   yarg\n");
    PRINTF("work  - use work stylesheet (default).  Others are:\n");
    PRINTF("     boxed   axes\n");
    PRINTF("wide      - toggle wide lines in all line drawings\n");
    PRINTF("linetype  - cycle through 5 line types in all line drawings\n");
    PRINTF("closed    - toggle closed lines in plg tests\n");
    PRINTF("rays      - toggle ray arrows in plg tests\n");
    PRINTF("smooth    - cycle through 3 hcp smoothnesses in plg tests\n");
    PRINTF("limits, xmin, xmax, ymin, ymax, nice=1/0, square=1/0,\n");
    PRINTF("        restrict=1/0   - set limits (default limits,e,e,e,e)\n");
    PRINTF("logxy, 1/0, 1/0   - set/reset log scaling for x,y axes\n");
    PRINTF("animate  - toggle animation mode\n");
    PRINTF("marks    - toggle occasional curve markers for line drawings\n");
    PRINTF("legends  - toggle dumping of legends into hcp file\n");
    PRINTF("\n    Movies to test high level Gist functions--\n");
    PRINTF("plg  - test GdLines\n");
    PRINTF("plm  - test GdMesh\n");
    PRINTF("plf  - test GdFillMesh\n");
    PRINTF("plv  - test GdVectors\n");
    PRINTF("plc  - test GdContours\n");
    PRINTF("\n    Single frame tests of high level Gist functions--\n");
    PRINTF("plg1  - test GdLines\n");
    PRINTF("plm1  - test GdMesh\n");
    PRINTF("plf1  - test GdFillMesh\n");
    PRINTF("plc1  - test GdContours\n");
    PRINTF("pli   - test GdCells\n");
    PRINTF("plt   - test GdText\n");
    PRINTF("pls   - test GdLines in polymarker mode\n");
    PRINTF("pldj  - test GdDisjoint\n");
    PRINTF("txin  - test GdText     txout  - 2nd test of GdText\n");
    PRINTF("calib - text calibration frame\n");
    PRINTF("cfont - toggle font for calibration frame\n");
    PRINTF("\n    Tests of high level Gist control functions--\n");
    PRINTF("fma    - frame advance\n");
    PRINTF("hcp    - send current frame to hardcopy file\n");
    PRINTF("hcpon  - send every frame to hardcopy file\n");
    PRINTF("hcpoff - stop sending every frame to hardcopy file\n");
    PRINTF("redraw - redraw the X window\n");
    PRINTF("\nquit  - exit bench program (just q works too)\n");

  } else {
    PRINTF("  NO ACTION TAKEN -- type help for (long) list of commands\n");
  }

  return;
}

static int waitTest= 0;
static void window_exposed(void);
static int no_key = 0;

int
on_idle(void)
{
  if (noKeyboard) {
    if (noKeyTest[no_key]) {
      char line[24];
      strcpy(line, noKeyTest[no_key++]);
      on_stdin(line);
    } else {
      p_quit();
    }
    return 1;
  }
  if (waitTest || !prompt_issued) {
    GhBeforeWait();
    p_stdout("bench> ");
    prompt_issued = 1;
  }
  return 0;
}

static char *default_path = "~/gist:~/Gist:../g";

int
on_launch(int argc, char *argv[])
{
  int i, j;
  int usePS= (argc>1 && argv[1][0]!='0');

  GpReal ypage= 1.033461;  /* 11.0 inches in NDC */
  GpReal xpage= 0.798584;  /* 8.5 inches in NDC */
  GpReal meshSize;

  noKeyboard = waitTest = (argc>1 && argv[1][0]=='-');

  gistPathDefault = default_path;

  p_stdinit(&on_stdin);
  p_idler(&on_idle);
  p_handler(&on_exception);
  p_quitter(&on_quit);

#ifndef NO_XLIB
  g_initializer(&argc, argv);
  engine= GpFXEngine("Gist bench test", 0, DPI, argc>2? argv[2] : 0);

  hcpDefault= usePS? GpPSEngine("Gist bench test", 0, 0, "bench.ps") :
                     GpCGMEngine("Gist bench test", 0, 0, "bench.cgm");
#else
  engine= GpCGMEngine("Gist bench test", 0, 0, "bench.cgm");

  hcpDefault= GpPSEngine("Gist bench test", 0, 0, "bench.ps");
#endif

  ghDevices[0].drawing= GdNewDrawing("work.gs");
  ghDevices[0].display= engine;
  ghDevices[0].hcp= 0;
  ghDevices[0].doLegends= 0;

  actualColors= GpReadPalette(engine, PALETTE_FILE, &palette, 240);
  GpSetPalette(hcpDefault, palette, actualColors);

  /*  GpActivate(engine);  */
  GhSetPlotter(0);

  meshTrans.window= unitBox;

  /* Set up viewport to be 512-by-512 pixels at 100 dpi,
     centered for portrait mode.  */
  meshSize= ONE_INCH*5.125;
  meshTrans.viewport.xmin= 0.5*(xpage-meshSize);
  meshTrans.viewport.xmax= 0.5*(xpage+meshSize);
  meshTrans.viewport.ymin= ypage-0.5*(xpage+meshSize);
  meshTrans.viewport.ymax= ypage-0.5*(xpage-meshSize);

  /* Initialize mesh arrays */
  mesh.iMax= N;
  mesh.jMax= N;
  mesh.x= xmesh= (GpReal *)p_malloc(sizeof(GpReal)*N*N);
  mesh.y= ymesh= (GpReal *)p_malloc(sizeof(GpReal)*N*N);
  mesh.reg= ireg= (int *)p_malloc(sizeof(int)*(N*N+N+1));
  mesh.triangle= 0;
  cmesh= (GpColor *)p_malloc(sizeof(GpColor)*N*N);
  z= (GpReal *)p_malloc(sizeof(GpReal)*N*N);
  u= (GpReal *)p_malloc(sizeof(GpReal)*N*N);
  v= (GpReal *)p_malloc(sizeof(GpReal)*N*N);
  j= 0;
  for (i=0 ; i<N*N+N+1 ; i++) {
    if (j==0 || i<N || i>=N*N) ireg[i]= 0;
    else ireg[i]= 1;
    j++;
    if (j==N) j= 0;
  }

  /* Initialize cosines and sines for perturbation amplitude */
  for (i=0 ; i<N ; i++) cosine[i]= cos(i*PI/(N-1.0));
  for (i=0 ; i<NFRAMES ; i++) amplitude[i]= AMPLITUDE*sin(i*DPHI);

  /* Initialize contour level function */
  PRINTF2("Initializing %dx%d contour function, STANDBY...", N, N);
  for (i=0 ; i<N*N ; i++) z[i]= 0.0;
  TenseZ();
  for (i=0 ; i<300 ; i++) { RelaxZ(); TenseZ(); }
  PRINTF("  DONE\n");

  /* Initialize vector function */
  for (i=0 ; i<N*N ; i++) { u[i]= z[i];  v[N*N-1-i]= z[i]; }

  if (waitTest || noKeyboard) gist_expose_wait(engine, window_exposed);

  return 0;
}

static void
window_exposed(void)
{
  waitTest = 0;
}

int
on_quit(void)
{

  GpDeactivate(hcpDefault);
  GpKillEngine(hcpDefault);
  hcpDefault= 0;

  GpDeactivate(engine);
  GpKillEngine(engine);
  ghDevices[0].display= 0;

  return 0;
}

void
on_exception(int signal, char *errmsg)
{
  PRINTF2("\n\nbench: signal %d, msg = %s\n\n", signal,
          errmsg? errmsg : "<none>");
  p_quit();
}
