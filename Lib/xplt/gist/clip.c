/*
    CLIP.C
    Routines to perform floating point clipping operations.

    $Id$

    Interface routines:

    *** include "clip.h" to define xClip, yClip ***
    int nc;
    ClipSetup(xmin, xmax, ymin, ymax);

       *** to draw an open curve, use ***
    if (ClipBegin(x, y, n, 0)) DrawPolyline(x, y, n);
    else while (nc=ClipMore()) DrawPolyline(xClip, yClip, nc);

       *** to draw a closed curve, use ***
    if (ClipBegin(x, y, n, 1)) DrawPolygon(x, y, n);
    else while (nc=ClipMore()) DrawPolyline(xClip, yClip, nc);

       *** to draw a set of unconnected points, use ***
    if (nc=ClipPoints(x, y, n)) DrawPolymarker(xClip, yClip, nc);

       *** to draw a closed filled curve, use ***
    if (nc=ClipFilled(x, y, n)) DrawFilledPolygon(xClip, yClip, nc);

       *** to draw a disjoint segments, use ***
    if (nc=ClipDisjoint(x0, y0, x1, y1, n))
       DrawDisjoint(xClip, yClip, xClip1, yClip1, nc);
 */
/*    Copyright (c) 1994.  The Regents of the University of California.
                    All rights reserved.  */

#include "clip.h"

extern void (*GmFree)(void *);
extern void *(*GmMalloc)(long);

/*-----------------------------------------------------------------------*/
/* State variables */

/* Global variables can be set to (x,y) or (xwork,ywork) as needed. */
GpReal *xClip, *yClip;
GpReal *xClip1, *yClip1;  /* for disjoint segments */

/* xmin, xmax, ymin, ymax -- clip bounds (set by ClipSetup)
   x, y -- list of points to be clipped
   n -- number of (x,y)
   i -- current index into (x,y) (state info for ClipMore)
   xwork, ywork -- workspace (ClipBegin, ClipFilled, ClipPoints, ClipFreeWS)
   nwork, maxwork -- current and maximum number of points in WS
   side -- 0 means out-of-bounds left   (state info for ClipMore)
           1 means out-of-bounds bottom
	   2 means out-of-bounds right
	   3 means out-of bounds top
   wind -- winding number for out-of-bounds excursions.  Increments for
           each ccw corner crossing, decrements for each cw crossing.
	   When combined with side, this allows corners to be inserted
	   in ClipFilled algorithm.
   xclose, yclose -- space for closure segment of closed curve
   xc, yc -- space for the four corners for ClipFilled
 */
static const GpReal *x, *y;
static GpReal xmin, xmax, ymin, ymax;
static int n, i, closed, nclosed;
static GpReal *xwork, *ywork;
static int nwork, maxwork=0;
static int side, wind;
static GpReal xclose[2], yclose[2], xc[4], yc[4];

static int FirstScan(void);
static int FindEntry(GpReal* xx, GpReal* yy);
static int FindExit(void);
static void Copy1stN(GpReal* xx, GpReal* yy, int nn);
static void SetupClosure(int ex);
static void WindCorners(int n);

/*-----------------------------------------------------------------------*/
/* Worker routines.
   FindEntry and FindExit are the most important routines.  The hope is
   that most of the time will be spent in the FindExit loop and the
   four loops in FindEntry.  The rest of the code is not too carefully
   optimized. */

static int FirstScan(void)
{
  /* Scan (x,y) until outside xmin,xmax,ymin,ymax, returning 1 if never
     outside bounds, 0 if scan stopped.  Sets side to value
     at point (x[i],y[i]) if 0 returned; this is first point outside
     clip bounds. */
  for (i=0 ; i<n ; i++) {
    if (x[i]<xmin) {side=0; return 0;}
    if (y[i]<ymin) {side=1; return 0;}
    if (x[i]>xmax) {side=2; return 0;}
    if (y[i]>ymax) {side=3; return 0;}
  }
  return 1;
}

static int FindEntry(GpReal* xx, GpReal* yy)
{
  /* Given exit side flag, scan from (x[i],y[i]) to
     (x[n],y[n]) for first entry into xmin,xmax,ymin,ymax.  Return 0
     if no such point, else return 1 and set (xx,yy) to entry point.
     On exit, if return 0, side set to continue FindEntry after SetupClosure,
     but if return 1, side is set as on entry.
     The winding number, wind, is incremented for counterclockwise
     corner crossings, decremented for clockwise corner crossings. */

  if (side==0) goto left;
  if (side==1) goto bottom;
  if (side==2) goto right;
  if (side==3) goto top;

  return 0; /* no exit point */

  /* wind algorithm: Assume FindEntry is called when you have
     just emerged from the initial side.  Before each transition
     to another side (or entry), increment or decrement wind to
     refect how many corners you would wrap if you entered on
     the new side.
     The assumption that you just emerged from the initial side
     is not true for the first point -- the side set by FirstScan.  */

 left:
  /* Get back in bounds if possible. */
  for (i++ ; i<n ; i++) if (x[i]>=xmin) break;
  if (i>=n) {side=0; return 0;}

  /* Clip to find entry point. */
  if (y[i-1]<ymin) {
    if (y[i]<ymin) {wind++; goto bottom;}
    *yy= (xmin-x[i])*(y[i]-y[i-1])/(x[i]-x[i-1])+y[i];
    if (*yy<ymin) {
      wind++;
      *xx= (ymin-y[i])*(x[i]-x[i-1])/(y[i]-y[i-1])+x[i];
      if (*xx>xmax) {wind++; goto right;}
      *yy= ymin;
      return 1;
    }
    if (*yy>ymax) {wind--; goto top;}
    *xx= xmin;
    return 1;
  } else if (y[i-1]>ymax) {
    if (y[i]>ymax) {wind--; goto top;}
    *yy= (xmin-x[i])*(y[i]-y[i-1])/(x[i]-x[i-1])+y[i];
    if (*yy>ymax) {
      wind--;
      *xx= (ymax-y[i])*(x[i]-x[i-1])/(y[i]-y[i-1])+x[i];
      if (*xx>xmax) {wind--; goto right;}
      *yy= ymax;
      return 1;
    }
    if (*yy<ymin) {wind++; goto bottom;}
    *xx= xmin;
    return 1;
  } else {
    *yy= (xmin-x[i])*(y[i]-y[i-1])/(x[i]-x[i-1])+y[i];
    if (*yy<ymin) {wind++; goto bottom;}
    if (*yy>ymax) {wind--; goto top;}
    *xx= xmin;
    return 1;
  }

 bottom:
  /* Get back in bounds if possible. */
  for (i++ ; i<n ; i++) if (y[i]>=ymin) break;
  if (i>=n) {side=1; return 0;}

  /* Clip to find entry point. */
  if (x[i-1]<xmin) {
    if (x[i]<xmin) {wind--; goto left;}
    *xx= (ymin-y[i])*(x[i]-x[i-1])/(y[i]-y[i-1])+x[i];
    if (*xx<xmin) {
      wind--;
      *yy= (xmin-x[i])*(y[i]-y[i-1])/(x[i]-x[i-1])+y[i];
      if (*yy>ymax) {wind--; goto top;}
      *xx= xmin;
      return 1;
    }
    if (*xx>xmax) {wind++; goto right;}
    *yy= ymin;
    return 1;
  } else if (x[i-1]>xmax) {
    if (x[i]>xmax) {wind++; goto right;}
    *xx= (ymin-y[i])*(x[i]-x[i-1])/(y[i]-y[i-1])+x[i];
    if (*xx>xmax) {
      wind++;
      *yy= (xmax-x[i])*(y[i]-y[i-1])/(x[i]-x[i-1])+y[i];
      if (*yy>ymax) {wind++; goto top;}
      *xx= xmax;
      return 1;
    }
    if (*xx<xmin) {wind--; goto left;}
    *yy= ymin;
    return 1;
  } else {
    *xx= (ymin-y[i])*(x[i]-x[i-1])/(y[i]-y[i-1])+x[i];
    if (*xx<xmin) {wind--; goto left;}
    if (*xx>xmax) {wind++; goto right;}
    *yy= ymin;
    return 1;
  }

 right:
  /* Get back in bounds if possible. */
  for (i++ ; i<n ; i++) if (x[i]<=xmax) break;
  if (i>=n) {side=2; return 0;}

  /* Clip to find entry point. */
  if (y[i-1]<ymin) {
    if (y[i]<ymin) {wind--; goto bottom;}
    *yy= (xmax-x[i])*(y[i]-y[i-1])/(x[i]-x[i-1])+y[i];
    if (*yy<ymin) {
      wind--;
      *xx= (ymin-y[i])*(x[i]-x[i-1])/(y[i]-y[i-1])+x[i];
      if (*xx<xmin) {wind--; goto left;}
      *yy= ymin;
      return 1;
    }
    if (*yy>ymax) {wind++; goto top;}
    *xx= xmax;
    return 1;
  } else if (y[i-1]>ymax) {
    if (y[i]>ymax) {wind++; goto top;}
    *yy= (xmax-x[i])*(y[i]-y[i-1])/(x[i]-x[i-1])+y[i];
    if (*yy>ymax) {
      wind++;
      *xx= (ymax-y[i])*(x[i]-x[i-1])/(y[i]-y[i-1])+x[i];
      if (*xx<xmin) {wind++; goto left;}
      *yy= ymax;
      return 1;
    }
    if (*yy<ymin) {wind--; goto bottom;}
    *xx= xmax;
    return 1;
  } else {
    *yy= (xmax-x[i])*(y[i]-y[i-1])/(x[i]-x[i-1])+y[i];
    if (*yy<ymin) {wind--; goto bottom;}
    if (*yy>ymax) {wind++; goto top;}
    *xx= xmax;
    return 1;
  }

 top:
  /* Get back in bounds if possible. */
  for (i++ ; i<n ; i++) if (y[i]<=ymax) break;
  if (i>=n) {side=3; return 0;}

  /* Clip to find entry point. */
  if (x[i-1]<xmin) {
    if (x[i]<xmin) {wind++; goto left;}
    *xx= (ymax-y[i])*(x[i]-x[i-1])/(y[i]-y[i-1])+x[i];
    if (*xx<xmin) {
      wind++;
      *yy= (xmin-x[i])*(y[i]-y[i-1])/(x[i]-x[i-1])+y[i];
      if (*yy<ymin) {wind++; goto bottom;}
      *xx= xmin;
      return 1;
    }
    if (*xx>xmax) {wind--; goto right;}
    *yy= ymax;
    return 1;
  } else if (x[i-1]>xmax) {
    if (x[i]>xmax) {wind--; goto right;}
    *xx= (ymax-y[i])*(x[i]-x[i-1])/(y[i]-y[i-1])+x[i];
    if (*xx>xmax) {
      wind--;
      *yy= (xmax-x[i])*(y[i]-y[i-1])/(x[i]-x[i-1])+y[i];
      if (*yy<ymin) {wind--; goto bottom;}
      *xx= xmax;
      return 1;
    }
    if (*xx<xmin) {wind--; goto left;}
    *yy= ymax;
    return 1;
  } else {
    *xx= (ymax-y[i])*(x[i]-x[i-1])/(y[i]-y[i-1])+x[i];
    if (*xx<xmin) {wind++; goto left;}
    if (*xx>xmax) {wind--; goto right;}
    *yy= ymax;
    return 1;
  }
}

static int FindExit(void)
{
  /* Return 1 if found an exit, 0 if none.  In either case, on exit
     there will be nwork in-bounds points in (xwork,ywork).  If found,
     must leave side set properly for following FindEntry. */
  int offL, offB, offR, offT;

  for ( ; i<n ; i++) {
    offL= x[i]<xmin;
    offB= y[i]<ymin;
    offR= !offL && x[i]>xmax;
    offT= !offB && y[i]>ymax;
    if (offL || offB || offR || offT) break;
    xwork[nwork]= x[i];
    ywork[nwork]= y[i];
    nwork++;
  }
  if (i>=n) return 0;

  if (offL) {
    ywork[nwork]= (xmin-x[i])*(y[i]-y[i-1])/(x[i]-x[i-1])+y[i];
    if (offB && ywork[nwork]<ymin) {
      xwork[nwork]= (ymin-y[i])*(x[i]-x[i-1])/(y[i]-y[i-1])+x[i];
      ywork[nwork++]= ymin;
      side= 1;
      return 1;
    } else if (offT && ywork[nwork]>ymax) {
      xwork[nwork]= (ymax-y[i])*(x[i]-x[i-1])/(y[i]-y[i-1])+x[i];
      ywork[nwork++]= ymax;
      side= 3;
      return 1;
    } else {
      xwork[nwork++]= xmin;
      side= 0;
      return 1;
    }
  } else if (offB) {
    xwork[nwork]= (ymin-y[i])*(x[i]-x[i-1])/(y[i]-y[i-1])+x[i];
    if (offR && xwork[nwork]>xmax) {
      ywork[nwork]= (xmax-x[i])*(y[i]-y[i-1])/(x[i]-x[i-1])+y[i];
      xwork[nwork++]= xmax;
      side= 2;
      return 1;
    } else {
      ywork[nwork++]= ymin;
      side= 1;
      return 1;
    }
  } else if (offR) {
    ywork[nwork]= (xmax-x[i])*(y[i]-y[i-1])/(x[i]-x[i-1])+y[i];
    if (offT && ywork[nwork]>ymax) {
      xwork[nwork]= (ymax-y[i])*(x[i]-x[i-1])/(y[i]-y[i-1])+x[i];
      ywork[nwork++]= ymax;
      side= 3;
      return 1;
    } else {
      xwork[nwork++]= xmax;
      side= 2;
      return 1;
    }
  } else /* if (offT) */ {
    xwork[nwork]= (ymax-y[i])*(x[i]-x[i-1])/(y[i]-y[i-1])+x[i];
    ywork[nwork++]= ymax;
    side= 3;
    return 1;
  }
}

static void Copy1stN(GpReal* xx, GpReal* yy, int nn)
{
  for (i=0 ; i<nn ; i++) {
    xx[i]= x[i];
    yy[i]= y[i];
  }
}

static void SetupClosure(int ex)
{
  xclose[0]= x[n-1];
  yclose[0]= y[n-1];
  xclose[1]= x[0];
  yclose[1]= y[0];
  x= xclose;
  y= yclose;
  n= 2;
  i= ex;   /* 0 if FindEntry to be called, 1 if FindExit to be called */
  closed= 0;
}

static void WindCorners(int n)
{
  if (n) {
    if (n>0) {
      while (n--) {
	if (side>3) side= 0;
	xwork[nwork]= xc[side];
	ywork[nwork++]= yc[side++];
      }
    } else {
      while (n++) {
	if (--side<0) side= 3;
	xwork[nwork]= xc[side];
	ywork[nwork++]= yc[side];
      }
    }
  }
}

/*-----------------------------------------------------------------------*/
/* Interface routines */

void ClipFreeWS(void)
{
  if (maxwork) {
    GmFree(xwork); GmFree(ywork);
    maxwork= 0;
  }
}

void ClipSetup(GpReal xmn, GpReal xmx, GpReal ymn, GpReal ymx)
{
  if (xmn<=xmx) {
    xmin= xmn;
    xmax= xmx;
  } else {
    xmin= xmx;
    xmax= xmn;
  }
  if (ymn<=ymx) {
    ymin= ymn;
    ymax= ymx;
  } else {
    ymin= ymx;
    ymax= ymn;
  }
}

int ClipBegin(const GpReal* xx, const GpReal* yy, int nn, int clsd)
{
  x= xx;
  y= yy;
  n= nn>1? nn:0;

  /* Make quick check to see if all points are in bounds.  This backfires
     if some points aren't, since they will have to be fetched a second
     time to be copied into the workspace.  However, it gives best possible
     performance if the points are all in bounds. */
  if (FirstScan()) return 1;
  closed= clsd;
  wind= 0;

  /* Get workspace if necessary.
     Worst case requires n points for open curves, n+1 points for closed. */
  if (n+1>maxwork) {
    ClipFreeWS();
    maxwork= n+256;
    xwork= (GpReal*)GmMalloc(sizeof(GpReal)*maxwork);
    ywork= (GpReal*)GmMalloc(sizeof(GpReal)*maxwork);
  }
  xClip= xwork;
  yClip= ywork;

  if (!closed) {
    /* For open curves, copy initial unclipped points into workspace. */
    Copy1stN(xwork, ywork, i);
    nwork= i;
    /* FindExit() will have to discover (x[i],y[i]) is out-of-bounds,
       but this and the second fetch to do the preceding copy are the
       only duplications of effort. */
  } else {
    /* For closed curves, initial unclipped points must be connected to
       final points.  Therefore, remember how many there are. */
    nwork= 0;
    nclosed= i; /* ie- index of 1st out-of-bounds point */
  }
  return 0;
}

int ClipMore(void)
{
  int nw;
  if (i>=n) return 0;

  /* nwork>0 if and only if this is the first pass through an open
     curve whose first nwork points are unclipped. */
  if (nwork==0) {
    if (FindEntry(xwork,ywork)) nwork= 1;
    else if (!closed) return 0;
    else {
      /* There are 4 possiblities for the end of a closed curve:
         point 0 in-bounds or out-of-bounds (nclosed 0 or not)
	 and point n-1 in-bounds or out-of-bounds (FindExit or
	 FindEntry eventually fails).  Since FindEntry has failed,
	 these two cases both have point n-1 out-of-bounds. */
      if (nclosed) {
	/* If nclosed>0, we know the final segment will re-enter, and
	   the final batch of points can be inserted, leaving space
	   for the first point. */
	int side0= side; /* save current FindEntry side */
	Copy1stN(xwork+1, ywork+1, nclosed); /* sets i=nclosed */
	nwork= nclosed+1;
	FindExit();      /* exits immediately, incrementing nwork */
	SetupClosure(0);
	side= side0;     /* restore FindEntry side */
	FindEntry(xwork, ywork);

      } else {
	/* If the first point is out-of-bounds, the final segment may
	   still re-enter, though it must immediately exit. */
	SetupClosure(0);
	if (!FindEntry(xwork, ywork)) return 0;
	nwork= 1;
	FindExit();
      }
      i= n;              /* ensure prompt exit next pass */
      return nwork;
    }
  }

  if (!FindExit() && closed) {
    /* Continue with other two closed endcases.  Since FindExit has
       failed, these two both have point n-1 in-bounds. */
    if (nclosed) {
      /* first point was also inbounds */
      Copy1stN(xwork+nwork, ywork+nwork, nclosed); /* sets i=nclosed */
      nwork+= nclosed;
    } else {
      /* first point was out-of-bounds */
      SetupClosure(1);
    }
    FindExit(); /* Clip final segment. */
    i= n;              /* ensure prompt exit next pass */
    return nwork;
  }

  nw= nwork;
  nwork= 0;
  return nw; /* must process this batch of points */
}

int ClipPoints(const GpReal* xx, const GpReal* yy, int nn)
{
  x= xx;
  y= yy;
  n= nn>0? nn:0;

  /* Make quick check to see if all points are in bounds.  This backfires
     if some points aren't, since they will have to be fetched a second
     time to be copied into the workspace.  However, it gives best possible
     performance if the points are all in bounds. */
  if (FirstScan()) {
    xClip= (GpReal *)x;
    yClip= (GpReal *)y;
    return nn;
  }

  /* Get workspace if necessary.
     Worst case requires n-1 points (at least one clipped). */
  if (n-1>maxwork) {
    ClipFreeWS();
    maxwork= n+256;
    xwork= (GpReal*)GmMalloc(sizeof(GpReal)*maxwork);
    ywork= (GpReal*)GmMalloc(sizeof(GpReal)*maxwork);
  }
  xClip= xwork;
  yClip= ywork;

  i= nwork= 0;
  if (y[i]<ymin) goto bottom;
  if (y[i]>ymax) goto top;
  if (x[i]>xmax) goto right;
  if (x[i]>=xmin) goto inbounds;

  left:
    for (i++ ; i<n ; i++) if (x[i]>=xmin) break;
    if (i>=n) return nwork;
    if (x[i]>xmax) goto right;
    if (y[i]>ymax) goto top;
    if (y[i]>=ymin) goto inbounds;

  bottom:
    for (i++ ; i<n ; i++) if (y[i]>=ymin) break;
    if (i>=n) return nwork;
    if (y[i]>ymax) goto top;
    if (x[i]<xmin) goto left;
    if (x[i]<=xmax) goto inbounds;

  right:
    for (i++ ; i<n ; i++) if (x[i]<=xmax) break;
    if (i>=n) return nwork;
    if (x[i]<xmin) goto left;
    if (y[i]<ymin) goto bottom;
    if (y[i]<=ymax) goto inbounds;

  top:
    for (i++ ; i<n ; i++) if (y[i]<=ymax) break;
    if (i>=n) return nwork;
    if (y[i]<ymin) goto bottom;
    if (x[i]<xmin) goto left;
    if (x[i]>xmax) goto right;

  inbounds:
    xwork[nwork]= x[i];
    ywork[nwork]= y[i];
    nwork++;
    for (i++ ; i<n ; i++) {
      if (x[i]<xmin) goto left;
      if (y[i]<ymin) goto bottom;
      if (x[i]>xmax) goto right;
      if (y[i]>ymax) goto top;
      xwork[nwork]= x[i];
      ywork[nwork]= y[i];
      nwork++;
    }

  return nwork;
}

int ClipFilled(const GpReal* xx, const GpReal* yy, int nn)
{
  GpReal xe, ye;
  int wind0, side0, sideWind, reentry;

  x= xx;
  y= yy;
  n= nn>1? nn:0;

  /* Make quick check to see if all points are in bounds.  This backfires
     if some points aren't, since they will have to be fetched a second
     time to be copied into the workspace.  However, it gives best possible
     performance if the points are all in bounds. */
  if (FirstScan()) {
    xClip= (GpReal *)x;
    yClip= (GpReal *)y;
    return n;
  }

  /* Get workspace if necessary.
     Proceding from one point to the next can cause at most 3 points
     to be placed in the workspace: 1 corner point from the direction
     change plus 2 clips as the new edge cuts through the clip
     rectangle.  The worst case figures therefore look like multipointed
     stars with all of their vertices outside the boundaries and each
     edge cutting through the clip rectangle.  I don't think that a
     full 3*n points is possible, but I can't put a stricter upper
     bound on the worst possible case (e.g.- 11 is worst case for
     n=4).  For arbitrary n, I don't believe the worst case departs
     from 3*n by any significant factor. */
  if (3*n>maxwork) {
    ClipFreeWS();
    maxwork= 3*n+256;
    xwork= (GpReal*)GmMalloc(sizeof(GpReal)*maxwork);
    ywork= (GpReal*)GmMalloc(sizeof(GpReal)*maxwork);
  }
  xClip= xwork;
  yClip= ywork;

  /* Store corners in case out-of-bounds points wind around. */
  xc[0]= xmin;  yc[0]= ymin;
  xc[1]= xmax;  yc[1]= ymin;
  xc[2]= xmax;  yc[2]= ymax;
  xc[3]= xmin;  yc[3]= ymax;

  wind= 0;
  nclosed= i;
  if (nclosed) {
    /* Copy initial unclipped points into workspace. */
    Copy1stN(xwork, ywork, nclosed); /* sets i=nclosed */
    nwork= nclosed;
  } else {
    /* Find first entry point. */
    sideWind= side;
    if (!FindEntry(xwork, ywork)) {
      if (!wind) return 0;
      side0= sideWind;
      wind0= 0;
      nwork= 0;
      goto chkcls;
    }
    nwork= 1;
    wind0= wind; /* save for closure calculation */
    side0= side;
  }

  for (;;) {
    /* Find the next exit point. */
    if (!FindExit()) {
      /* Do two of the possible 4 endcases.  Since FindExit has
	 failed, these two both have point n-1 in-bounds. */
      if (nclosed) {
	/* first point was also inbounds */
	xwork[nwork]= x[0];
	ywork[nwork++]= y[0];
      } else {
	/* first point was out-of-bounds */
	SetupClosure(1);
	FindExit(); /* Clip final segment. */
	if (side!=side0) {
	  /* FirstScan began on a different side from the final exit.
	     The difference can be at most one side -- adjust for it.  */
	  if (((side+1)&3)==side0) wind++;
	  else wind--;
	}
	if (wind) WindCorners(wind);
      }
      break;
    }

    /* Find the next entry point. */
    sideWind= side;
    wind= 0;
    if (!FindEntry(&xe, &ye)) {
      /* Do other 2 endcases.  Since FindEntry has failed,
	 these two cases both have point n-1 out-of-bounds, and
	 are distinguished by whether point 0 is out-of-bounds. */
    chkcls:
      SetupClosure(0); /* preserves wind from previous FindEntry */
      reentry= FindEntry(&xe, &ye);
      if (!reentry) {
	if (side!=side0) {
	  /* FirstScan began on a different side from the final exit.
	     The difference can be at most one side -- adjust for it.  */
	  if (((side+1)&3)==side0) wind0++;
	  else wind0--;
	}
	wind+= wind0;
	if (wind) { side= sideWind;    WindCorners(wind); }
	break;
      }
      if (wind) { side= sideWind;    WindCorners(wind); }
      xwork[nwork]= xe;
      ywork[nwork++]= ye;
      if (nclosed) break; /* If point 0 inbounds, nothing more to do. */
      /* Hardest case is both n-1 and 0 points out-of-bounds, but
	 the closure segment connecting them reenters. */
      FindExit(); /* guaranteed immediate exit */
      if (side!=side0) {
	/* FirstScan began on a different side from the final exit.
	   The difference can be at most one side -- adjust for it.  */
	if (((side+1)&3)==side0) wind0++;
	else wind0--;
      }
      if (wind0) WindCorners(wind0);
      break;
    }

    /* Check for corner crossings. */
    if (wind) { side= sideWind;    WindCorners(wind); }

    xwork[nwork]= xe;
    ywork[nwork++]= ye;
  }

  return nwork;
}

int ClipDisjoint(const GpReal* x0, const GpReal* y0,
		 const GpReal* x1, const GpReal* y1, int nn)
{
  /* This routine is quite far from optimal...  */
  GpReal xx[2], yy[2];
  int notClipped, seg, nsegs= 0;

  if (2*nn+2>maxwork) {
    ClipFreeWS();
    maxwork= 2*nn+256;
    xwork= (GpReal*)GmMalloc(sizeof(GpReal)*maxwork);
    ywork= (GpReal*)GmMalloc(sizeof(GpReal)*maxwork);
  }

  for (seg=0 ; seg<nn ; seg++) {
    xx[0]= x0[seg];   yy[0]= y0[seg];
    xx[1]= x1[seg];   yy[1]= y1[seg];
    if ((notClipped= ClipBegin(xx, yy, 2, 0))) {
      xClip= xx;
      yClip= yy;
    }
    if (notClipped || ClipMore()) {
      xwork[2   +nsegs]= xClip[0];  ywork[2   +nsegs]= yClip[0];
      xwork[2+nn+nsegs]= xClip[1];  ywork[2+nn+nsegs]= yClip[1];
      nsegs++;
    }
  }

  xClip= xwork+2;
  yClip= ywork+2;
  xClip1= xwork+nn+2;
  yClip1= ywork+nn+2;

  return nsegs;
}

int ClipTest(const GpReal* xx, const GpReal* yy, int nn, int clsd,
	     const GpReal* box)
{
  GpReal xns=xmin, yns=ymin, xxs=xmax, yxs=ymax;
  GpReal xe, ye;
  xmin= box[0];  ymin= box[1];  xmax= box[2];  ymax= box[3];
  x= xx;  y= yy;  n= nn;

  if (!FirstScan()) {
    if (!FindEntry(&xe, &ye)) {
      if (clsd) {
	SetupClosure(0);
	if (!FindEntry(&xe, &ye)) i= nn+1;
	else i= nn;
      }
    }
  }

  xmin= xns;  ymin= yns;  xmax= xxs;  ymax= yxs;
  return i;
}
