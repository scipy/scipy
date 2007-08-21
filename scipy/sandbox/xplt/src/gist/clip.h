/*
    CLIP.H (ANSI C version)
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

#ifndef CLIP_H
#define CLIP_H

#ifndef GIST_H
#ifndef SINGLE_P
typedef double GpReal;
#else
typedef float GpReal;
#endif
#endif

extern GpReal *xClip, *yClip, *xClip1, *yClip1;

extern void ClipFreeWS(void);

extern void ClipSetup(GpReal xmn, GpReal xmx, GpReal ymn, GpReal ymx);

extern int ClipBegin(const GpReal* xx, const GpReal* yy, int nn, int clsd);

extern int ClipMore(void);

extern int ClipPoints(const GpReal* xx, const GpReal* yy, int nn);

extern int ClipFilled(const GpReal* xx, const GpReal* yy, int nn);

extern int ClipDisjoint(const GpReal* x0, const GpReal* y0,
                        const GpReal* x1, const GpReal* y1, int nn);

extern int ClipTest(const GpReal* xx, const GpReal* yy, int nn, int clsd,
                    const GpReal* box);

#endif
