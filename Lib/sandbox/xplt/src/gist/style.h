/*
 *     STYLE.H
 *
 *         Declare functions and structs in style.c
 *
 */

/*    Copyright (c) 1994.  The Regents of the University of California.
 *                        All rights reserved.  */



#include "draw.h"

typedef struct GfakeSystem GfakeSystem;
struct GfakeSystem {
  double viewport[4];    /* [xmin,xmax,ymin,ymax] in NDC coordinates */
  GaTickStyle ticks;     /* tick style for this coordinate system */
  char *legend;          /* e.g.- "System 0" or "System 1", p_malloc */
};

extern int raw_style(long nsys, int *landscape,
		     GfakeSystem *systems, GeLegendBox *legends);
