/*
   STYLE.C
   Set/get details of graphics style for the get_style, set_style  functions.

   $Id$
 */
/*    Copyright (c) 1996.  The Regents of the University of California.
                    All rights reserved.  */

#include "hlevel.h"
#include "pstdlib.h"
#include "style.h"

/* This function is based on raw_style in yorick.
   If nsys==0, this is a query operation which fills in the input
               data arrays (if the pointers are non-zero)
   otherwise, systems is systems[nsys], and the operation is to
               set the values specified in the data arrays
   Unless an error occurs, the return value is the number of coordinate
   systems; for queries the routine must be called twice, first with
   systems==0 to retrieve the number of systems, then with a large enough
   systems to hold the returned values. If an error occurs, the function
   returns -1 instead of the number of coordinate systems. */
extern int raw_style(long nsys, int *landscape,
		     GfakeSystem *systems, GeLegendBox *legends);

int raw_style(long nsys, int *landscape,
	      GfakeSystem *systems, GeLegendBox *legends)
{
  int nsy= GhGetPlotter();
  Drauing *drawing= (nsy>=0 && nsy<8)? ghDevices[nsy].drawing : 0;
  GeSystem *sys= drawing? drawing->systems : 0;

  if (!nsys) {
    /* query operation */
    if (!drawing) return 0;
    if (landscape) *landscape= drawing->landscape;
    nsy= drawing->nSystems;
    if (systems && nsy>0) {
      int i;
      for (i=0 ; i<nsy ; i++,sys=(GeSystem *)sys->el.next) {
	if (systems[i].legend) {
	  char *legend= systems[i].legend;
	  systems[i].legend= 0;
	  p_free(legend);
	}
	systems[i].legend= p_strcpy(sys->el.legend);
	systems[i].viewport[0]= sys->trans.viewport.xmin;
	systems[i].viewport[1]= sys->trans.viewport.xmax;
	systems[i].viewport[2]= sys->trans.viewport.ymin;
	systems[i].viewport[3]= sys->trans.viewport.ymax;
	/* lazy today -- use ANSI struct assignment */
	systems[i].ticks= sys->ticks;
      }
    }
    if (legends) {
      /* lazy today -- use ANSI struct assignment */
      legends[0]= drawing->legends[0];
      legends[1]= drawing->legends[1];
    }
    return nsy;

  } else {
    /* set new style operation */
    int i;
    extern void GdKillSystems(void);  /* defined in draw.c */
    GpBox vp;

    if (!landscape || !systems || !legends) return -1;

    /* don't clobber the current display list(s) unless the
       number of coordinate systems has changed */
    nsy= drawing? drawing->nSystems : 0;
    if (nsy != nsys) GdKillSystems();

    for (i=0 ; i<nsys ; i++) {
      gistD.hidden= 0;
      gistD.legend= systems[i].legend;
      vp.xmin= systems[i].viewport[0];
      vp.xmax= systems[i].viewport[1];
      vp.ymin= systems[i].viewport[2];
      vp.ymax= systems[i].viewport[3];
      if (nsy==nsys) {
	GdSetSystem(i+1);
	/* don't bother with the legend */
	gistD.trans.viewport.xmin= vp.xmin;
	gistD.trans.viewport.xmax= vp.xmax;
	gistD.trans.viewport.ymin= vp.ymin;
	gistD.trans.viewport.ymax= vp.ymax;
	/* lazy today -- use ANSI struct assignment */
	gistD.ticks= systems[i].ticks;
	GdSetPort();
      } else if (GdNewSystem(&vp, &systems[i].ticks)<0) {
	gistD.legend= 0;
	return -1;
      }
      gistD.legend= 0;
    }
    if (nsys && nsy==nsys) GdSetSystem(1);

    for (i=0 ; i<2 ; i++)
      GdLegendBox(i, legends[i].x, legends[i].y,
		  legends[i].dx, legends[i].dy,
		  &legends[i].textStyle, legends[i].nchars,
		  legends[i].nlines, legends[i].nwrap);

    GdLandscape(*landscape);

    return (int)nsys;
  }
}

