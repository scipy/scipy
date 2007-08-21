/*
 * CGMIN.H
 *
 * $Id$
 *
 * Declare the CGM reader/echoer for GIST.
 *
 */
/*    Copyright (c) 1994.  The Regents of the University of California.
                    All rights reserved.  */

#ifndef CGMIN_H
#define CGMIN_H

#include "gist.h"

extern int OpenCGM(char *file);
extern int ReadCGM(int *mPage, int *nPage, int *sPage, int nPageGroups);
extern int CGMRelative(int offset);
extern void CGMinfo(void);
extern int CatalogCGM(void);

extern void Warning(char *general, char *particular);
extern int amBatch, cgmLandscape;
extern Engine *outEngines[8];
extern int outTypes[8];

extern int bg0fg1;

#endif
