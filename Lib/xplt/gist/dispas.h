/*
 * DISPAS.H
 *
 * $Id$
 *
 * Declare dispatcher routines for ordinary file i/o streams
 *
 */
/*    Copyright (c) 1994.  The Regents of the University of California.
                    All rights reserved.  */

#ifndef DISPAS_H
#define DISPAS_H

#include <stdio.h>

/* When input arrives for file, DispatchEvents will call the
   Dispatch method, unless this would cause a SIGTTIN signal.  */
extern int AddFDispatcher(FILE *file,
			  int (*Dispatch)(FILE *file, void *context),
			  void *context);

extern void RemoveFDispatcher(FILE *file);

/* need dispat.h for DispatchEvents declaration */
#include "dispat.h"

#endif
