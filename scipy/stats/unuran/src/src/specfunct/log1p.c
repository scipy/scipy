/* _unur_log1p
 *
 * Compute log(1+x)
 * 
 * Replacement for missing C99 function log1p
 *
 * Copied and renamed from gsl_log1p into _unur_log1p
 * by Josef Leydold, Wed Jun 29 16:39:10 CEST 2005
 */

/* sys/log1p.c
 * 
 * Copyright (C) 1996, 1997, 1998, 1999, 2000 Gerard Jungman, Brian Gough
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or (at
 * your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
 */

#include <config.h>
#include <math.h>
#include "unur_specfunct_source.h"

#if !HAVE_DECL_LOG1P
double _unur_log1p (double x)
{
  volatile double y;
  y = 1 + x;
  return log(y) - ((y-1)-x)/y ;  /* cancels errors with IEEE arithmetic */
}
#endif
