/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   FILE: unur_struct.h                                                     *
 *                                                                           *
 *   PURPOSE:                                                                *
 *         declares structures for distribution, parameter, and generator    *
 *         objects.                                                          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   Copyright (c) 2000-2006 Wolfgang Hoermann and Josef Leydold             *
 *   Department of Statistics and Mathematics, WU Wien, Austria              *
 *                                                                           *
 *   This program is free software; you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation; either version 2 of the License, or       *
 *   (at your option) any later version.                                     *
 *                                                                           *
 *   This program is distributed in the hope that it will be useful,         *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of          *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the           *
 *   GNU General Public License for more details.                            *
 *                                                                           *
 *   You should have received a copy of the GNU General Public License       *
 *   along with this program; if not, write to the                           *
 *   Free Software Foundation, Inc.,                                         *
 *   59 Temple Place, Suite 330, Boston, MA 02111-1307, USA                  *
 *                                                                           *
 *****************************************************************************/

/*---------------------------------------------------------------------------*/
#ifndef UNUR_STRUCT_H_SEEN
#define UNUR_STRUCT_H_SEEN
/*---------------------------------------------------------------------------*/

/*****************************************************************************/
/**  Basic header files                                                     **/
/*****************************************************************************/

/*****************************************************************************/
/**  UNU.RAN objects                                                        **/
/*****************************************************************************/

struct unur_distr;    /* distribution object      */
struct unur_par;      /* parameters for generator */
struct unur_gen;      /* generator object         */

/*****************************************************************************/
/**  Generic functions                                                      **/
/*****************************************************************************/

/*---------------------------------------------------------------------------*/
/* Generic functions                                                          */

typedef double UNUR_FUNCT_GENERIC  (double  x, void *params);
typedef double UNUR_FUNCT_VGENERIC (double *x, void *params);

/* for univariate functions with optional parameter array */
struct unur_funct_generic {
  UNUR_FUNCT_GENERIC *f;
  void *params;
};

/* for multivariate functions with optional parameter array */
struct unur_funct_vgeneric {
  UNUR_FUNCT_VGENERIC *f;
  void *params;
};

/*****************************************************************************/
/**  Auxiliary tools                                                        **/
/*****************************************************************************/

#include <utils/slist_struct.h>
#include <utils/string_struct.h>

/*****************************************************************************/
/**  Declaration for parser                                                 **/
/*****************************************************************************/

#include <parser/functparser_struct.h>

/*****************************************************************************/
/**  URNG (uniform random number generator) objects                         **/
/*****************************************************************************/

#include <urng/urng_struct.h>

/*****************************************************************************/
/**  Distribution objects                                                   **/
/*****************************************************************************/

#include <distr/distr_struct.h>

/*****************************************************************************/
/**  Parameter and generators objects                                       **/
/*****************************************************************************/

#include <methods/x_gen_struct.h>

/*---------------------------------------------------------------------------*/
#endif  /* UNUR_STRUCT_H_SEEN */
/*---------------------------------------------------------------------------*/
