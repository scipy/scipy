/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   FILE: unur_typedefs.h                                                   *
 *                                                                           *
 *   PURPOSE:                                                                *
 *         typedefs and names of globally used structures                    *
 *         (not already defined in unuran_config.h)                          *
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
#ifndef UNUR_TYPEDEFS_H_SEEN
#define UNUR_TYPEDEFS_H_SEEN
/*---------------------------------------------------------------------------*/
/* UNURAN objects                                                            */

struct unur_distr;                       /* distribution object              */
typedef struct unur_distr UNUR_DISTR;

struct unur_par;                         /* parameters for generator         */
typedef struct unur_par   UNUR_PAR;

struct unur_gen;                         /* generator object                 */
typedef struct unur_gen   UNUR_GEN;

struct unur_urng;                        /* uniform random number generator  */
typedef struct unur_urng  UNUR_URNG;
/* comment out the following definition when type 'UNUR_URNG' is changed !!  */
/* You have to define  _unur_call_urng(urng) and _unur_call_reset(urng)      */
/* in file 'src/urng/urng_source.h' then!                                    */
#define UNUR_URNG_UNURAN 1

/*---------------------------------------------------------------------------*/
/* functions for continuous univariate PDF, CDF, and their derivatives       */

typedef double UNUR_FUNCT_CONT  (double x, const struct unur_distr *distr);
typedef double UNUR_FUNCT_DISCR (int x, const struct unur_distr *distr);
typedef int    UNUR_IFUNCT_DISCR(double x, const struct unur_distr *distr);

/*---------------------------------------------------------------------------*/
/* functions for continuous multivariate PDF, CDF, and their gradients       */

/* Remark: we cannot use a const pointer to distr as some data are           */
/* computed "on the fly" when they are needed and stored in the              */
/* distribution object.                                                      */

typedef double UNUR_FUNCT_CVEC (const double *x, struct unur_distr *distr);
typedef int    UNUR_VFUNCT_CVEC(double *result, const double *x, struct unur_distr *distr);
typedef double UNUR_FUNCTD_CVEC(const double *x, int coord, struct unur_distr *distr);

/*---------------------------------------------------------------------------*/
/* structures for auxiliary tools                                            */

struct unur_slist;         /* structure for simple list                      */

/*---------------------------------------------------------------------------*/
/* error handler                                                             */

typedef void UNUR_ERROR_HANDLER( const char *objid, const char *file, int line, 
				 const char *errortype, int unur_errno, const char *reason );

/*---------------------------------------------------------------------------*/
#endif  /* UNUR_TYPEDEFS_H_SEEN */
/*---------------------------------------------------------------------------*/
