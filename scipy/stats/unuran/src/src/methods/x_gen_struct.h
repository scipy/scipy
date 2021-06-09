/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   FILE: x_gen_struct.h                                                    *
 *                                                                           *
 *   PURPOSE:                                                                *
 *         declares structures for parameter and generator objects.          *
 *                                                                           *
 *   USAGE:                                                                  *
 *         only included in unur_struct.h                                    *
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
/* types for sampling routines                                               */

/* for univariate continuous distribution */
typedef double UNUR_SAMPLING_ROUTINE_CONT(struct unur_gen *gen);

/* for univariate discrete distribution */
typedef int UNUR_SAMPLING_ROUTINE_DISCR(struct unur_gen *gen);

/* for multivariate continuous distribution */
typedef int UNUR_SAMPLING_ROUTINE_CVEC(struct unur_gen *gen, double *vec);


/*---------------------------------------------------------------------------*/
/* parameter objects                                                         */

struct unur_par {
  void *datap;                /* pointer to data for method                  */
  size_t s_datap;             /* size of data structure                      */

  struct unur_gen* (*init)(struct unur_par *par);

  unsigned method;            /* indicates method and generator to be used   */
  unsigned variant;           /* indicates variant of method                 */
  unsigned set;               /* stores which parameters have been changed   */

  UNUR_URNG *urng;            /* pointer to uniform random number generator  */
  UNUR_URNG *urng_aux;        /* pointer to second (auxiliary) uniform RNG   */

  const struct unur_distr *distr;  /* pointer to distribution object         */
  int distr_is_privatecopy;   /* whether the distribution object has to be
				 copied into the generator object (TRUE) or
				 just the pointer to the given (external)
				 distribution object (FALSE).

				 Notice: The UNU.RAN design assumes that the
				 generator object keeps its own private copy.
				 However, in some cases it can be useful 
				 to avoid making this copy, e.g. when only
				 a single random variate is required.

				 HOWEVER, this must be used with extreme CARE!

				 When the distrubtion object is changed or
				 freed then the generator object does not work 
				 any more or (even worse) produces garbage.
			      */

  unsigned debug;             /* debugging flags                             */
#ifdef UNUR_COOKIES
  unsigned cookie;            /* magic cookie                                */
#endif
};


/*---------------------------------------------------------------------------*/
/* generator objects                                                         */

struct unur_gen { 
  void *datap;                /* pointer to data for method                  */
  
  union {
    UNUR_SAMPLING_ROUTINE_CONT  *cont;
    UNUR_SAMPLING_ROUTINE_DISCR *discr;
    UNUR_SAMPLING_ROUTINE_CVEC  *cvec;
    UNUR_SAMPLING_ROUTINE_CVEC  *matr;
  } sample;                   /* pointer to sampling routine                 */
  
  UNUR_URNG *urng;            /* pointer to uniform random number generator  */
  UNUR_URNG *urng_aux;        /* pointer to second (auxiliary) uniform RNG   */

  struct unur_distr *distr;   /* distribution object                         */
  int distr_is_privatecopy;   /* whether the distribution object was 
				 copied into the generator object (TRUE) or
				 just the pointer to the given (external)
				 distribution object (FALSE).

				 Notice: The UNU.RAN design assumes that the
				 generator object keeps its own private copy.
				 However, in some cases it can be useful 
				 to avoid making this copy, e.g. when only
				 a single random variate is required.

				 HOWEVER, this must be used with extreme CARE!

				 When the distrubtion object is changed or
				 freed then the generator object does not work 
				 any more or (even worse) produces garbage.
			      */


  unsigned method;            /* indicates method and generator to be used   */
  unsigned variant;           /* indicates variant of method                 */
  unsigned set;               /* stores which parameters have been changed   */
  unsigned status;            /* status of generator object                  */
  
  char *genid;                /* identifier for generator                    */

  struct unur_gen *gen_aux;   /* pointer to auxiliary generator object       */
  struct unur_gen **gen_aux_list; /* list of pointers to auxiliary generator objects */
  int n_gen_aux_list;         /* length of this list                         */

  size_t s_datap;             /* size of data structure                      */
  unsigned debug;             /* debugging flags                             */

  void (*destroy)(struct unur_gen *gen); /* pointer to destructor            */ 
  struct unur_gen* (*clone)(const struct unur_gen *gen ); /* clone generator */
  int (*reinit)(struct unur_gen *gen); /* pointer to reinit routine          */ 

#ifdef UNUR_ENABLE_INFO
  struct unur_string *infostr; /* pointer to info string                     */
  void (*info)(struct unur_gen *gen, int help); /* routine for creating info string */
#endif

#ifdef UNUR_COOKIES
  unsigned cookie;            /* magic cookie                                */
#endif
};

/*---------------------------------------------------------------------------*/
