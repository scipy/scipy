/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   FILE: x_gen_source.h                                                    *
 *                                                                           *
 *   PURPOSE:                                                                *
 *         defines macros and function prototypes for handling               *
 *         generator objects.                                                *
 *                                                                           *
 *   USAGE:                                                                  *
 *         only included in source_unuran.h                                  *
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
/* Invoke generators (macros to avoid function calls)                        */  

#define _unur_init(par)               (par)->init(par)

#define _unur_sample_discr(gen)       (gen)->sample.discr(gen)
#define _unur_sample_cont(gen)        (gen)->sample.cont(gen)
#define _unur_sample_vec(gen,vector)  (gen)->sample.cvec(gen,vector)

#define _unur_free(gen)               do {if(gen) (gen)->destroy(gen);} while(0)

/*---------------------------------------------------------------------------*/
/* get type of transformation method                                         */

#define _unur_gen_is_discr(gen) ( ((gen)->distr->type == UNUR_DISTR_DISCR) ? 1 : 0 )
#define _unur_gen_is_cont(gen)  ( ((gen)->distr->type == UNUR_DISTR_CONT)  ? 1 : 0 )
#define _unur_gen_is_vec(gen)   ( ((gen)->distr->type == UNUR_DISTR_CVEC)  ? 1 : 0 )

/*---------------------------------------------------------------------------*/
/* aux routine when no sampling routine is available                         */

int _unur_sample_discr_error( struct unur_gen *gen );
double _unur_sample_cont_error( struct unur_gen *gen );
int _unur_sample_cvec_error( struct unur_gen *gen, double *vec );
int _unur_sample_matr_error( struct unur_gen *gen, double *mat );

/*---------------------------------------------------------------------------*/
/* create, copy and free parameter object                                    */

/* create an empty parameter object with data structure of size 's'          */ 
struct unur_par *_unur_par_new( size_t s );

struct unur_par *_unur_par_clone( const struct unur_par *par );

/* free memory allocated by parameter obejct                                 */
#define _unur_par_free(par)  do {free((par)->datap); free(par);} while(0)

/*---------------------------------------------------------------------------*/
/* create (new) generic generator object                                     */

struct unur_gen *_unur_generic_create( struct unur_par *par, size_t s );

/*---------------------------------------------------------------------------*/
/* copy (clone) generator objects                                            */

struct unur_gen *_unur_generic_clone( const struct unur_gen *gen, const char *type );

#define _unur_gen_clone(gen)    ((gen)->clone(gen))

/*---------------------------------------------------------------------------*/
/* free generic generator object                                             */

void _unur_generic_free( struct unur_gen *gen );

/*---------------------------------------------------------------------------*/
/* check for type of generator object                                        */

int _unur_gen_is_inversion ( struct unur_gen *gen );

/*---------------------------------------------------------------------------*/
/* set and clone arrays of generator objects                                 */

struct unur_gen **_unur_gen_list_set( struct unur_gen *gen, int n_gen_list );
/* set all entries in list to same generator object                          */
/* IMPORTANT: Be careful when using this call. When the resulting array      */
/*   is stored in some multivariate generator object then 'gen' _must not_   */
/*   be used any more after this call!                                       */

struct unur_gen **_unur_gen_list_clone( struct unur_gen **gen_list, int n_gen_list );
/* clone list of generator objects                                           */

void _unur_gen_list_free( struct unur_gen **gen_list, int n_gen_list );
/* free list of generator objects                                            */

/*---------------------------------------------------------------------------*/
