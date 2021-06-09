/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   FILE: slist.h                                                           *
 *                                                                           *
 *   PURPOSE:                                                                *
 *         defines function prototypes for simple list                       *
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
#ifndef SLIST_H_SEEN
#define SLIST_H_SEEN
/*---------------------------------------------------------------------------*/
/* Not part of manual!                                                       */
/*---------------------------------------------------------------------------*/
/*                                                                           */
/* A simple list can be used to store an arbitrary numbers of pointers       */
/* to allocated memory in a list.                                            */
/*                                                                           */
/* IMPORTANT: These elements must be allocated via (c|m|re)alloc()!!         */
/*                                                                           */
/*---------------------------------------------------------------------------*/

struct unur_slist *_unur_slist_new( void );
/*---------------------------------------------------------------------------*/
/* Make new simple list.                                                     */
/*---------------------------------------------------------------------------*/

int _unur_slist_append( struct unur_slist *slist, void *element );
/*---------------------------------------------------------------------------*/
/* Append pointer to element to simple list.                                 */
/*---------------------------------------------------------------------------*/

int _unur_slist_length( const struct unur_slist *slist );
/*---------------------------------------------------------------------------*/
/* Get length if list (number of list entries).                              */
/*---------------------------------------------------------------------------*/

void *_unur_slist_get( const struct unur_slist *slist, int n );
/*---------------------------------------------------------------------------*/
/* Get pointer to n-th element.                                              */
/*---------------------------------------------------------------------------*/

void *_unur_slist_replace( struct unur_slist *slist, int n, void *element );
/*---------------------------------------------------------------------------*/
/* Replace (existing) pointer to n-th element by 'element'.                  */
/*---------------------------------------------------------------------------*/

void _unur_slist_free( struct unur_slist *slist );
/*---------------------------------------------------------------------------*/
/* Free all elements and list in simple list.                                */
/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*/
#endif  /* SLIST_H_SEEN */
/*---------------------------------------------------------------------------*/
