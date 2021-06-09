/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   FILE: functparser_struct.h                                              *
 *                                                                           *
 *   PURPOSE:                                                                *
 *         declares structures for function parser                           *
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
/* Structure for function tree                                               */

struct ftreenode { 
  char            *symbol;      /* name of token                             */
  int             token;        /* location of token in list of symbols      */
  int             type;         /* type of token (e.g. S_ADD_OP)             */
  double          val;          /* value of constant or (and)
				   value of node during evalution of tree    */
  struct ftreenode *left;       /* pointer to left branch/leave of node      */
  struct ftreenode *right;      /* pointer to right branch/leave of node     */

#ifdef UNUR_COOKIES
  unsigned cookie;              /* magic cookie                              */
#endif
}; 


/*---------------------------------------------------------------------------*/

