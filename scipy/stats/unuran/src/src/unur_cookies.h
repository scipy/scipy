/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   FILE: source_cookies.h                                                  *
 *                                                                           *
 *   PURPOSE:                                                                *
 *         defines magic cookies.                                            *
 *                                                                           *
 *   USAGE:                                                                  *
 *         internal header file                                              *
 *          included only in source_unuran.h                                 *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   Copyright (c) 2000-2010 Wolfgang Hoermann and Josef Leydold             *
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
#ifndef SOURCE_COOKIES_H_SEEN
#define SOURCE_COOKIES_H_SEEN
/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*/
#ifdef UNUR_COOKIES    /* use magic cookies                                  */
/*---------------------------------------------------------------------------*/
/* name of cookies                                                           */

/* generators discrete distributions */
#define CK_DARI_PAR      0x00000010u
#define CK_DARI_GEN      0x00000011u
#define CK_DAU_PAR       0x00000020u
#define CK_DAU_GEN       0x00000021u
#define CK_DGT_PAR       0x00000030u
#define CK_DGT_GEN       0x00000031u
#define CK_DSROU_PAR     0x00000040u
#define CK_DSROU_GEN     0x00000041u
#define CK_DSS_PAR       0x00000050u
#define CK_DSS_GEN       0x00000051u

/* generators continuous distributions */
#define CK_AROU_PAR      0x00100010u
#define CK_AROU_GEN      0x00100011u
#define CK_AROU_SEG      0x00100012u
#define CK_ARS_PAR       0x001000d0u
#define CK_ARS_GEN       0x001000d1u
#define CK_ARS_IV        0x001000d2u
#define CK_HINV_PAR      0x00100020u
#define CK_HINV_GEN      0x00100021u
#define CK_HINV_IV       0x00100022u
#define CK_HRB_PAR       0x00100030u
#define CK_HRB_GEN       0x00100031u
#define CK_HRD_PAR       0x00100040u
#define CK_HRD_GEN       0x00100041u
#define CK_HRI_PAR       0x00100050u
#define CK_HRI_GEN       0x00100051u
#define CK_ITDR_PAR      0x00100080u
#define CK_ITDR_GEN      0x00100081u
#define CK_NINV_PAR      0x00100060u
#define CK_NINV_GEN      0x00100061u
#define CK_NROU_PAR      0x00100070u
#define CK_NROU_GEN      0x00100071u
#define CK_PINV_PAR      0x00100130u
#define CK_PINV_GEN      0x00100131u
#define CK_PINV_IV       0x00100132u
#define CK_SROU_PAR      0x00100090u
#define CK_SROU_GEN      0x00100091u
#define CK_SSR_PAR       0x001000a0u
#define CK_SSR_GEN       0x001000a1u
#define CK_TABL_PAR      0x001000b0u
#define CK_TABL_GEN      0x001000b1u
#define CK_TABL_IV       0x001000b2u
#define CK_TDR_PAR       0x001000c0u
#define CK_TDR_GEN       0x001000c1u
#define CK_TDR_IV        0x001000c2u
#define CK_UNIF_PAR      0x001000e0u
#define CK_UNIF_GEN      0x001000e1u
#define CK_UTDR_PAR      0x001000f0u
#define CK_UTDR_GEN      0x001000f1u
#define CK_EMPK_PAR      0x00100100u
#define CK_EMPK_GEN      0x00100101u
#define CK_EMPL_PAR      0x00100110u
#define CK_EMPL_GEN      0x00100111u
#define CK_HIST_PAR      0x00100120u
#define CK_HIST_GEN      0x00100121u
/* meta generators */
#define CK_MIXT_PAR      0x00110130u
#define CK_MIXT_GEN      0x00110131u

/* generators multivariate continuous distributions */
#define CK_MVTDR_PAR     0x00200010u
#define CK_MVTDR_GEN     0x00200010u
#define CK_VMT_PAR       0x00200020u
#define CK_VMT_GEN       0x00200021u
#define CK_VEMPK_PAR     0x00200030u
#define CK_VEMPK_GEN     0x00200031u
#define CK_VNROU_PAR     0x00200040u
#define CK_VNROU_GEN     0x00200041u
#define CK_VAROU_PAR     0x00200050u
#define CK_VAROU_GEN     0x00200051u
#define CK_NORTA_PAR     0x00200060u
#define CK_NORTA_GEN     0x00200061u
/* Markov chain samplers */
#define CK_BALL_PAR      0x00201010u
#define CK_BALL_GEN      0x00201011u
#define CK_GIBBS_PAR     0x00201020u
#define CK_GIBBS_GEN     0x00201021u
#define CK_HITRO_PAR     0x00201030u
#define CK_HITRO_GEN     0x00201031u
#define CK_WALK_PAR      0x00201040u
#define CK_WALK_GEN      0x00201041u

/* generators for random matrices */
#define CK_MCORR_PAR     0x00400010u
#define CK_MCORR_GEN     0x00400011u

/* special generators */
#define CK_CSTD_PAR      0x10000010u
#define CK_CSTD_GEN      0x10000011u
#define CK_DSTD_PAR      0x10000020u
#define CK_DSTD_GEN      0x10000021u
#define CK_MVSTD_PAR     0x10000040u
#define CK_MVSTD_GEN     0x10000041u

#define CK_CEXT_PAR      0x10000110u
#define CK_CEXT_GEN      0x10000111u
#define CK_DEXT_PAR      0x10000120u
#define CK_DEXT_GEN      0x10000121u

#define CK_MBLOCK        0xf0000001u

/* misc methods */
#define CK_AUTO_PAR      0xa0000010u

/* distribution objects */
#define CK_DISTR         0xe0000000u
#define CK_DISTR_CONT    0xe0000001u
#define CK_DISTR_CEMP    0xe0000002u
#define CK_DISTR_MATR    0xe0000003u
#define CK_DISTR_CVEC    0xe0000004u
#define CK_DISTR_CVEMP   0xe0000005u
#define CK_DISTR_DISCR   0xe0000006u

#define CK_SPECIALGEN_CONT 0xd00001u

/* URNG (uniform random number generator objects) */
#define CK_URNG          0xb0000001u

/* function parser */
#define CK_FSTR_PDATA    0xa0000001u
#define CK_FSTR_TNODE    0xa0000002u

/* auxiliary tools */
#define CK_SLIST         0xa0000004u

/*---------------------------------------------------------------------------*/
/* macros for dealing with magic cookies                                     */

/* set magic cookies */
#define COOKIE_SET(ptr,ck)           (ptr)->cookie=(ck)

#define COOKIE_SET_ARRAY(ptr,ck,n) { \
  register int i; \
  for (i=0;i<(n);i++) \
    ((ptr)+(i))->cookie=(ck); \
}

/* check magic cookies */
#define COOKIE_CHECK(ptr,ck,rval) \
  if((ptr)->cookie!=(ck)) { \
    _unur_error_cookies(__FILE__,__LINE__,(ptr)->cookie, (ck)); \
    return rval; \
  }

/* clear magic cookies (set to 0u). To be used in connection with free() to  */
/* to detect access to freed object.                                         */                     
#define COOKIE_CLEAR(ptr)           (ptr)->cookie=0u

/*---------------------------------------------------------------------------*/
#else                               /* do not use magic cookies              */
/*---------------------------------------------------------------------------*/

#define COOKIE_SET(ptr,ck) 
#define COOKIE_SET_ARRAY(ptr,ck,n)
#define COOKIE_CHECK(ptr,ck,rval) 
#define COOKIE_CLEAR(ptr)

/*---------------------------------------------------------------------------*/
#endif 
/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*/
#endif  /* SOURCE_COOKIES_H_SEEN */
/*---------------------------------------------------------------------------*/
