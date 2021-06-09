/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   FILE: unur_stddistr.h                                                   *
 *                                                                           *
 *   PURPOSE:                                                                *
 *         defines identifiers for standard distributions                    *
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
#ifndef UNUR_STDDISTR_H_SEEN
#define UNUR_STDDISTR_H_SEEN
/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*/
/* indentifiers for standard distributions                                   */

#define UNUR_DISTR_STD         (0x00000001u)  /* flag for UNU.RAN standard distribution  */

enum {
  UNUR_DISTR_GENERIC          = 0x00000000u,

  UNUR_DISTR_CORDER           = 0x00000010u,  /* order statistics                        */
  UNUR_DISTR_CXTRANS          = 0x00000020u,  /* distribution of transformed RV          */
  UNUR_DISTR_CONDI            = 0x00000030u,  /* full conditional distribution           */

  /**                                              pdf   cdf   mode  area  gen   doc     */
  UNUR_DISTR_BETA             = 0x00000101u,  /*    X     X     X     X     X     .      */
  UNUR_DISTR_CAUCHY           = 0x00000201u,  /*    X     X     X     X     X     .      */
  UNUR_DISTR_CHI              = 0x00000301u,  /*    X     X     X     X     X     .      */
  UNUR_DISTR_CHISQUARE        = 0x00000401u,  /*    X     X     X     X     .     .      */
  UNUR_DISTR_EPANECHNIKOV     = 0x00000501u,  /*    .     .     .     .     .     .      */
  UNUR_DISTR_EXPONENTIAL      = 0x00000601u,  /*    X     X     X     X     X     .      */
  UNUR_DISTR_EXTREME_I        = 0x00000701u,  /*    X     X     X     X     X     .      */
  UNUR_DISTR_EXTREME_II       = 0x00000801u,  /*    X     X     X     X     X     .      */
  UNUR_DISTR_F                = 0x00000901u,  /*    X     X     X     X     X     .      */
  UNUR_DISTR_GAMMA            = 0x00000a01u,  /*    X     X     X     X     X     .      */
  UNUR_DISTR_GHYP             = 0x00002401u,  /*    X     .     .     .     .     .      */
  UNUR_DISTR_GIG              = 0x00000b01u,  /*    X     .     X     .     X     .      */
  UNUR_DISTR_GIG2             = 0x00002201u,  /*    X     .     X     .     X     .      */
  UNUR_DISTR_HYPERBOLIC       = 0x00002301u,  /*    X     .     X     .     .     .      */
  UNUR_DISTR_IG               = 0x00002101u,  /*    X     X     X     .     .     .      */
  UNUR_DISTR_LAPLACE          = 0x00000c01u,  /*    X     X     X     X     X     .      */
  UNUR_DISTR_LOGISTIC         = 0x00000d01u,  /*    X     X     X     X     X     .      */
  UNUR_DISTR_LOGNORMAL        = 0x00000e01u,  /*    X     X     .     X     .     .      */
  UNUR_DISTR_LOMAX            = 0x00000f01u,  /*    X     X     X     X     .     .      */
  UNUR_DISTR_NORMAL           = 0x00001001u,  /*    X     X     X     X     X     .      */
   UNUR_DISTR_GAUSSIAN        = 0x00001001u,  /*    same as NORMAL                       */
  UNUR_DISTR_PARETO           = 0x00001101u,  /*    X     X     X     X     .     .      */
  UNUR_DISTR_POWEREXPONENTIAL = 0x00001201u,  /*    X     X     X     X     X     .      */
  UNUR_DISTR_RAYLEIGH         = 0x00001301u,  /*    X     X     X     X     .     .      */
  UNUR_DISTR_SLASH            = 0x00001401u,  /*    X     .     X     X     X     .      */
  UNUR_DISTR_STUDENT          = 0x00001501u,  /*    X    (X)    X     X     X     .      */
  UNUR_DISTR_TRIANGULAR       = 0x00001601u,  /*    X     X     X     X     X     .      */
  UNUR_DISTR_UNIFORM          = 0x00002001u,  /*    X     X     X     X     .     .      */
   UNUR_DISTR_BOXCAR          = 0x00002001u,  /*    same as UNIFORM                      */
  UNUR_DISTR_WEIBULL          = 0x00001801u,  /*    X     X     X     X     X     .      */

  UNUR_DISTR_BURR_I           = 0x0000b001u,  /*    .     X     .     .     X     .      */
  UNUR_DISTR_BURR_II          = 0x0000b101u,  /*    .     X     .     .     X     .      */
  UNUR_DISTR_BURR_III         = 0x0000b201u,  /*    .     X     .     .     X     .      */
  UNUR_DISTR_BURR_IV          = 0x0000b301u,  /*    .     X     .     .     X     .      */
  UNUR_DISTR_BURR_V           = 0x0000b401u,  /*    .     X     .     .     X     .      */
  UNUR_DISTR_BURR_VI          = 0x0000b501u,  /*    .     X     .     .     X     .      */
  UNUR_DISTR_BURR_VII         = 0x0000b601u,  /*    .     X     .     .     X     .      */
  UNUR_DISTR_BURR_VIII        = 0x0000b701u,  /*    .     X     .     .     X     .      */
  UNUR_DISTR_BURR_IX          = 0x0000b801u,  /*    .     X     .     .     X     .      */
  UNUR_DISTR_BURR_X           = 0x0000b901u,  /*    .     X     .     .     X     .      */
  UNUR_DISTR_BURR_XI          = 0x0000ba01u,  /*    .     X     .     .     .     .      */
  UNUR_DISTR_BURR_XII         = 0x0000bb01u,  /*    .     X     .     .     X     .      */

  /**                                              pmf   cdf   mode  sum   gen   doc     */
  UNUR_DISTR_BINOMIAL         = 0x00010001u,  /*    .     .     .     .     .     .      */
  UNUR_DISTR_GEOMETRIC        = 0x00020001u,  /*    X     X     X     X     X     .      */
  UNUR_DISTR_HYPERGEOMETRIC   = 0x00030001u,  /*    .     .     .     .     .     .      */
  UNUR_DISTR_LOGARITHMIC      = 0x00040001u,  /*    X     .     X     X     X     .      */
  UNUR_DISTR_NEGATIVEBINOMIAL = 0x00050001u,  /*    X     .     X     X     X     .      */
  UNUR_DISTR_POISSON          = 0x00060001u,  /*    X     X     X     X     X     .      */
  UNUR_DISTR_ZIPF             = 0x00070001u,  /*    X     .     X     .     X     .      */

  /**                                              pdf                           doc     */
  UNUR_DISTR_MCAUCHY          = 0x01000001u,  /*    X                             .      */
  UNUR_DISTR_MNORMAL          = 0x02000001u,  /*    X                             .      */
  UNUR_DISTR_MSTUDENT         = 0x03000001u,  /*    X                             .      */
  UNUR_DISTR_MEXPONENTIAL     = 0x04000001u,  /*    X                             .      */
  UNUR_DISTR_COPULA           = 0x05000001u,  /*    .                             .      */

  /**                                                                        */
  UNUR_DISTR_MCORRELATION     = 0x10000001u   /*                                         */

};

/*---------------------------------------------------------------------------*/
#endif  /* UNUR_STDDISTR_H_SEEN */
/*---------------------------------------------------------------------------*/
