/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   FILE: cemp.h                                                            *
 *                                                                           *
 *   PURPOSE:                                                                *
 *         function prototypes for manipulating distribution objects of      *
 *         type  CEMP  (continuous empirical univariate distribution)        *
 *                                                                           *
 *   USAGE:                                                                  *
 *         only included in unuran.h                                         *
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

/* 
   =NODEX   CEMP   Continuous empirical univariate distributions

   =UP Distribution_objects [20]

   =DESCRIPTION
      Empirical univariate distributions are derived from observed data.
      There are two ways to create such a generator object:
      @enumerate
      @item
         By a list of @emph{raw data} by means of a
         unur_distr_cemp_set_data() call.
      @item
         By a @emph{histogram} (i.e. preprocessed data) by means of a
	 unur_distr_cemp_set_hist() call.
      @end enumerate
      How these data are used to sample from the empirical distribution
      depends from the chosen generation method.  

   =END
*/

/*---------------------------------------------------------------------------*/

/* 
   Routines for handling empirical univariate continuous distributions (CEMP).
*/

/* =ROUTINES */

UNUR_DISTR *unur_distr_cemp_new( void );
/* 
   Create a new (empty) object for empirical univariate continuous distribution.
*/

/* ==DOC
   @subsubheading Essential parameters
*/

int unur_distr_cemp_set_data( UNUR_DISTR *distribution, const double *sample, int n_sample );
/* 
   Set observed sample for empirical distribution.
*/

int unur_distr_cemp_read_data( UNUR_DISTR *distribution, const char *filename );
/* 
   Read data from file @file{filename}.
   It reads the first number from each line. 
   Numbers are parsed by means of the C standard routine @command{strtod}.
   Lines that do not start with @code{+}, @code{-}, @code{.}, or a
   digit are ignored. (Beware of lines starting with a blank!)

   In case of an error (file cannot be opened, invalid string for
   double in line) no data are copied into the distribution object
   and an error code is returned.
*/

int unur_distr_cemp_get_data( const UNUR_DISTR *distribution, const double **sample );
/* 
   Get number of samples and set pointer @var{sample} to array of
   observations. If no sample has been given, an error code 
   is returned and @code{sample} is set to NULL.

   @emph{Important:} Do @strong{not} change the entries in @var{sample}!
*/

int unur_distr_cemp_set_hist( UNUR_DISTR *distribution, const double *prob, int n_prob, double xmin, double xmax );
/* 
   Set a histogram with bins of equal width. @var{prob} is an array
   of length @var{n_prob} that contains the probabilities for the bins
   (in ascending order). @var{xmin} and @var{xmax} give the lower and
   upper bound of the histogram, respectively. The bins are assumed to
   have equal width.

   @emph{Remark:} This is shortcut for calling
   unur_distr_cemp_set_hist_prob() and
   unur_distr_cemp_set_hist_domain().

   @emph{Notice:} All sampling methods either use raw data or histogram.
   It is possible to set both types of data; however, it is not
   checked whether the given histogran corresponds to possibly given
   raw data.
*/

int unur_distr_cemp_set_hist_prob( UNUR_DISTR *distribution, const double *prob, int n_prob );
/* 
   Set probabilities of a histogram with @var{n_prob} bins.
   Hence @var{prob} must be an array of length @var{n_prob} that
   contains the probabilities for the bins in ascending order.
   It is important also to set the location of the bins either
   with a unur_distr_cemp_set_hist_domain() for bins of equal
   width or unur_distr_cemp_set_hist_bins() when the bins have
   different width.

   @emph{Notice:} All sampling methods either use raw data or histogram.
   It is possible to set both types of data; however, it is not
   checked whether the given histogram corresponds to possibly given
   raw data.
*/

int unur_distr_cemp_set_hist_domain( UNUR_DISTR *distribution, double xmin, double xmax );
/* 
   Set a domain of a histogram with bins of equal width.
   @var{xmin} and @var{xmax} give the lower and upper bound of the
   histogram, respectively.
*/

int unur_distr_cemp_set_hist_bins( UNUR_DISTR *distribution, const double *bins, int n_bins );
/* 
   Set location of bins of a histogram with @var{n_bins} bins.
   Hence @var{bins} must be an array of length @var{n_bins}.
   The domain of the @var{distribution} is automatically set by 
   this call and overrides any calls to
   unur_distr_cemp_set_hist_domain().

   @emph{Important:}
   The probabilities of the bins of the @var{distribution} must be
   already be set by a unur_distr_cemp_set_hist_prob() 
   (or a unur_distr_cemp_set_hist() call) and the value of 
   @var{n_bins} must equal @var{n_prob}@code{+1} from the
   corresponding value of the respective call.
*/


/* =END */

/*---------------------------------------------------------------------------*/
