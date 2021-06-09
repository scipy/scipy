/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   file:      quantiles.c                                                  *
 *                                                                           *
 *   compute estimate for quartiles and mean of samples                      *
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
 *****************************************************************************
 *                                                                           *
 *   REFERENCES:                                                             *
 *   [1] R. Jain, I. Chlamtac (1985): The P^2 Algorithm for Dynamic          *
 *       Calculation of Quantiles and Histograms Without Storing             *
 *       Observations, Comm. ACM 28(19), pp. 1076-1085.                      *
 *                                                                           *
 *****************************************************************************/

/*---------------------------------------------------------------------------*/

#include <unur_source.h>
#include <methods/unur_methods_source.h>
#include <methods/x_gen_source.h>
#include "unuran_tests.h"

/*---------------------------------------------------------------------------*/
static char test_name[] = "Quantiles";
/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*/

int
unur_test_quartiles( UNUR_GEN *gen, double *q0 ,double *q1, double *q2, double *q3, double *q4,
		     int samplesize, int verbosity, FILE *out )
     /*----------------------------------------------------------------------*/
     /*  compute estimate for quartiles and mean of samples.                 */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen        ... pointer to generator object                         */
     /*   q0         ... minimum                                             */
     /*   q1         ... 1st quartile (25%)                                  */
     /*   q2         ... median (2nd quartile, 50%)                          */
     /*   q3         ... 3rd quartile (75%)                                  */
     /*   q4         ... maximum                                             */
     /*   samplesize ... sample size                                         */
     /*   verbosity  ... verbosity level, 0 = no output, 1 = output          */
     /*   out        ... output stream                                       */
     /*                                                                      */
     /* return:                                                              */
     /*   UNUR_SUCCESS ... on success                                        */
     /*   error code   ... on error                                          */
     /*----------------------------------------------------------------------*/
{
  double x = 0.;
  int n;

  double  p = 0.5; /* p-quantile -- median        */
  double  h[5];    /* marker height               */
  int     pos[5];  /* marker position             */
  double  dpos[5]; /* desired position            */
  int     i, j;    /* loop variables              */
  int     sgnd;    /* sign of d                   */
  double  tmp;
  double  d;       /* (dp-p) at first, then sign  */


  /* check parameter */
  _unur_check_NULL(test_name, gen, UNUR_ERR_NULL);

  /* type of distribution */
  if (! ( ((gen->method & UNUR_MASK_TYPE) == UNUR_METH_DISCR) ||
	  ((gen->method & UNUR_MASK_TYPE) == UNUR_METH_CONT) )) {
    _unur_error(test_name,UNUR_ERR_GENERIC,"dont know how to compute quartiles for distribution");
    return UNUR_ERR_GENERIC;
  }

  /* sample size >= 10 */
  if (samplesize < 10) 
    samplesize = 10;

  /* sampling */
  /* estimate quartiles using algorithm [1]. */
  
  for (n=0; n<samplesize; n++) {

    /* which type of distribution */
    switch (gen->method & UNUR_MASK_TYPE) {
    case UNUR_METH_DISCR:
      x = _unur_sample_discr(gen); break;
    case UNUR_METH_CONT:
      x = _unur_sample_cont(gen); break;
    }

    /* ******************************************* */
    /* PP-algorithm                                */


    if ( n==0 ){  /* BEGIN OF INITIALIZATION */

      h[n]    = x;
      pos[n]  = n;

      /* initialize desired positions   */
      dpos[0] = 0.;
      dpos[1] = p/2.;
      dpos[2] = 4.*p;
      dpos[3] = 2. + 2.*p;
      dpos[4] = 4.;

    }                  /* end of n==0   */
    else if ( n<4 ){   /* 0<n<4         */

      h[n]   = x;
      pos[n] = n;

    }                  /* end of 0<n<4  */
    else if ( n==4 ){

      h[n] = x;
      pos[n] = n;

      /* sort the first five elements */
      for (i=0; i<4; i++){
	for (j=0; j<4-i; j++){
	  if ( h[j] > h[j+1] ){  /* swap */
	    tmp    = h[j];
	    h[j]   = h[j+1];
	    h[j+1] = tmp;
	  }
	}
      }  /* end of sorting */

    }      /* end of n==4  -- END OF INITIALIZATION */
    else{  /* n > 4  */


      /* adjust height of extreme markers */
      if ( x < h[0] )
	h[0] = x;
      if ( x > h[4] )
	h[4] = x;


      /* adjust position of all markers */
      for (i=1; i<4; i++){
	if ( x < h[i] ){
	  pos[i]++;
	}
      } /* end of for */
      pos[4]++; /* will always be at the last positon */


      /* adjust desired position of all markers */
      dpos[1] = (double) n * p/2.;
      dpos[2] = (double) n * p;
      dpos[3] = (double) n * (1.+p)/2.;
      dpos[4] = (double) n;


      /* adjust height of markers if necessary */
      for (i=1; i<4; i++){

	d = dpos[i] - pos[i];  /* discrepancy of position*/
    
	/* decide on update of markers */
        if ( (d >=  1.  &&  pos[i+1]-pos[i] >  1) ||
             (d <= -1.  &&  pos[i-1]-pos[i] < -1)   ){

	  /* markers are shifted by at most 1 */
	  sgnd = (d < 0.) ? -1 : 1;
	  d = (double) sgnd;

	  /* try parabolic formula */
          tmp = h[i] + d/(pos[i+1]-pos[i-1]) * 
	       ( (pos[i]-pos[i-1]+d) * (h[i+1]-h[i])/(pos[i+1]-pos[i]) +
                 (pos[i+1]-pos[i]-d) * (h[i]-h[i-1])/(pos[i]-pos[i-1]) );


	  /* use linear formula -- parabolic model failed */
	  if ( tmp <= h[i-1] || tmp >= h[i+1]){
	    tmp = h[i] + d*(h[i]-h[i+sgnd])/(pos[i]-pos[i+sgnd]);
	  }

	  h[i]    = tmp;      /* update height   */
	  pos[i] += sgnd;     /* update position */
       
	} /* end if -- marker i updated */

      } /* end for -- all markers updated  */


/* 
      for (i=0;i<5;i++)
	fprintf(out,"%6.2f ", h[i]);
      fprintf(out,"\n");
*/

    } /* end n>4 */

    /* end of PP-algorithm                         */
    /* ******************************************* */
 
  }
  
  /* copy results */
  *q0 = h[0];
  *q1 = h[1];
  *q2 = h[2];
  *q3 = h[3];
  *q4 = h[4];

  /* now print results */
  if (verbosity) {
    fprintf(out,"\nQuartiles:\n");
    fprintf(out,"\tmin = \t%6.5g\n",*q0);
    fprintf(out,"\t25%% =\t%6.5g\n",*q1);
    fprintf(out,"\t50%% =\t%6.5g\n",*q2);
    fprintf(out,"\t75%% =\t%6.5g\n",*q3);
    fprintf(out,"\tmax = \t%6.5g\n",*q4);
  }

  return UNUR_SUCCESS;

} /* end of unur_test_quartiles() */

/*---------------------------------------------------------------------------*/





