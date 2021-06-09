/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   FILE: distr_struct.h                                                    *
 *                                                                           *
 *   PURPOSE:                                                                *
 *         declares structures for distributions                             *
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
/* define object for univariate continuous distribution                      */

struct unur_distr_cont {

  UNUR_FUNCT_CONT *pdf;         /* pointer to PDF                            */
  UNUR_FUNCT_CONT *dpdf;        /* pointer to derivative of PDF              */
  UNUR_FUNCT_CONT *cdf;         /* pointer to CDF                            */
  UNUR_FUNCT_CONT *invcdf;      /* pointer to inverse of CDF                 */
  UNUR_FUNCT_CONT *logpdf;      /* pointer to logPDF                         */
  UNUR_FUNCT_CONT *dlogpdf;     /* pointer to derivative of logPDF           */
  UNUR_FUNCT_CONT *logcdf;      /* pointer to logCDF                         */
  UNUR_FUNCT_CONT *hr;          /* pointer to hazard rate                    */

  double norm_constant;         /* (log of) normalization constant for PDF   */

  double params[UNUR_DISTR_MAXPARAMS];  /* parameters of the PDF             */
  int    n_params;              /* number of parameters of the PDF           */

  double *param_vecs[UNUR_DISTR_MAXPARAMS];  /* parameter vectors            */
  int    n_param_vec[UNUR_DISTR_MAXPARAMS]; /* lengths of the parameter vecs */

  double mode;                  /* location of mode                          */
  double center;                /* location of center                        */
  double area;                  /* area below PDF                            */
  double domain[2];             /* boundary of domain                        */
  double trunc[2];              /* boundary of truncated domain              */

  struct ftreenode *pdftree;    /* pointer to function tree for PDF          */
  struct ftreenode *dpdftree;   /* pointer to function tree for dPDF         */
  struct ftreenode *logpdftree; /* pointer to function tree for logPDF       */
  struct ftreenode *dlogpdftree;/* pointer to function tree for dlogPDF      */
  struct ftreenode *cdftree;    /* pointer to function tree for CDF          */
  struct ftreenode *logcdftree; /* pointer to function tree for logCDF       */
  struct ftreenode *hrtree;     /* pointer to function tree for hazard rate  */

  int (*set_params)(struct unur_distr *distr, const double *params, int n_params );
                                /* function for setting parameters and domain*/
  int (*upd_mode)(struct unur_distr *distr);
                                /* function for computing mode               */
  int (*upd_area)(struct unur_distr *distr);
                                /* function for computing area               */

  int  (*init)(struct unur_par *par,struct unur_gen *gen);
                                /* pointer to special init routine           */
};

/*---------------------------------------------------------------------------*/
/* define object for multivariate continuous distribution                    */

struct unur_distr_cvec {

  UNUR_FUNCT_CVEC *pdf;         /* pointer to PDF                            */
  UNUR_VFUNCT_CVEC *dpdf;       /* pointer to gradient of PDF                */
  UNUR_FUNCTD_CVEC *pdpdf;      /* pointer to partial derivative of PDF      */
  UNUR_FUNCT_CVEC *logpdf;      /* pointer to logPDF                         */
  UNUR_VFUNCT_CVEC *dlogpdf;    /* pointer to gradient of logPDF             */
  UNUR_FUNCTD_CVEC *pdlogpdf;   /* pointer to partial derivative of logPDF   */

  double *mean;                 /* mean vector of distribution               */
  double *covar;                /* covariance matrix of distribution         */

  double *cholesky;             /* cholesky factor of covariance matrix      */
  double *covar_inv;            /* inverse of covariance matrix              */

  double *rankcorr;             /* rank correlation matrix                   */
  double *rk_cholesky;          /* cholesky factor of rank corr. matrix      */

  struct unur_distr **marginals; /* array of pointers to marginal distributions */

  double params[UNUR_DISTR_MAXPARAMS];  /* parameters of the PDF             */
  int    n_params;              /* number of parameters of the PDF           */
  
  double *param_vecs[UNUR_DISTR_MAXPARAMS];  /* parameter vectors            */
  int    n_param_vec[UNUR_DISTR_MAXPARAMS]; /* lengths of the parameter vecs */
  
  double norm_constant;         /* (log of) normalization constant for PDF   */

  double *mode;                 /* location of mode                          */
  double *center;               /* location of center                        */
  double volume;                /* volume below PDF                          */
  double *domainrect;           /* (rectangular) domain of domain:            
				   the lower left and the upper right vertices
				   are stored componentwise, i.e.
				   l[0], r[0], l[1], r[1], l[2], r[2], ...   */

#ifdef USE_DEPRECATED_CODE
  struct unur_distr **stdmarginals; /* array of pointers to standardized marginal distributions */
#endif

  int (*upd_mode)(struct unur_distr *distr);
                                /* function for computing mode               */
  int (*upd_volume)(struct unur_distr *distr);
                                /* function for computing volume             */

  int (*init)(struct unur_gen *gen);
                                /* pointer to special init routine           */
};

/*---------------------------------------------------------------------------*/
/* define object for univariate discrete distribution                        */

struct unur_distr_discr {
  /* (finite) probability vector */
  double *pv;                   /* pointer to probability vector             */
  int     n_pv;                 /* length of probability vector              */

  /* probability mass function */
  UNUR_FUNCT_DISCR  *pmf;       /* pointer to probability mass function      */
  UNUR_FUNCT_DISCR  *cdf;       /* pointer to CDF                            */
  UNUR_IFUNCT_DISCR *invcdf;    /* pointer to inverse of CDF                 */

  double params[UNUR_DISTR_MAXPARAMS];  /* parameters of the PMF             */
  int    n_params;              /* number of parameters of the PMF           */

  double norm_constant;         /* (log of) normalization constant for PMF   */

  int    mode;                  /* location of mode                          */
  double sum;                   /* sum over PMF                              */

  int (*set_params)(struct unur_distr *distr, const double *params, int n_params );
                                /* function for setting parameters and domain*/
  int (*upd_mode)(struct unur_distr *distr);
                                /* function for computing mode               */
  int (*upd_sum)(struct unur_distr *distr);
                                /* function for computing sum                */

  /* other parameters */
  int domain[2];                /* boundary of domain                        */
  int trunc[2];                 /* boundary of truncated domain              */
  /** trunc[] not supported yet **/

  struct ftreenode *pmftree;    /* pointer to function tree for PMF          */
  struct ftreenode *cdftree;    /* pointer to function tree for CDF          */

  int  (*init)(struct unur_par *par,struct unur_gen *gen);
                                /* pointer to special init routine           */
};

/*---------------------------------------------------------------------------*/
/* define object for empirical univariate constinuous distribution           */
/* (given by empirical sample)                                               */

struct unur_distr_cemp {
  /* raw data */
  int     n_sample;             /* length of sample probability vector       */
  double *sample;               /* pointer to sample                         */
  /* histogram */
  int     n_hist;               /* number of bins in histogram               */
  double *hist_prob;            /* probabilities for bins                    */
  double  hmin, hmax;           /* lower and upper bound for histograms with
				   bins of equal length                      */
  double *hist_bins;            /* boundary between bins (of different length) */
};

/*---------------------------------------------------------------------------*/
/* define object for empirical mulitvariate constinuous distribution         */
/* (given by empirical sample)                                               */

struct unur_distr_cvemp {
  double *sample;              /* pointer to sample                          */
  int    n_sample;             /* length of sample probability vector        */
};

/*---------------------------------------------------------------------------*/
/* define object for matrix distribution                                     */

struct unur_distr_matr {

  int n_rows;                   /* number of rows                            */
  int n_cols;                   /* number of columns                         */

  int  (*init)(struct unur_par *par,struct unur_gen *gen);
                                /* pointer to special init routine           */
};

/*---------------------------------------------------------------------------*/
/* define distribution object                                                */

struct unur_distr {
  union {             
    struct unur_distr_cont  cont;   /* univariate continuous distribution    */
    struct unur_distr_matr  matr;   /* matrix distribution                   */
    struct unur_distr_cvec  cvec;   /* multivariate continuous distribution  */
    struct unur_distr_discr discr;  /* univariate discrete distribution      */
    struct unur_distr_cemp  cemp;   /* empirical univ. cont. distr. (sample) */
    struct unur_distr_cvemp cvemp;  /* empir. multiv. cont. distr. (sample)  */
  } data;                           /* data for distribution                 */

  unsigned type;                    /* type of distribution                  */
  unsigned id;                      /* identifier for distribution           */
  const char *name;                 /* name of distribution                  */
  char *name_str;                   /* string for storing user name of distr */
  int dim;                          /* number of components of random vector */

  unsigned set;                     /* indicate changed parameters           */

  const void *extobj;               /* pointer to an object for additional data */

  struct unur_distr *base;          /* pointer to distribution object for
				       derived distribution 
				       (e.g. order statistics)               */

  void (*destroy)(struct unur_distr *distr); /* pointer to destructor        */
  struct unur_distr* (*clone)(const struct unur_distr *distr ); /* clone     */

#ifdef UNUR_COOKIES
  unsigned cookie;                  /* magic cookie                          */
#endif
};

/*---------------------------------------------------------------------------*/
