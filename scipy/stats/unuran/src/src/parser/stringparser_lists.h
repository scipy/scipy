/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   FILE:      stringparser_lists.ch.in / stringparser_lists.ch             *
 *                                                                           *
 *                                                                           *
 *   DESCRIPTION:                                                            *
 *                                                                           *
 *   Switch lists for string parser.                                         *
 *   (See file stringparser.c)                                               *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   This file is parsed by the perl script make_stringparser.pl which       *
 *   replaces the                                                            *
 *        =INPUT keyword                                                     *
 *   tags with search lists created from the information found within the    *
 *   header files of the source code.                                        *
 *   These lists (implemented via 'if' rules together with switch lists to   *
 *   with the first letter of the set calls as simple hash function) are     *
 *   used to call the corresponding ..._set and ..._new calls for the        *
 *   the keywords found in the string.                                       *
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


/*****************************************************************************/
/**  Distributions                                                          **/
/*****************************************************************************/

struct unur_distr *
_unur_str_distr_new( char *distribution )
     /*----------------------------------------------------------------------*/
     /* get new distribution object                                          */
     /*                                                                      */
     /* parameters:                                                          */
     /*   distribution ... string that contains distribution name            */
     /*   params       ... string that contains list of parameters           */
     /*                                                                      */
     /* return:                                                              */
     /*   distribution object (pointer to structure)                         */
     /*                                                                      */
     /* error:                                                               */
     /*   return NULL                                                        */
     /*----------------------------------------------------------------------*/
{
  struct unur_distr *distr = NULL;    /* pointer to distribution object      */
  char distr_unknown;

  char *name;                         /* pointer to name of distribution     */
  char *params;                       /* pointer to parameter list of distr  */

  double *darray = NULL;              /* array of arguments for distribution */
  int n_darray = 0;                   /* size of array of parameters         */

  /* name of distribution */
  name = distribution;

  /* get parameter list */
  params = strchr(distribution,'(');
  if (params != NULL) {
    *params = '\0';                   /* terminate key string        */
    ++params;                         /* set pointer to value string */
  }

  /* get parameter list */
  n_darray = _unur_parse_dlist(params, &darray );

#ifdef UNUR_ENABLE_LOGGING
  /* write info into LOG file */
  if (_unur_default_debugflag & UNUR_DEBUG_SETUP)
    _unur_str_debug_distr(1,name,darray,n_darray);
#endif

  /* mark distribution as unknown (this is a very ugly hack) */
  distr = (struct unur_distr *) &distr_unknown;

	 switch (*distribution) {
	 case 'b':
		 if ( !strcmp( distribution, "beta") ) {
			 distr = unur_distr_beta (darray,n_darray);
			 break;
		 }
		 if ( !strcmp( distribution, "binomial") ) {
			 distr = unur_distr_binomial (darray,n_darray);
			 break;
		 }
		 break;
	 case 'c':
		 if ( !strcmp( distribution, "cauchy") ) {
			 distr = unur_distr_cauchy (darray,n_darray);
			 break;
		 }
		 if ( !strcmp( distribution, "chi") ) {
			 distr = unur_distr_chi (darray,n_darray);
			 break;
		 }
		 if ( !strcmp( distribution, "chisquare") ) {
			 distr = unur_distr_chisquare (darray,n_darray);
			 break;
		 }
		 break;
	 case 'e':
		 if ( !strcmp( distribution, "exponential") ) {
			 distr = unur_distr_exponential (darray,n_darray);
			 break;
		 }
		 if ( !strcmp( distribution, "extremei") ) {
			 distr = unur_distr_extremeI (darray,n_darray);
			 break;
		 }
		 if ( !strcmp( distribution, "extremeii") ) {
			 distr = unur_distr_extremeII (darray,n_darray);
			 break;
		 }
		 break;
	 case 'f':
		 if ( !strcmp( distribution, "f") ) {
			 distr = unur_distr_F (darray,n_darray);
			 break;
		 }
		 break;
	 case 'g':
		 if ( !strcmp( distribution, "gamma") ) {
			 distr = unur_distr_gamma (darray,n_darray);
			 break;
		 }
		 if ( !strcmp( distribution, "geometric") ) {
			 distr = unur_distr_geometric (darray,n_darray);
			 break;
		 }
		 if ( !strcmp( distribution, "gig") ) {
			 distr = unur_distr_gig (darray,n_darray);
			 break;
		 }
		 if ( !strcmp( distribution, "gig2") ) {
			 distr = unur_distr_gig2 (darray,n_darray);
			 break;
		 }
		 break;
	 case 'h':
		 if ( !strcmp( distribution, "hyperbolic") ) {
			 distr = unur_distr_hyperbolic (darray,n_darray);
			 break;
		 }
		 if ( !strcmp( distribution, "hypergeometric") ) {
			 distr = unur_distr_hypergeometric (darray,n_darray);
			 break;
		 }
		 break;
	 case 'i':
		 if ( !strcmp( distribution, "ig") ) {
			 distr = unur_distr_ig (darray,n_darray);
			 break;
		 }
		 break;
	 case 'l':
		 if ( !strcmp( distribution, "laplace") ) {
			 distr = unur_distr_laplace (darray,n_darray);
			 break;
		 }
		 if ( !strcmp( distribution, "logarithmic") ) {
			 distr = unur_distr_logarithmic (darray,n_darray);
			 break;
		 }
		 if ( !strcmp( distribution, "logistic") ) {
			 distr = unur_distr_logistic (darray,n_darray);
			 break;
		 }
		 if ( !strcmp( distribution, "lognormal") ) {
			 distr = unur_distr_lognormal (darray,n_darray);
			 break;
		 }
		 if ( !strcmp( distribution, "lomax") ) {
			 distr = unur_distr_lomax (darray,n_darray);
			 break;
		 }
		 break;
	 case 'n':
		 if ( !strcmp( distribution, "negativebinomial") ) {
			 distr = unur_distr_negativebinomial (darray,n_darray);
			 break;
		 }
		 if ( !strcmp( distribution, "normal") ) {
			 distr = unur_distr_normal (darray,n_darray);
			 break;
		 }
		 break;
	 case 'p':
		 if ( !strcmp( distribution, "pareto") ) {
			 distr = unur_distr_pareto (darray,n_darray);
			 break;
		 }
		 if ( !strcmp( distribution, "poisson") ) {
			 distr = unur_distr_poisson (darray,n_darray);
			 break;
		 }
		 if ( !strcmp( distribution, "powerexponential") ) {
			 distr = unur_distr_powerexponential (darray,n_darray);
			 break;
		 }
		 break;
	 case 'r':
		 if ( !strcmp( distribution, "rayleigh") ) {
			 distr = unur_distr_rayleigh (darray,n_darray);
			 break;
		 }
		 break;
	 case 's':
		 if ( !strcmp( distribution, "slash") ) {
			 distr = unur_distr_slash (darray,n_darray);
			 break;
		 }
		 if ( !strcmp( distribution, "student") ) {
			 distr = unur_distr_student (darray,n_darray);
			 break;
		 }
		 break;
	 case 't':
		 if ( !strcmp( distribution, "triangular") ) {
			 distr = unur_distr_triangular (darray,n_darray);
			 break;
		 }
		 break;
	 case 'u':
		 if ( !strcmp( distribution, "uniform") ) {
			 distr = unur_distr_uniform (darray,n_darray);
			 break;
		 }
		 break;
	 case 'w':
		 if ( !strcmp( distribution, "weibull") ) {
			 distr = unur_distr_weibull (darray,n_darray);
			 break;
		 }
	 }

	 /* get pointer to generic distribution object */
	 if (distr == (struct unur_distr *) &distr_unknown) { 
		 do {
			 if ( !strcmp( distribution, "cemp") ) {
				 distr = unur_distr_cemp_new();
				 break;
			 }
			 if ( !strcmp( distribution, "cont") ) {
				 distr = unur_distr_cont_new();
				 break;
			 }
			 if ( !strcmp( distribution, "discr") ) {
				 distr = unur_distr_discr_new();
				 break;
			 }
		 } while (0);
	 }


   if (distr == (struct unur_distr *) &distr_unknown) {
     /* unknown method */
     _unur_error_unknown(distribution,"distribution");
     distr = NULL;
   }
   else if (distr == NULL) {
     /* invalid data for chosen method */
     _unur_error_invalid(distribution,"distribution");
   }

   /* clear memory */
   if (darray) free(darray);

   /* return result */
   return distr;

} /* end of _unur_str_distr_new() */

/*---------------------------------------------------------------------------*/

int
_unur_str_distr_set( UNUR_DISTR **ptr_distr, const char *key, char *value )
     /*----------------------------------------------------------------------*/
     /* set parameters for distribution                                      */
     /*                                                                      */
     /* it also makes a distribution object for an order statistics for      */
     /* the given distribution when the key word "orderstatistics" occurs.   */
     /* Thus the distribution object itself might be changed.                */
     /*                                                                      */
     /* parameters:                                                          */
     /*   ptr_distr ... holds pointer to distribution object                 */
     /*   key       ... string that contains key for parameter               */
     /*   value     ... string that contains list of arguments for parameter */
     /*                                                                      */
     /* return:                                                              */
     /*   UNUR_SUCCESS ... on success                                        */
     /*   error code   ... on error                                          */
     /*                                                                      */
     /* error:                                                               */
     /*   return error code                                                  */
     /*----------------------------------------------------------------------*/
{
  int result;      /* result of UNU.RAN set call */

  /* derefence pointer to distribution object */
  struct unur_distr *distr = *ptr_distr;

  /* storing arguments of set calls: */
  char type_args[MAX_SET_ARGS+1];  /* array containing info about arguments */
  char *args[MAX_SET_ARGS+1];      /* pointers to argument strings */

  /* tokenize argument string */
  if (_unur_str_set_args( value, type_args, args, MAX_SET_ARGS ) < 0) {
    /* error */
    return UNUR_ERR_STR_SYNTAX;
  }

  /* set result indicator to unknown */
  result = UNUR_ERR_STR_UNKNOWN;

  /* find and execute set call */

	 switch (distr->type) {
	 case UNUR_DISTR_CEMP:
		 switch (*key) {
		 case 'd':
			 if ( !strcmp(key, "data") ) {
				 /* n = 2; type = Di:  UNUR_DISTR *distribution, const double *sample, int n_sample */
				 result = _unur_str_distr_set_Di(distr,key,type_args,args,unur_distr_cemp_set_data);
				 break;
			 }
			 break;
		 case 'h':
			 if ( !strcmp(key, "hist_bins") ) {
				 /* n = 2; type = Di:  UNUR_DISTR *distribution, const double *bins, int n_bins */
				 result = _unur_str_distr_set_Di(distr,key,type_args,args,unur_distr_cemp_set_hist_bins);
				 break;
			 }
			 if ( !strcmp(key, "hist_domain") ) {
				 /* n = 2; type = dd:  UNUR_DISTR *distribution, double xmin, double xmax */
				 result = _unur_str_distr_set_dd(distr,key,type_args,args,unur_distr_cemp_set_hist_domain);
				 break;
			 }
			 if ( !strcmp(key, "hist_prob") ) {
				 /* n = 2; type = Di:  UNUR_DISTR *distribution, const double *prob, int n_prob */
				 result = _unur_str_distr_set_Di(distr,key,type_args,args,unur_distr_cemp_set_hist_prob);
				 break;
			 }
		 }
		 break;
	 case UNUR_DISTR_CONT:
		 switch (*key) {
		 case 'c':
			 if ( !strcmp(key, "cdf") ) {
				 /* n = 1; type = C:  UNUR_DISTR *distribution, const char *cdfstr */
				 result = _unur_str_distr_set_C(distr,key,type_args,args,unur_distr_cont_set_cdfstr);
				 break;
			 }
			 if ( !strcmp(key, "cdfstr") ) {
				 /* n = 1; type = C:  UNUR_DISTR *distribution, const char *cdfstr */
				 result = _unur_str_distr_set_C(distr,key,type_args,args,unur_distr_cont_set_cdfstr);
				 break;
			 }
			 if ( !strcmp(key, "center") ) {
				 /* n = 1; type = d:  UNUR_DISTR *distribution, double center */
				 result = _unur_str_distr_set_d(distr,key,type_args,args,unur_distr_cont_set_center);
				 break;
			 }
			 break;
		 case 'd':
			 if ( !strcmp(key, "domain") ) {
				 /* n = 2; type = dd:  UNUR_DISTR *distribution, double left, double right */
				 result = _unur_str_distr_set_dd(distr,key,type_args,args,unur_distr_cont_set_domain);
				 break;
			 }
			 break;
		 case 'h':
			 if ( !strcmp(key, "hr") ) {
				 /* n = 1; type = C:  UNUR_DISTR *distribution, const char *hrstr */
				 result = _unur_str_distr_set_C(distr,key,type_args,args,unur_distr_cont_set_hrstr);
				 break;
			 }
			 if ( !strcmp(key, "hrstr") ) {
				 /* n = 1; type = C:  UNUR_DISTR *distribution, const char *hrstr */
				 result = _unur_str_distr_set_C(distr,key,type_args,args,unur_distr_cont_set_hrstr);
				 break;
			 }
			 break;
		 case 'l':
			 if ( !strcmp(key, "logcdf") ) {
				 /* n = 1; type = C:  UNUR_DISTR *distribution, const char *logcdfstr */
				 result = _unur_str_distr_set_C(distr,key,type_args,args,unur_distr_cont_set_logcdfstr);
				 break;
			 }
			 if ( !strcmp(key, "logcdfstr") ) {
				 /* n = 1; type = C:  UNUR_DISTR *distribution, const char *logcdfstr */
				 result = _unur_str_distr_set_C(distr,key,type_args,args,unur_distr_cont_set_logcdfstr);
				 break;
			 }
			 if ( !strcmp(key, "logpdf") ) {
				 /* n = 1; type = C:  UNUR_DISTR *distribution, const char *logpdfstr */
				 result = _unur_str_distr_set_C(distr,key,type_args,args,unur_distr_cont_set_logpdfstr);
				 break;
			 }
			 if ( !strcmp(key, "logpdfstr") ) {
				 /* n = 1; type = C:  UNUR_DISTR *distribution, const char *logpdfstr */
				 result = _unur_str_distr_set_C(distr,key,type_args,args,unur_distr_cont_set_logpdfstr);
				 break;
			 }
			 break;
		 case 'm':
			 if ( !strcmp(key, "mode") ) {
				 /* n = 1; type = d:  UNUR_DISTR *distribution, double mode */
				 result = _unur_str_distr_set_d(distr,key,type_args,args,unur_distr_cont_set_mode);
				 break;
			 }
			 break;
		 case 'p':
			 if ( !strcmp(key, "pdf") ) {
				 /* n = 1; type = C:  UNUR_DISTR *distribution, const char *pdfstr */
				 result = _unur_str_distr_set_C(distr,key,type_args,args,unur_distr_cont_set_pdfstr);
				 break;
			 }
			 if ( !strcmp(key, "pdfarea") ) {
				 /* n = 1; type = d:  UNUR_DISTR *distribution, double area */
				 result = _unur_str_distr_set_d(distr,key,type_args,args,unur_distr_cont_set_pdfarea);
				 break;
			 }
			 if ( !strcmp(key, "pdfparams") ) {
				 /* n = 2; type = Di:  UNUR_DISTR *distribution, const double *params, int n_params */
				 result = _unur_str_distr_set_Di(distr,key,type_args,args,unur_distr_cont_set_pdfparams);
				 break;
			 }
			 if ( !strcmp(key, "pdfstr") ) {
				 /* n = 1; type = C:  UNUR_DISTR *distribution, const char *pdfstr */
				 result = _unur_str_distr_set_C(distr,key,type_args,args,unur_distr_cont_set_pdfstr);
				 break;
			 }
		 }
		 break;
	 case UNUR_DISTR_DISCR:
		 switch (*key) {
		 case 'c':
			 if ( !strcmp(key, "cdf") ) {
				 /* n = 1; type = C:  UNUR_DISTR *distribution, const char *cdfstr */
				 result = _unur_str_distr_set_C(distr,key,type_args,args,unur_distr_discr_set_cdfstr);
				 break;
			 }
			 if ( !strcmp(key, "cdfstr") ) {
				 /* n = 1; type = C:  UNUR_DISTR *distribution, const char *cdfstr */
				 result = _unur_str_distr_set_C(distr,key,type_args,args,unur_distr_discr_set_cdfstr);
				 break;
			 }
			 break;
		 case 'd':
			 if ( !strcmp(key, "domain") ) {
				 /* n = 2; type = ii:  UNUR_DISTR *distribution, int left, int right */
				 result = _unur_str_distr_set_ii(distr,key,type_args,args,unur_distr_discr_set_domain);
				 break;
			 }
			 break;
		 case 'm':
			 if ( !strcmp(key, "mode") ) {
				 /* n = 1; type = i:  UNUR_DISTR *distribution, int mode */
				 result = _unur_str_distr_set_i(distr,key,type_args,args,unur_distr_discr_set_mode);
				 break;
			 }
			 break;
		 case 'p':
			 if ( !strcmp(key, "pmf") ) {
				 /* n = 1; type = C:  UNUR_DISTR *distribution, const char *pmfstr */
				 result = _unur_str_distr_set_C(distr,key,type_args,args,unur_distr_discr_set_pmfstr);
				 break;
			 }
			 if ( !strcmp(key, "pmfparams") ) {
				 /* n = 2; type = Di:  UNUR_DISTR *distribution, const double *params, int n_params */
				 result = _unur_str_distr_set_Di(distr,key,type_args,args,unur_distr_discr_set_pmfparams);
				 break;
			 }
			 if ( !strcmp(key, "pmfstr") ) {
				 /* n = 1; type = C:  UNUR_DISTR *distribution, const char *pmfstr */
				 result = _unur_str_distr_set_C(distr,key,type_args,args,unur_distr_discr_set_pmfstr);
				 break;
			 }
			 if ( !strcmp(key, "pmfsum") ) {
				 /* n = 1; type = d:  UNUR_DISTR *distribution, double sum */
				 result = _unur_str_distr_set_d(distr,key,type_args,args,unur_distr_discr_set_pmfsum);
				 break;
			 }
			 if ( !strcmp(key, "pv") ) {
				 /* n = 2; type = Di:  UNUR_DISTR *distribution, const double *pv, int n_pv */
				 result = _unur_str_distr_set_Di(distr,key,type_args,args,unur_distr_discr_set_pv);
				 break;
			 }
		 }
		 break;
	 }

	 /* set calls for all distribution types */
	 if (result == UNUR_ERR_STR_UNKNOWN) {
		 switch (*key) {
		 case 'n':
			 if ( !strcmp(key, "name") ) {
				 /* n = 1; type = C:  UNUR_DISTR *distribution, const char *name */
				 result = _unur_str_distr_set_C(distr,key,type_args,args,unur_distr_set_name);
				 break;
			 }
		 }
	 }


	 /* Ignored set commands: */
 	 /* int unur_distr_cont_set_pdf( UNUR_DISTR *distribution, UNUR_FUNCT_CONT *pdf );
		 n = 1; type = ?	 */
	 /* int unur_distr_cont_set_dpdf( UNUR_DISTR *distribution, UNUR_FUNCT_CONT *dpdf );
		 n = 1; type = ?	 */
	 /* int unur_distr_cont_set_cdf( UNUR_DISTR *distribution, UNUR_FUNCT_CONT *cdf );
		 n = 1; type = ?	 */
	 /* int unur_distr_cont_set_logpdf( UNUR_DISTR *distribution, UNUR_FUNCT_CONT *logpdf );
		 n = 1; type = ?	 */
	 /* int unur_distr_cont_set_dlogpdf( UNUR_DISTR *distribution, UNUR_FUNCT_CONT *dlogpdf );
		 n = 1; type = ?	 */
	 /* int unur_distr_cont_set_logcdf( UNUR_DISTR *distribution, UNUR_FUNCT_CONT *logcdf );
		 n = 1; type = ?	 */
	 /* int unur_distr_cont_set_hr( UNUR_DISTR *distribution, UNUR_FUNCT_CONT *hazard );
		 n = 1; type = ?	 */
	 /* int unur_distr_discr_set_pmf( UNUR_DISTR *distribution, UNUR_FUNCT_DISCR *pmf );
		 n = 1; type = ?	 */
	 /* int unur_distr_discr_set_cdf( UNUR_DISTR *distribution, UNUR_FUNCT_DISCR *cdf );
		 n = 1; type = ?	 */
	 /* int unur_distr_set_extobj( UNUR_DISTR *distribution, const void *extobj );
		 n = 1; type = ?	 */


	 /* Subsituted set commands: */
 	 /* int unur_distr_cont_set_pdfstr( UNUR_DISTR *distribution, const char *pdfstr );
		 n = 1; type = C	 */
	 /* int unur_distr_cont_set_cdfstr( UNUR_DISTR *distribution, const char *cdfstr );
		 n = 1; type = C	 */
	 /* int unur_distr_cont_set_logpdfstr( UNUR_DISTR *distribution, const char *logpdfstr );
		 n = 1; type = C	 */
	 /* int unur_distr_cont_set_logcdfstr( UNUR_DISTR *distribution, const char *logcdfstr );
		 n = 1; type = C	 */
	 /* int unur_distr_cont_set_hrstr( UNUR_DISTR *distribution, const char *hrstr );
		 n = 1; type = C	 */
	 /* int unur_distr_discr_set_pmfstr( UNUR_DISTR *distribution, const char *pmfstr );
		 n = 1; type = C	 */
	 /* int unur_distr_discr_set_cdfstr( UNUR_DISTR *distribution, const char *cdfstr );
		 n = 1; type = C	 */


	 /* Unsupported set commands: */
 	/* int unur_distr_cemp_set_hist( UNUR_DISTR *distribution, const double *prob, int n_prob, double xmin, double xmax );
		 n = 4; type = Didd	 */
	/* int unur_distr_cont_set_invcdf( UNUR_DISTR *distribution, UNUR_FUNCT_CONT *invcdf );
		 n = 1; type = ?	 */
	/* int unur_distr_cont_set_pdfparams_vec( UNUR_DISTR *distribution, int par, const double *param_vec, int n_param_vec );
		 n = 3; type = iDi	 */
	/* int unur_distr_discr_set_invcdf( UNUR_DISTR *distribution, UNUR_IFUNCT_DISCR *invcdf );
		 n = 1; type = ?	 */


   /* special keyword */
   if (result == UNUR_ERR_STR_UNKNOWN)
     if (distr->type == UNUR_DISTR_CONT) {
       if ( !strcmp(key, "orderstatistics") ) {
	 /* make order statistics and replace distribution object */
	 *ptr_distr = _unur_str_distr_make_os (distr, key, type_args, args);
	 result = (*ptr_distr == NULL) ? UNUR_ERR_STR_SYNTAX : UNUR_SUCCESS;
       }
     }
  
  if (result == UNUR_ERR_STR_UNKNOWN) {
    /* unknown parameter */
    _unur_error_unknown(key,"parameter for given distribution");
    return UNUR_ERR_STR_UNKNOWN;
  }
  else if (result != UNUR_SUCCESS) {
    /* invalid data for parameter */
    _unur_error_invalid(key,"set call");
    return UNUR_ERR_STR_SYNTAX;
  }
  
  return UNUR_SUCCESS;
  
} /* end of _unur_str_distr_set() */

/*---------------------------------------------------------------------------*/

/*****************************************************************************/
/**  Methods                                                                **/
/*****************************************************************************/

/*---------------------------------------------------------------------------*/

struct unur_par *
_unur_str_par_new( const char *method, const UNUR_DISTR *distr )
     /*----------------------------------------------------------------------*/
     /* get new parameter object for method                                  */
     /*                                                                      */
     /* parameters:                                                          */
     /*   method ... string that contains method name                        */
     /*   distr  ... pointer to distribution object                          */
     /*                                                                      */
     /* return:                                                              */
     /*   default parameters (pointer to structure)                          */
     /*                                                                      */
     /* error:                                                               */
     /*   return NULL                                                        */
     /*----------------------------------------------------------------------*/
{
  struct unur_par *par = NULL;
  char method_unknown;
  
#ifdef UNUR_ENABLE_LOGGING
  /* write info into LOG file */
  if (_unur_default_debugflag & UNUR_DEBUG_SETUP)
    _unur_str_debug_string(1,"method",method);
#endif

  /* mark method as unknown (this is a very ugly hack) */
  par = (struct unur_par *) &method_unknown;

	 switch (*method) {
	 case 'a':
		 if ( !strcmp( method, "arou") ) {
			 par = unur_arou_new(distr);
			 break;
		 }
		 if ( !strcmp( method, "ars") ) {
			 par = unur_ars_new(distr);
			 break;
		 }
		 if ( !strcmp( method, "auto") ) {
			 par = unur_auto_new(distr);
			 break;
		 }
		 break;
	 case 'c':
		 if ( !strcmp( method, "cstd") ) {
			 par = unur_cstd_new(distr);
			 break;
		 }
		 break;
	 case 'd':
		 if ( !strcmp( method, "dari") ) {
			 par = unur_dari_new(distr);
			 break;
		 }
		 if ( !strcmp( method, "dau") ) {
			 par = unur_dau_new(distr);
			 break;
		 }
		 if ( !strcmp( method, "dgt") ) {
			 par = unur_dgt_new(distr);
			 break;
		 }
		 if ( !strcmp( method, "dsrou") ) {
			 par = unur_dsrou_new(distr);
			 break;
		 }
		 if ( !strcmp( method, "dss") ) {
			 par = unur_dss_new(distr);
			 break;
		 }
		 if ( !strcmp( method, "dstd") ) {
			 par = unur_dstd_new(distr);
			 break;
		 }
		 break;
	 case 'e':
		 if ( !strcmp( method, "empk") ) {
			 par = unur_empk_new(distr);
			 break;
		 }
		 if ( !strcmp( method, "empl") ) {
			 par = unur_empl_new(distr);
			 break;
		 }
		 break;
	 case 'g':
		 if ( !strcmp( method, "gibbs") ) {
			 par = unur_gibbs_new(distr);
			 break;
		 }
		 break;
	 case 'h':
		 if ( !strcmp( method, "hinv") ) {
			 par = unur_hinv_new(distr);
			 break;
		 }
		 if ( !strcmp( method, "hist") ) {
			 par = unur_hist_new(distr);
			 break;
		 }
		 if ( !strcmp( method, "hitro") ) {
			 par = unur_hitro_new(distr);
			 break;
		 }
		 if ( !strcmp( method, "hrb") ) {
			 par = unur_hrb_new(distr);
			 break;
		 }
		 if ( !strcmp( method, "hrd") ) {
			 par = unur_hrd_new(distr);
			 break;
		 }
		 if ( !strcmp( method, "hri") ) {
			 par = unur_hri_new(distr);
			 break;
		 }
		 break;
	 case 'i':
		 if ( !strcmp( method, "itdr") ) {
			 par = unur_itdr_new(distr);
			 break;
		 }
		 break;
	 case 'm':
		 if ( !strcmp( method, "mcorr") ) {
			 par = unur_mcorr_new(distr);
			 break;
		 }
		 if ( !strcmp( method, "mvstd") ) {
			 par = unur_mvstd_new(distr);
			 break;
		 }
		 if ( !strcmp( method, "mvtdr") ) {
			 par = unur_mvtdr_new(distr);
			 break;
		 }
		 break;
	 case 'n':
		 if ( !strcmp( method, "ninv") ) {
			 par = unur_ninv_new(distr);
			 break;
		 }
		 if ( !strcmp( method, "norta") ) {
			 par = unur_norta_new(distr);
			 break;
		 }
		 if ( !strcmp( method, "nrou") ) {
			 par = unur_nrou_new(distr);
			 break;
		 }
		 break;
	 case 'p':
		 if ( !strcmp( method, "pinv") ) {
			 par = unur_pinv_new(distr);
			 break;
		 }
		 break;
	 case 's':
		 if ( !strcmp( method, "srou") ) {
			 par = unur_srou_new(distr);
			 break;
		 }
		 if ( !strcmp( method, "ssr") ) {
			 par = unur_ssr_new(distr);
			 break;
		 }
		 break;
	 case 't':
		 if ( !strcmp( method, "tabl") ) {
			 par = unur_tabl_new(distr);
			 break;
		 }
		 if ( !strcmp( method, "tdr") ) {
			 par = unur_tdr_new(distr);
			 break;
		 }
		 break;
	 case 'u':
		 if ( !strcmp( method, "unif") ) {
			 par = unur_unif_new(distr);
			 break;
		 }
		 if ( !strcmp( method, "utdr") ) {
			 par = unur_utdr_new(distr);
			 break;
		 }
		 break;
	 case 'v':
		 if ( !strcmp( method, "vempk") ) {
			 par = unur_vempk_new(distr);
			 break;
		 }
		 if ( !strcmp( method, "vnrou") ) {
			 par = unur_vnrou_new(distr);
			 break;
		 }
	 }

   if (par == (struct unur_par *) &method_unknown) {
     /* unknown method */
     _unur_error_unknown(method,"method");
     par = NULL;
   }
   else if (par == NULL) {
     /* invalid data for chosen method */
     _unur_error_invalid(method,"method");
   }

   return par;
} /* end of _unur_str_par_new() */

/*---------------------------------------------------------------------------*/

int
_unur_str_par_set( UNUR_PAR *par, const char *key, char *value, struct unur_slist *mlist )
     /*----------------------------------------------------------------------*/
     /* set parameters for method                                            */
     /*                                                                      */
     /* parameters:                                                          */
     /*   par   ... pointer to parameter for building generator object       */
     /*   key   ... string that contains key for parameter                   */
     /*   value ... string that contains list of arguments for parameter     */
     /*   mlist ... list of allocated memory blocks                          */
     /*                                                                      */
     /* return:                                                              */
     /*   UNUR_SUCCESS ... on success                                        */
     /*   error code   ... on error                                          */
     /*----------------------------------------------------------------------*/
{
  int result = 0;                  /* result of UNU.RAN set call (0 or 1) */

  /* storing arguments of set calls: */
  char type_args[MAX_SET_ARGS+1];  /* array containing info about arguments */
  char *args[MAX_SET_ARGS+1];      /* pointers to argument strings */

  /* tokenize argument string */
  if (_unur_str_set_args( value, type_args, args, MAX_SET_ARGS ) < 0) {
    /* error */
    return UNUR_ERR_STR_SYNTAX;
  }

  /* set result indicator to unknown */
  result = UNUR_ERR_STR_UNKNOWN;

  /* find and execute set call */
	 switch (par->method) {
	 case UNUR_METH_AROU:
		 switch (*key) {
		 case 'c':
			 if ( !strcmp(key, "cpoints") ) {
				 /* n = 2; type = iD:  UNUR_PAR *parameters, int n_stp, const double *stp */
				 result = _unur_str_par_set_iD(par,key,type_args,args,unur_arou_set_cpoints,mlist);
				 break;
			 }
			 break;
		 case 'd':
			 if ( !strcmp(key, "darsfactor") ) {
				 /* n = 1; type = d:  UNUR_PAR *parameters, double factor */
				 result = _unur_str_par_set_d(par,key,type_args,args,unur_arou_set_darsfactor);
				 break;
			 }
			 break;
		 case 'g':
			 if ( !strcmp(key, "guidefactor") ) {
				 /* n = 1; type = d:  UNUR_PAR *parameters, double factor */
				 result = _unur_str_par_set_d(par,key,type_args,args,unur_arou_set_guidefactor);
				 break;
			 }
			 break;
		 case 'm':
			 if ( !strcmp(key, "max_segments") ) {
				 /* n = 1; type = i:  UNUR_PAR *parameters, int max_segs */
				 result = _unur_str_par_set_i(par,key,type_args,args,unur_arou_set_max_segments);
				 break;
			 }
			 if ( !strcmp(key, "max_sqhratio") ) {
				 /* n = 1; type = d:  UNUR_PAR *parameters, double max_ratio */
				 result = _unur_str_par_set_d(par,key,type_args,args,unur_arou_set_max_sqhratio);
				 break;
			 }
			 break;
		 case 'p':
			 if ( !strcmp(key, "pedantic") ) {
				 /* n = 1; type = i:  UNUR_PAR *parameters, int pedantic */
				 result = _unur_str_par_set_i(par,key,type_args,args,unur_arou_set_pedantic);
				 break;
			 }
			 break;
		 case 'u':
			 if ( !strcmp(key, "usecenter") ) {
				 /* n = 1; type = i:  UNUR_PAR *parameters, int usecenter */
				 result = _unur_str_par_set_i(par,key,type_args,args,unur_arou_set_usecenter);
				 break;
			 }
			 if ( !strcmp(key, "usedars") ) {
				 /* n = 1; type = i:  UNUR_PAR *parameters, int usedars */
				 result = _unur_str_par_set_i(par,key,type_args,args,unur_arou_set_usedars);
				 break;
			 }
			 break;
		 case 'v':
			 if ( !strcmp(key, "verify") ) {
				 /* n = 1; type = i:  UNUR_PAR *parameters, int verify */
				 result = _unur_str_par_set_i(par,key,type_args,args,unur_arou_set_verify);
				 break;
			 }
		 }
		 break;
	 case UNUR_METH_ARS:
		 switch (*key) {
		 case 'c':
			 if ( !strcmp(key, "cpoints") ) {
				 /* n = 2; type = iD:  UNUR_PAR *parameters, int n_cpoints, const double *cpoints */
				 result = _unur_str_par_set_iD(par,key,type_args,args,unur_ars_set_cpoints,mlist);
				 break;
			 }
			 break;
		 case 'm':
			 if ( !strcmp(key, "max_intervals") ) {
				 /* n = 1; type = i:  UNUR_PAR *parameters, int max_ivs */
				 result = _unur_str_par_set_i(par,key,type_args,args,unur_ars_set_max_intervals);
				 break;
			 }
			 if ( !strcmp(key, "max_iter") ) {
				 /* n = 1; type = i:  UNUR_PAR *parameters, int max_iter */
				 result = _unur_str_par_set_i(par,key,type_args,args,unur_ars_set_max_iter);
				 break;
			 }
			 break;
		 case 'p':
			 if ( !strcmp(key, "pedantic") ) {
				 /* n = 1; type = i:  UNUR_PAR *parameters, int pedantic */
				 result = _unur_str_par_set_i(par,key,type_args,args,unur_ars_set_pedantic);
				 break;
			 }
			 break;
		 case 'r':
			 if ( !strcmp(key, "reinit_ncpoints") ) {
				 /* n = 1; type = i:  UNUR_PAR *parameters, int ncpoints */
				 result = _unur_str_par_set_i(par,key,type_args,args,unur_ars_set_reinit_ncpoints);
				 break;
			 }
			 if ( !strcmp(key, "reinit_percentiles") ) {
				 /* n = 2; type = iD:  UNUR_PAR *parameters, int n_percentiles, const double *percentiles */
				 result = _unur_str_par_set_iD(par,key,type_args,args,unur_ars_set_reinit_percentiles,mlist);
				 break;
			 }
			 break;
		 case 'v':
			 if ( !strcmp(key, "verify") ) {
				 /* n = 1; type = i:  UNUR_PAR *parameters, int verify */
				 result = _unur_str_par_set_i(par,key,type_args,args,unur_ars_set_verify);
				 break;
			 }
		 }
		 break;
	 case UNUR_METH_AUTO:
		 switch (*key) {
		 case 'l':
			 if ( !strcmp(key, "logss") ) {
				 /* n = 1; type = i:  UNUR_PAR *parameters, int logss */
				 result = _unur_str_par_set_i(par,key,type_args,args,unur_auto_set_logss);
				 break;
			 }
		 }
		 break;
	 case UNUR_METH_CSTD:
		 switch (*key) {
		 case 'v':
			 if ( !strcmp(key, "variant") ) {
				 /* n = 1; type = u:  UNUR_PAR *parameters, unsigned variant */
				 result = _unur_str_par_set_u(par,key,type_args,args,unur_cstd_set_variant);
				 break;
			 }
		 }
		 break;
	 case UNUR_METH_DARI:
		 switch (*key) {
		 case 'c':
			 if ( !strcmp(key, "cpfactor") ) {
				 /* n = 1; type = d:  UNUR_PAR *parameters, double cp_factor */
				 result = _unur_str_par_set_d(par,key,type_args,args,unur_dari_set_cpfactor);
				 break;
			 }
			 break;
		 case 's':
			 if ( !strcmp(key, "squeeze") ) {
				 /* n = 1; type = i:  UNUR_PAR *parameters, int squeeze */
				 result = _unur_str_par_set_i(par,key,type_args,args,unur_dari_set_squeeze);
				 break;
			 }
			 break;
		 case 't':
			 if ( !strcmp(key, "tablesize") ) {
				 /* n = 1; type = i:  UNUR_PAR *parameters, int size */
				 result = _unur_str_par_set_i(par,key,type_args,args,unur_dari_set_tablesize);
				 break;
			 }
			 break;
		 case 'v':
			 if ( !strcmp(key, "verify") ) {
				 /* n = 1; type = i:  UNUR_PAR *parameters, int verify */
				 result = _unur_str_par_set_i(par,key,type_args,args,unur_dari_set_verify);
				 break;
			 }
		 }
		 break;
	 case UNUR_METH_DAU:
		 switch (*key) {
		 case 'u':
			 if ( !strcmp(key, "urnfactor") ) {
				 /* n = 1; type = d:  UNUR_PAR *parameters, double factor */
				 result = _unur_str_par_set_d(par,key,type_args,args,unur_dau_set_urnfactor);
				 break;
			 }
		 }
		 break;
	 case UNUR_METH_DGT:
		 switch (*key) {
		 case 'g':
			 if ( !strcmp(key, "guidefactor") ) {
				 /* n = 1; type = d:  UNUR_PAR *parameters, double factor */
				 result = _unur_str_par_set_d(par,key,type_args,args,unur_dgt_set_guidefactor);
				 break;
			 }
			 break;
		 case 'v':
			 if ( !strcmp(key, "variant") ) {
				 /* n = 1; type = u:  UNUR_PAR *parameters, unsigned variant */
				 result = _unur_str_par_set_u(par,key,type_args,args,unur_dgt_set_variant);
				 break;
			 }
		 }
		 break;
	 case UNUR_METH_DSROU:
		 switch (*key) {
		 case 'c':
			 if ( !strcmp(key, "cdfatmode") ) {
				 /* n = 1; type = d:  UNUR_PAR *parameters, double Fmode */
				 result = _unur_str_par_set_d(par,key,type_args,args,unur_dsrou_set_cdfatmode);
				 break;
			 }
			 break;
		 case 'v':
			 if ( !strcmp(key, "verify") ) {
				 /* n = 1; type = i:  UNUR_PAR *parameters, int verify */
				 result = _unur_str_par_set_i(par,key,type_args,args,unur_dsrou_set_verify);
				 break;
			 }
		 }
		 break;
	 case UNUR_METH_DSTD:
		 switch (*key) {
		 case 'v':
			 if ( !strcmp(key, "variant") ) {
				 /* n = 1; type = u:  UNUR_PAR *parameters, unsigned variant */
				 result = _unur_str_par_set_u(par,key,type_args,args,unur_dstd_set_variant);
				 break;
			 }
		 }
		 break;
	 case UNUR_METH_EMPK:
		 switch (*key) {
		 case 'b':
			 if ( !strcmp(key, "beta") ) {
				 /* n = 1; type = d:  UNUR_PAR *parameters, double beta */
				 result = _unur_str_par_set_d(par,key,type_args,args,unur_empk_set_beta);
				 break;
			 }
			 break;
		 case 'k':
			 if ( !strcmp(key, "kernel") ) {
				 /* n = 1; type = u:  UNUR_PAR *parameters, unsigned kernel*/
				 result = _unur_str_par_set_u(par,key,type_args,args,unur_empk_set_kernel);
				 break;
			 }
			 break;
		 case 'p':
			 if ( !strcmp(key, "positive") ) {
				 /* n = 1; type = i:  UNUR_PAR *parameters, int positive */
				 result = _unur_str_par_set_i(par,key,type_args,args,unur_empk_set_positive);
				 break;
			 }
			 break;
		 case 's':
			 if ( !strcmp(key, "smoothing") ) {
				 /* n = 1; type = d:  UNUR_PAR *parameters, double smoothing */
				 result = _unur_str_par_set_d(par,key,type_args,args,unur_empk_set_smoothing);
				 break;
			 }
			 break;
		 case 'v':
			 if ( !strcmp(key, "varcor") ) {
				 /* n = 1; type = i:  UNUR_PAR *parameters, int varcor */
				 result = _unur_str_par_set_i(par,key,type_args,args,unur_empk_set_varcor);
				 break;
			 }
		 }
		 break;
	 case UNUR_METH_GIBBS:
		 switch (*key) {
		 case 'b':
			 if ( !strcmp(key, "burnin") ) {
				 /* n = 1; type = i:  UNUR_PAR *parameters, int burnin */
				 result = _unur_str_par_set_i(par,key,type_args,args,unur_gibbs_set_burnin);
				 break;
			 }
			 break;
		 case 'c':
			 if ( !strcmp(key, "c") ) {
				 /* n = 1; type = d:  UNUR_PAR *parameters, double c */
				 result = _unur_str_par_set_d(par,key,type_args,args,unur_gibbs_set_c);
				 break;
			 }
			 break;
		 case 't':
			 if ( !strcmp(key, "thinning") ) {
				 /* n = 1; type = i:  UNUR_PAR *parameters, int thinning */
				 result = _unur_str_par_set_i(par,key,type_args,args,unur_gibbs_set_thinning);
				 break;
			 }
			 break;
		 case 'v':
			 if ( !strcmp(key, "variant_coordinate") ) {
				 /* n = 0; type = :  UNUR_PAR *parameters */
				 result = _unur_str_par_set_void(par,key,type_args,args,unur_gibbs_set_variant_coordinate);
				 break;
			 }
			 if ( !strcmp(key, "variant_random_direction") ) {
				 /* n = 0; type = :  UNUR_PAR *parameters */
				 result = _unur_str_par_set_void(par,key,type_args,args,unur_gibbs_set_variant_random_direction);
				 break;
			 }
		 }
		 break;
	 case UNUR_METH_HINV:
		 switch (*key) {
		 case 'b':
			 if ( !strcmp(key, "boundary") ) {
				 /* n = 2; type = dd:  UNUR_PAR *parameters, double left, double right */
				 result = _unur_str_par_set_dd(par,key,type_args,args,unur_hinv_set_boundary);
				 break;
			 }
			 break;
		 case 'c':
			 if ( !strcmp(key, "cpoints") ) {
				 /* n = 2; type = Di:  UNUR_PAR *parameters, const double *stp, int n_stp */
				 result = _unur_str_par_set_Di(par,key,type_args,args,unur_hinv_set_cpoints,mlist);
				 break;
			 }
			 break;
		 case 'g':
			 if ( !strcmp(key, "guidefactor") ) {
				 /* n = 1; type = d:  UNUR_PAR *parameters, double factor */
				 result = _unur_str_par_set_d(par,key,type_args,args,unur_hinv_set_guidefactor);
				 break;
			 }
			 break;
		 case 'm':
			 if ( !strcmp(key, "max_intervals") ) {
				 /* n = 1; type = i:  UNUR_PAR *parameters, int max_ivs */
				 result = _unur_str_par_set_i(par,key,type_args,args,unur_hinv_set_max_intervals);
				 break;
			 }
			 break;
		 case 'o':
			 if ( !strcmp(key, "order") ) {
				 /* n = 1; type = i:  UNUR_PAR *parameters, int order*/
				 result = _unur_str_par_set_i(par,key,type_args,args,unur_hinv_set_order);
				 break;
			 }
			 break;
		 case 'u':
			 if ( !strcmp(key, "u_resolution") ) {
				 /* n = 1; type = d:  UNUR_PAR *parameters, double u_resolution*/
				 result = _unur_str_par_set_d(par,key,type_args,args,unur_hinv_set_u_resolution);
				 break;
			 }
		 }
		 break;
	 case UNUR_METH_HITRO:
		 switch (*key) {
		 case 'a':
			 if ( !strcmp(key, "adaptive_multiplier") ) {
				 /* n = 1; type = d:  UNUR_PAR *parameters, double factor */
				 result = _unur_str_par_set_d(par,key,type_args,args,unur_hitro_set_adaptive_multiplier);
				 break;
			 }
			 break;
		 case 'b':
			 if ( !strcmp(key, "burnin") ) {
				 /* n = 1; type = i:  UNUR_PAR *parameters, int burnin */
				 result = _unur_str_par_set_i(par,key,type_args,args,unur_hitro_set_burnin);
				 break;
			 }
			 break;
		 case 'r':
			 if ( !strcmp(key, "r") ) {
				 /* n = 1; type = d:  UNUR_PAR *parameters, double r */
				 result = _unur_str_par_set_d(par,key,type_args,args,unur_hitro_set_r);
				 break;
			 }
			 break;
		 case 't':
			 if ( !strcmp(key, "thinning") ) {
				 /* n = 1; type = i:  UNUR_PAR *parameters, int thinning */
				 result = _unur_str_par_set_i(par,key,type_args,args,unur_hitro_set_thinning);
				 break;
			 }
			 break;
		 case 'u':
			 if ( !strcmp(key, "use_adaptiveline") ) {
				 /* n = 1; type = i:  UNUR_PAR *parameters, int adaptive */
				 result = _unur_str_par_set_i(par,key,type_args,args,unur_hitro_set_use_adaptiveline);
				 break;
			 }
			 if ( !strcmp(key, "use_adaptiverectangle") ) {
				 /* n = 1; type = i:  UNUR_PAR *parameters, int adaptive */
				 result = _unur_str_par_set_i(par,key,type_args,args,unur_hitro_set_use_adaptiverectangle);
				 break;
			 }
			 if ( !strcmp(key, "use_boundingrectangle") ) {
				 /* n = 1; type = i:  UNUR_PAR *parameters, int rectangle */
				 result = _unur_str_par_set_i(par,key,type_args,args,unur_hitro_set_use_boundingrectangle);
				 break;
			 }
			 break;
		 case 'v':
			 if ( !strcmp(key, "v") ) {
				 /* n = 1; type = d:  UNUR_PAR *parameters, double vmax */
				 result = _unur_str_par_set_d(par,key,type_args,args,unur_hitro_set_v);
				 break;
			 }
			 if ( !strcmp(key, "variant_coordinate") ) {
				 /* n = 0; type = :  UNUR_PAR *parameters */
				 result = _unur_str_par_set_void(par,key,type_args,args,unur_hitro_set_variant_coordinate);
				 break;
			 }
			 if ( !strcmp(key, "variant_random_direction") ) {
				 /* n = 0; type = :  UNUR_PAR *parameters */
				 result = _unur_str_par_set_void(par,key,type_args,args,unur_hitro_set_variant_random_direction);
				 break;
			 }
		 }
		 break;
	 case UNUR_METH_HRB:
		 switch (*key) {
		 case 'u':
			 if ( !strcmp(key, "upperbound") ) {
				 /* n = 1; type = d:  UNUR_PAR *parameters, double upperbound */
				 result = _unur_str_par_set_d(par,key,type_args,args,unur_hrb_set_upperbound);
				 break;
			 }
			 break;
		 case 'v':
			 if ( !strcmp(key, "verify") ) {
				 /* n = 1; type = i:  UNUR_PAR *parameters, int verify */
				 result = _unur_str_par_set_i(par,key,type_args,args,unur_hrb_set_verify);
				 break;
			 }
		 }
		 break;
	 case UNUR_METH_HRD:
		 switch (*key) {
		 case 'v':
			 if ( !strcmp(key, "verify") ) {
				 /* n = 1; type = i:  UNUR_PAR *parameters, int verify */
				 result = _unur_str_par_set_i(par,key,type_args,args,unur_hrd_set_verify);
				 break;
			 }
		 }
		 break;
	 case UNUR_METH_HRI:
		 switch (*key) {
		 case 'p':
			 if ( !strcmp(key, "p0") ) {
				 /* n = 1; type = d:  UNUR_PAR *parameters, double p0 */
				 result = _unur_str_par_set_d(par,key,type_args,args,unur_hri_set_p0);
				 break;
			 }
			 break;
		 case 'v':
			 if ( !strcmp(key, "verify") ) {
				 /* n = 1; type = i:  UNUR_PAR *parameters, int verify */
				 result = _unur_str_par_set_i(par,key,type_args,args,unur_hri_set_verify);
				 break;
			 }
		 }
		 break;
	 case UNUR_METH_ITDR:
		 switch (*key) {
		 case 'c':
			 if ( !strcmp(key, "cp") ) {
				 /* n = 1; type = d:  UNUR_PAR *parameters, double cp */
				 result = _unur_str_par_set_d(par,key,type_args,args,unur_itdr_set_cp);
				 break;
			 }
			 if ( !strcmp(key, "ct") ) {
				 /* n = 1; type = d:  UNUR_PAR *parameters, double ct */
				 result = _unur_str_par_set_d(par,key,type_args,args,unur_itdr_set_ct);
				 break;
			 }
			 break;
		 case 'v':
			 if ( !strcmp(key, "verify") ) {
				 /* n = 1; type = i:  UNUR_PAR *parameters, int verify */
				 result = _unur_str_par_set_i(par,key,type_args,args,unur_itdr_set_verify);
				 break;
			 }
			 break;
		 case 'x':
			 if ( !strcmp(key, "xi") ) {
				 /* n = 1; type = d:  UNUR_PAR *parameters, double xi */
				 result = _unur_str_par_set_d(par,key,type_args,args,unur_itdr_set_xi);
				 break;
			 }
		 }
		 break;
	 case UNUR_METH_MVTDR:
		 switch (*key) {
		 case 'b':
			 if ( !strcmp(key, "boundsplitting") ) {
				 /* n = 1; type = d:  UNUR_PAR *parameters, double boundsplitting */
				 result = _unur_str_par_set_d(par,key,type_args,args,unur_mvtdr_set_boundsplitting);
				 break;
			 }
			 break;
		 case 'm':
			 if ( !strcmp(key, "maxcones") ) {
				 /* n = 1; type = i:  UNUR_PAR *parameters, int maxcones */
				 result = _unur_str_par_set_i(par,key,type_args,args,unur_mvtdr_set_maxcones);
				 break;
			 }
			 break;
		 case 's':
			 if ( !strcmp(key, "stepsmin") ) {
				 /* n = 1; type = i:  UNUR_PAR *parameters, int stepsmin */
				 result = _unur_str_par_set_i(par,key,type_args,args,unur_mvtdr_set_stepsmin);
				 break;
			 }
			 break;
		 case 'v':
			 if ( !strcmp(key, "verify") ) {
				 /* n = 1; type = i:  UNUR_PAR *parameters, int verify */
				 result = _unur_str_par_set_i(par,key,type_args,args,unur_mvtdr_set_verify);
				 break;
			 }
		 }
		 break;
	 case UNUR_METH_NINV:
		 switch (*key) {
		 case 'm':
			 if ( !strcmp(key, "max_iter") ) {
				 /* n = 1; type = i:  UNUR_PAR *parameters, int max_iter */
				 result = _unur_str_par_set_i(par,key,type_args,args,unur_ninv_set_max_iter);
				 break;
			 }
			 break;
		 case 's':
			 if ( !strcmp(key, "start") ) {
				 /* n = 2; type = dd:  UNUR_PAR *parameters, double left, double right*/
				 result = _unur_str_par_set_dd(par,key,type_args,args,unur_ninv_set_start);
				 break;
			 }
			 break;
		 case 't':
			 if ( !strcmp(key, "table") ) {
				 /* n = 1; type = i: UNUR_PAR *parameters, int no_of_points*/
				 result = _unur_str_par_set_i(par,key,type_args,args,unur_ninv_set_table);
				 break;
			 }
			 break;
		 case 'u':
			 if ( !strcmp(key, "u_resolution") ) {
				 /* n = 1; type = d:  UNUR_PAR *parameters, double u_resolution*/
				 result = _unur_str_par_set_d(par,key,type_args,args,unur_ninv_set_u_resolution);
				 break;
			 }
			 if ( !strcmp(key, "usebisect") ) {
				 /* n = 0; type = :  UNUR_PAR *parameters */
				 result = _unur_str_par_set_void(par,key,type_args,args,unur_ninv_set_usebisect);
				 break;
			 }
			 if ( !strcmp(key, "usenewton") ) {
				 /* n = 0; type = :  UNUR_PAR *parameters */
				 result = _unur_str_par_set_void(par,key,type_args,args,unur_ninv_set_usenewton);
				 break;
			 }
			 if ( !strcmp(key, "useregula") ) {
				 /* n = 0; type = :  UNUR_PAR *parameters */
				 result = _unur_str_par_set_void(par,key,type_args,args,unur_ninv_set_useregula);
				 break;
			 }
			 break;
		 case 'x':
			 if ( !strcmp(key, "x_resolution") ) {
				 /* n = 1; type = d:  UNUR_PAR *parameters, double x_resolution*/
				 result = _unur_str_par_set_d(par,key,type_args,args,unur_ninv_set_x_resolution);
				 break;
			 }
		 }
		 break;
	 case UNUR_METH_NROU:
		 switch (*key) {
		 case 'c':
			 if ( !strcmp(key, "center") ) {
				 /* n = 1; type = d:  UNUR_PAR *parameters, double center */
				 result = _unur_str_par_set_d(par,key,type_args,args,unur_nrou_set_center);
				 break;
			 }
			 break;
		 case 'r':
			 if ( !strcmp(key, "r") ) {
				 /* n = 1; type = d:  UNUR_PAR *parameters, double r */
				 result = _unur_str_par_set_d(par,key,type_args,args,unur_nrou_set_r);
				 break;
			 }
			 break;
		 case 'u':
			 if ( !strcmp(key, "u") ) {
				 /* n = 2; type = dd:  UNUR_PAR *parameters, double umin, double umax */
				 result = _unur_str_par_set_dd(par,key,type_args,args,unur_nrou_set_u);
				 break;
			 }
			 break;
		 case 'v':
			 if ( !strcmp(key, "v") ) {
				 /* n = 1; type = d:  UNUR_PAR *parameters, double vmax */
				 result = _unur_str_par_set_d(par,key,type_args,args,unur_nrou_set_v);
				 break;
			 }
			 if ( !strcmp(key, "verify") ) {
				 /* n = 1; type = i:  UNUR_PAR *parameters, int verify */
				 result = _unur_str_par_set_i(par,key,type_args,args,unur_nrou_set_verify);
				 break;
			 }
		 }
		 break;
	 case UNUR_METH_PINV:
		 switch (*key) {
		 case 'b':
			 if ( !strcmp(key, "boundary") ) {
				 /* n = 2; type = dd:  UNUR_PAR *parameters, double left, double right */
				 result = _unur_str_par_set_dd(par,key,type_args,args,unur_pinv_set_boundary);
				 break;
			 }
			 break;
		 case 'k':
			 if ( !strcmp(key, "keepcdf") ) {
				 /* n = 1; type = i:  UNUR_PAR *parameters, int keepcdf*/
				 result = _unur_str_par_set_i(par,key,type_args,args,unur_pinv_set_keepcdf);
				 break;
			 }
			 break;
		 case 'm':
			 if ( !strcmp(key, "max_intervals") ) {
				 /* n = 1; type = i:  UNUR_PAR *parameters, int max_ivs */
				 result = _unur_str_par_set_i(par,key,type_args,args,unur_pinv_set_max_intervals);
				 break;
			 }
			 break;
		 case 'o':
			 if ( !strcmp(key, "order") ) {
				 /* n = 1; type = i:  UNUR_PAR *parameters, int order*/
				 result = _unur_str_par_set_i(par,key,type_args,args,unur_pinv_set_order);
				 break;
			 }
			 break;
		 case 's':
			 if ( !strcmp(key, "searchboundary") ) {
				 /* n = 2; type = ii:  UNUR_PAR *parameters, int left, int right */
				 result = _unur_str_par_set_ii(par,key,type_args,args,unur_pinv_set_searchboundary);
				 break;
			 }
			 if ( !strcmp(key, "smoothness") ) {
				 /* n = 1; type = i:  UNUR_PAR *parameters, int smoothness*/
				 result = _unur_str_par_set_i(par,key,type_args,args,unur_pinv_set_smoothness);
				 break;
			 }
			 break;
		 case 'u':
			 if ( !strcmp(key, "u_resolution") ) {
				 /* n = 1; type = d:  UNUR_PAR *parameters, double u_resolution*/
				 result = _unur_str_par_set_d(par,key,type_args,args,unur_pinv_set_u_resolution);
				 break;
			 }
			 if ( !strcmp(key, "use_upoints") ) {
				 /* n = 1; type = i:  UNUR_PAR *parameters, int use_upoints */
				 result = _unur_str_par_set_i(par,key,type_args,args,unur_pinv_set_use_upoints);
				 break;
			 }
			 if ( !strcmp(key, "usecdf") ) {
				 /* n = 0; type = :  UNUR_PAR *parameters */
				 result = _unur_str_par_set_void(par,key,type_args,args,unur_pinv_set_usecdf);
				 break;
			 }
			 if ( !strcmp(key, "usepdf") ) {
				 /* n = 0; type = :  UNUR_PAR *parameters */
				 result = _unur_str_par_set_void(par,key,type_args,args,unur_pinv_set_usepdf);
				 break;
			 }
		 }
		 break;
	 case UNUR_METH_SROU:
		 switch (*key) {
		 case 'c':
			 if ( !strcmp(key, "cdfatmode") ) {
				 /* n = 1; type = d:  UNUR_PAR *parameters, double Fmode */
				 result = _unur_str_par_set_d(par,key,type_args,args,unur_srou_set_cdfatmode);
				 break;
			 }
			 break;
		 case 'p':
			 if ( !strcmp(key, "pdfatmode") ) {
				 /* n = 1; type = d:  UNUR_PAR *parameters, double fmode */
				 result = _unur_str_par_set_d(par,key,type_args,args,unur_srou_set_pdfatmode);
				 break;
			 }
			 break;
		 case 'r':
			 if ( !strcmp(key, "r") ) {
				 /* n = 1; type = d:  UNUR_PAR *parameters, double r */
				 result = _unur_str_par_set_d(par,key,type_args,args,unur_srou_set_r);
				 break;
			 }
			 break;
		 case 'u':
			 if ( !strcmp(key, "usemirror") ) {
				 /* n = 1; type = i:  UNUR_PAR *parameters, int usemirror */
				 result = _unur_str_par_set_i(par,key,type_args,args,unur_srou_set_usemirror);
				 break;
			 }
			 if ( !strcmp(key, "usesqueeze") ) {
				 /* n = 1; type = i:  UNUR_PAR *parameters, int usesqueeze */
				 result = _unur_str_par_set_i(par,key,type_args,args,unur_srou_set_usesqueeze);
				 break;
			 }
			 break;
		 case 'v':
			 if ( !strcmp(key, "verify") ) {
				 /* n = 1; type = i:  UNUR_PAR *parameters, int verify */
				 result = _unur_str_par_set_i(par,key,type_args,args,unur_srou_set_verify);
				 break;
			 }
		 }
		 break;
	 case UNUR_METH_SSR:
		 switch (*key) {
		 case 'c':
			 if ( !strcmp(key, "cdfatmode") ) {
				 /* n = 1; type = d:  UNUR_PAR *parameters, double Fmode */
				 result = _unur_str_par_set_d(par,key,type_args,args,unur_ssr_set_cdfatmode);
				 break;
			 }
			 break;
		 case 'p':
			 if ( !strcmp(key, "pdfatmode") ) {
				 /* n = 1; type = d:  UNUR_PAR *parameters, double fmode */
				 result = _unur_str_par_set_d(par,key,type_args,args,unur_ssr_set_pdfatmode);
				 break;
			 }
			 break;
		 case 'u':
			 if ( !strcmp(key, "usesqueeze") ) {
				 /* n = 1; type = i:  UNUR_PAR *parameters, int usesqueeze */
				 result = _unur_str_par_set_i(par,key,type_args,args,unur_ssr_set_usesqueeze);
				 break;
			 }
			 break;
		 case 'v':
			 if ( !strcmp(key, "verify") ) {
				 /* n = 1; type = i:  UNUR_PAR *parameters, int verify */
				 result = _unur_str_par_set_i(par,key,type_args,args,unur_ssr_set_verify);
				 break;
			 }
		 }
		 break;
	 case UNUR_METH_TABL:
		 switch (*key) {
		 case 'a':
			 if ( !strcmp(key, "areafraction") ) {
				 /* n = 1; type = d:  UNUR_PAR *parameters, double fraction */
				 result = _unur_str_par_set_d(par,key,type_args,args,unur_tabl_set_areafraction);
				 break;
			 }
			 break;
		 case 'b':
			 if ( !strcmp(key, "boundary") ) {
				 /* n = 2; type = dd:  UNUR_PAR *parameters, double left, double right */
				 result = _unur_str_par_set_dd(par,key,type_args,args,unur_tabl_set_boundary);
				 break;
			 }
			 break;
		 case 'c':
			 if ( !strcmp(key, "cpoints") ) {
				 /* n = 2; type = iD:  UNUR_PAR *parameters, int n_cpoints, const double *cpoints */
				 result = _unur_str_par_set_iD(par,key,type_args,args,unur_tabl_set_cpoints,mlist);
				 break;
			 }
			 break;
		 case 'd':
			 if ( !strcmp(key, "darsfactor") ) {
				 /* n = 1; type = d:  UNUR_PAR *parameters, double factor */
				 result = _unur_str_par_set_d(par,key,type_args,args,unur_tabl_set_darsfactor);
				 break;
			 }
			 break;
		 case 'g':
			 if ( !strcmp(key, "guidefactor") ) {
				 /* n = 1; type = d:  UNUR_PAR *parameters, double factor */
				 result = _unur_str_par_set_d(par,key,type_args,args,unur_tabl_set_guidefactor);
				 break;
			 }
			 break;
		 case 'm':
			 if ( !strcmp(key, "max_intervals") ) {
				 /* n = 1; type = i:  UNUR_PAR *parameters, int max_ivs */
				 result = _unur_str_par_set_i(par,key,type_args,args,unur_tabl_set_max_intervals);
				 break;
			 }
			 if ( !strcmp(key, "max_sqhratio") ) {
				 /* n = 1; type = d:  UNUR_PAR *parameters, double max_ratio */
				 result = _unur_str_par_set_d(par,key,type_args,args,unur_tabl_set_max_sqhratio);
				 break;
			 }
			 break;
		 case 'n':
			 if ( !strcmp(key, "nstp") ) {
				 /* n = 1; type = i:  UNUR_PAR *parameters, int n_stp */
				 result = _unur_str_par_set_i(par,key,type_args,args,unur_tabl_set_nstp);
				 break;
			 }
			 break;
		 case 'p':
			 if ( !strcmp(key, "pedantic") ) {
				 /* n = 1; type = i:  UNUR_PAR *parameters, int pedantic */
				 result = _unur_str_par_set_i(par,key,type_args,args,unur_tabl_set_pedantic);
				 break;
			 }
			 break;
		 case 's':
			 if ( !strcmp(key, "slopes") ) {
				 /* n = 2; type = Di:  UNUR_PAR *parameters, const double *slopes, int n_slopes */
				 result = _unur_str_par_set_Di(par,key,type_args,args,unur_tabl_set_slopes,mlist);
				 break;
			 }
			 break;
		 case 'u':
			 if ( !strcmp(key, "usedars") ) {
				 /* n = 1; type = i:  UNUR_PAR *parameters, int usedars */
				 result = _unur_str_par_set_i(par,key,type_args,args,unur_tabl_set_usedars);
				 break;
			 }
			 if ( !strcmp(key, "useear") ) {
				 /* n = 1; type = i:  UNUR_PAR *parameters, int useear */
				 result = _unur_str_par_set_i(par,key,type_args,args,unur_tabl_set_useear);
				 break;
			 }
			 break;
		 case 'v':
			 if ( !strcmp(key, "variant_ia") ) {
				 /* n = 1; type = i:  UNUR_PAR *parameters, int use_ia */
				 result = _unur_str_par_set_i(par,key,type_args,args,unur_tabl_set_variant_ia);
				 break;
			 }
			 if ( !strcmp(key, "variant_splitmode") ) {
				 /* n = 1; type = u:  UNUR_PAR *parameters, unsigned splitmode */
				 result = _unur_str_par_set_u(par,key,type_args,args,unur_tabl_set_variant_splitmode);
				 break;
			 }
			 if ( !strcmp(key, "verify") ) {
				 /* n = 1; type = i:  UNUR_PAR *parameters, int verify */
				 result = _unur_str_par_set_i(par,key,type_args,args,unur_tabl_set_verify);
				 break;
			 }
		 }
		 break;
	 case UNUR_METH_TDR:
		 switch (*key) {
		 case 'c':
			 if ( !strcmp(key, "c") ) {
				 /* n = 1; type = d:  UNUR_PAR *parameters, double c */
				 result = _unur_str_par_set_d(par,key,type_args,args,unur_tdr_set_c);
				 break;
			 }
			 if ( !strcmp(key, "cpoints") ) {
				 /* n = 2; type = iD:  UNUR_PAR *parameters, int n_stp, const double *stp */
				 result = _unur_str_par_set_iD(par,key,type_args,args,unur_tdr_set_cpoints,mlist);
				 break;
			 }
			 break;
		 case 'd':
			 if ( !strcmp(key, "darsfactor") ) {
				 /* n = 1; type = d:  UNUR_PAR *parameters, double factor */
				 result = _unur_str_par_set_d(par,key,type_args,args,unur_tdr_set_darsfactor);
				 break;
			 }
			 break;
		 case 'g':
			 if ( !strcmp(key, "guidefactor") ) {
				 /* n = 1; type = d:  UNUR_PAR *parameters, double factor */
				 result = _unur_str_par_set_d(par,key,type_args,args,unur_tdr_set_guidefactor);
				 break;
			 }
			 break;
		 case 'm':
			 if ( !strcmp(key, "max_intervals") ) {
				 /* n = 1; type = i:  UNUR_PAR *parameters, int max_ivs */
				 result = _unur_str_par_set_i(par,key,type_args,args,unur_tdr_set_max_intervals);
				 break;
			 }
			 if ( !strcmp(key, "max_sqhratio") ) {
				 /* n = 1; type = d:  UNUR_PAR *parameters, double max_ratio */
				 result = _unur_str_par_set_d(par,key,type_args,args,unur_tdr_set_max_sqhratio);
				 break;
			 }
			 break;
		 case 'p':
			 if ( !strcmp(key, "pedantic") ) {
				 /* n = 1; type = i:  UNUR_PAR *parameters, int pedantic */
				 result = _unur_str_par_set_i(par,key,type_args,args,unur_tdr_set_pedantic);
				 break;
			 }
			 break;
		 case 'r':
			 if ( !strcmp(key, "reinit_ncpoints") ) {
				 /* n = 1; type = i:  UNUR_PAR *parameters, int ncpoints */
				 result = _unur_str_par_set_i(par,key,type_args,args,unur_tdr_set_reinit_ncpoints);
				 break;
			 }
			 if ( !strcmp(key, "reinit_percentiles") ) {
				 /* n = 2; type = iD:  UNUR_PAR *parameters, int n_percentiles, const double *percentiles */
				 result = _unur_str_par_set_iD(par,key,type_args,args,unur_tdr_set_reinit_percentiles,mlist);
				 break;
			 }
			 break;
		 case 'u':
			 if ( !strcmp(key, "usecenter") ) {
				 /* n = 1; type = i:  UNUR_PAR *parameters, int usecenter */
				 result = _unur_str_par_set_i(par,key,type_args,args,unur_tdr_set_usecenter);
				 break;
			 }
			 if ( !strcmp(key, "usedars") ) {
				 /* n = 1; type = i:  UNUR_PAR *parameters, int usedars */
				 result = _unur_str_par_set_i(par,key,type_args,args,unur_tdr_set_usedars);
				 break;
			 }
			 if ( !strcmp(key, "usemode") ) {
				 /* n = 1; type = i:  UNUR_PAR *parameters, int usemode */
				 result = _unur_str_par_set_i(par,key,type_args,args,unur_tdr_set_usemode);
				 break;
			 }
			 break;
		 case 'v':
			 if ( !strcmp(key, "variant_gw") ) {
				 /* n = 0; type = :  UNUR_PAR *parameters */
				 result = _unur_str_par_set_void(par,key,type_args,args,unur_tdr_set_variant_gw);
				 break;
			 }
			 if ( !strcmp(key, "variant_ia") ) {
				 /* n = 0; type = :  UNUR_PAR *parameters */
				 result = _unur_str_par_set_void(par,key,type_args,args,unur_tdr_set_variant_ia);
				 break;
			 }
			 if ( !strcmp(key, "variant_ps") ) {
				 /* n = 0; type = :  UNUR_PAR *parameters */
				 result = _unur_str_par_set_void(par,key,type_args,args,unur_tdr_set_variant_ps);
				 break;
			 }
			 if ( !strcmp(key, "verify") ) {
				 /* n = 1; type = i:  UNUR_PAR *parameters, int verify */
				 result = _unur_str_par_set_i(par,key,type_args,args,unur_tdr_set_verify);
				 break;
			 }
		 }
		 break;
	 case UNUR_METH_UTDR:
		 switch (*key) {
		 case 'c':
			 if ( !strcmp(key, "cpfactor") ) {
				 /* n = 1; type = d:  UNUR_PAR *parameters, double cp_factor */
				 result = _unur_str_par_set_d(par,key,type_args,args,unur_utdr_set_cpfactor);
				 break;
			 }
			 break;
		 case 'd':
			 if ( !strcmp(key, "deltafactor") ) {
				 /* n = 1; type = d:  UNUR_PAR *parameters, double delta */
				 result = _unur_str_par_set_d(par,key,type_args,args,unur_utdr_set_deltafactor);
				 break;
			 }
			 break;
		 case 'p':
			 if ( !strcmp(key, "pdfatmode") ) {
				 /* n = 1; type = d:  UNUR_PAR *parameters, double fmode */
				 result = _unur_str_par_set_d(par,key,type_args,args,unur_utdr_set_pdfatmode);
				 break;
			 }
			 break;
		 case 'v':
			 if ( !strcmp(key, "verify") ) {
				 /* n = 1; type = i:  UNUR_PAR *parameters, int verify */
				 result = _unur_str_par_set_i(par,key,type_args,args,unur_utdr_set_verify);
				 break;
			 }
		 }
		 break;
	 case UNUR_METH_VEMPK:
		 switch (*key) {
		 case 's':
			 if ( !strcmp(key, "smoothing") ) {
				 /* n = 1; type = d:  UNUR_PAR *parameters, double smoothing */
				 result = _unur_str_par_set_d(par,key,type_args,args,unur_vempk_set_smoothing);
				 break;
			 }
			 break;
		 case 'v':
			 if ( !strcmp(key, "varcor") ) {
				 /* n = 1; type = i:  UNUR_PAR *parameters, int varcor */
				 result = _unur_str_par_set_i(par,key,type_args,args,unur_vempk_set_varcor);
				 break;
			 }
		 }
		 break;
	 case UNUR_METH_VNROU:
		 switch (*key) {
		 case 'r':
			 if ( !strcmp(key, "r") ) {
				 /* n = 1; type = d:  UNUR_PAR *parameters, double r */
				 result = _unur_str_par_set_d(par,key,type_args,args,unur_vnrou_set_r);
				 break;
			 }
			 break;
		 case 'v':
			 if ( !strcmp(key, "v") ) {
				 /* n = 1; type = d:  UNUR_PAR *parameters, double vmax */
				 result = _unur_str_par_set_d(par,key,type_args,args,unur_vnrou_set_v);
				 break;
			 }
			 if ( !strcmp(key, "verify") ) {
				 /* n = 1; type = i:  UNUR_PAR *parameters, int verify */
				 result = _unur_str_par_set_i(par,key,type_args,args,unur_vnrou_set_verify);
				 break;
			 }
		 }
		 break;
	 }

	 /* Ignored set commands: */
 	 /* int unur_mixt_set_useinversion( UNUR_PAR *parameters, int useinv );
		 n = 1; type = i	 */


	 /* Unsupported set commands: */
 	 /* int unur_cext_set_init( UNUR_PAR *parameters, int (*init)(UNUR_GEN *gen) );
		 n = 1; type = I	 */
	 /* int unur_cext_set_sample( UNUR_PAR *parameters, double (*sample)(UNUR_GEN *gen) );
		 n = 1; type = D	 */
	 /* int unur_dext_set_init( UNUR_PAR *parameters, int (*init)(UNUR_GEN *gen) );
		 n = 1; type = I	 */
	 /* int unur_dext_set_sample( UNUR_PAR *parameters, int (*sample)(UNUR_GEN *gen) );
		 n = 1; type = I	 */
	 /* int unur_empk_set_kernelgen( UNUR_PAR *parameters, const UNUR_GEN *kernelgen, double alpha, double kernelvar );
		 n = 3; type = ?dd	 */
	 /* int unur_gibbs_set_startingpoint( UNUR_PAR *parameters, const double *x0 );
		 n = 1; type = D	 */
	 /* int unur_hitro_set_u( UNUR_PAR *parameters, const double *umin, const double *umax );
		 n = 2; type = DD	 */
	 /* int unur_hitro_set_startingpoint( UNUR_PAR *parameters, const double *x0 );
		 n = 1; type = D	 */
	 /* int unur_mcorr_set_eigenvalues( UNUR_PAR *par, const double *eigenvalues );
		 n = 1; type = D	 */
	 /* int unur_vnrou_set_u( UNUR_PAR *parameters, double *umin, double *umax );
		 n = 2; type = DD	 */


   if (result == UNUR_ERR_STR_UNKNOWN) {
     /* no valid parameter found */
     /* try extra keywords */
     if ( !strcmp(key, "debug") ) {
       /* n = 1; type = u:  UNUR_PAR *parameters, unsigned DEBUG */
       result = _unur_str_par_set_u(par,key,type_args,args,unur_set_debug);
     }
   }
  
  /* check result */
  if (result == UNUR_ERR_STR_UNKNOWN) {
    /* unknown parameter */
    _unur_error_unknown(key,"parameter for given method");
    return UNUR_ERR_STR_UNKNOWN;
  }
  else if (result != UNUR_SUCCESS) {
    /* invalid data for parameter */
    _unur_error_invalid(key,"set call");
    return UNUR_ERR_STR_SYNTAX;
  }
  
  return UNUR_SUCCESS;

} /* end of _unur_str_par_set() */

/*---------------------------------------------------------------------------*/

