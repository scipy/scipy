
/* Cephes module version 1.5
 *  This module defines the functions in the cephes and amos libraries as
 *   Numerical python ufunc objects so that they can operate on arbitrary 
 *   NumPy arrays with broadcasting and typecasting rules implemented.
 *  
 *  Copyright 1999  Travis E. Oliphant
 * Revisions 2002 (added functions from cdflib)
 */

#include "Python.h"
#include "Numeric/arrayobject.h"
#include "Numeric/ufuncobject.h" 
#include "ufunc_extras.h"
#include "abstract.h"
#include "cephes.h"
#include "amos_wrappers.h"
#include "toms_wrappers.h"
#include "cdf_wrappers.h"
#include "specfun_wrappers.h"
#include "c_misc/misc.h"
#ifdef macintosh
#include "mymath.h"
#else
#include <math.h>
#endif

/* Defined in mtherr in the cephes library */
extern int scipy_special_print_error_messages;
 
#include "cephes_doc.h"


static PyUFuncGenericFunction cephes1_functions[] = { NULL, NULL, };
static PyUFuncGenericFunction cephes1rc_functions[] = { NULL, NULL, NULL, NULL};
static PyUFuncGenericFunction cephes1_2_functions[] = { NULL, NULL, NULL, NULL,};
static PyUFuncGenericFunction cephes1_2c_functions[] = { NULL, NULL,};
static PyUFuncGenericFunction cephes1c_4_functions[] = { NULL, NULL, NULL, NULL };
static PyUFuncGenericFunction cephes1cp_4_functions[] = { NULL, NULL, NULL, NULL};
static PyUFuncGenericFunction cephes1cpb_4_functions[] = { NULL, NULL,};
static PyUFuncGenericFunction cephes2_functions[] = { NULL, NULL, };
static PyUFuncGenericFunction cephes2_2_functions[] = { NULL, NULL, };
static PyUFuncGenericFunction cephes2_4_functions[] = { NULL, NULL, };
static PyUFuncGenericFunction cephes2a_functions[] = { NULL, NULL, };
static PyUFuncGenericFunction cephes2c_functions[] = { NULL, NULL, NULL, NULL };
static PyUFuncGenericFunction cephes2cp_functions[] = { NULL, NULL, NULL, NULL, };
static PyUFuncGenericFunction cephes2cpp_functions[] = { NULL, NULL, };
static PyUFuncGenericFunction cephes3_functions[] = { NULL, NULL, NULL, NULL};
static PyUFuncGenericFunction cephes3a_functions[] = { NULL, NULL, };
static PyUFuncGenericFunction cephes3_2_functions[] = { NULL, NULL,};
static PyUFuncGenericFunction cephes4_functions[] = { NULL, NULL, NULL, NULL,};
static PyUFuncGenericFunction cephes4a_2_functions[] = { NULL, NULL, };
static PyUFuncGenericFunction cephes4_2_functions[] = { NULL, NULL, };
static PyUFuncGenericFunction cephes5_2_functions[] = { NULL, NULL, };

static PyUFuncGenericFunction cephes1c_functions[] = { NULL, };

static void * airy_data[] = { (void *)airy, (void *)airy, (void *)cairy_wrap, (void *)cairy_wrap,};
static void * airye_data[] = { (void *)cairy_wrap_e, (void *)cairy_wrap_e, };
static void * itairy_data[] = { (void *)itairy_wrap, (void *)itairy_wrap, };

static void * kelvin_data[] = { (void *)kelvin_wrap, (void *)kelvin_wrap,};
static void * ber_data[] = { (void *)ber_wrap, (void *)ber_wrap,};
static void * bei_data[] = { (void *)bei_wrap, (void *)bei_wrap,};
static void * ker_data[] = { (void *)ker_wrap, (void *)ker_wrap,};
static void * kei_data[] = { (void *)kei_wrap, (void *)kei_wrap,};
static void * berp_data[] = { (void *)berp_wrap, (void *)berp_wrap,};
static void * beip_data[] = { (void *)beip_wrap, (void *)beip_wrap,};
static void * kerp_data[] = { (void *)kerp_wrap, (void *)kerp_wrap,};
static void * keip_data[] = { (void *)keip_wrap, (void *)keip_wrap,};

static void * ellpj_data[] = { (void *)ellpj, (void *)ellpj,};

static void * exp1_data[] = { (void *)exp1_wrap, (void *)exp1_wrap, (void *)cexp1_wrap, (void *)cexp1_wrap,};
static void * expi_data[] = { (void *)expi_wrap, (void *)expi_wrap,};
static void * expn_data[] = { (void *)expn, (void *)expn, };
static void * jn_data[] = { (void *)jn, (void *)jn, };
static void * kn_data[] = { (void *)kn, (void *)kn, };

static void * pdtrc_data[] = { (void *)pdtrc, (void *)pdtrc, };
static void * pdtr_data[] = { (void *)pdtr, (void *)pdtr, };
static void * pdtri_data[] = { (void *)pdtri, (void *)pdtri, };

static void * fresnl_data[] = { (void *)fresnl, (void *)fresnl, 
                                (void *)cfresnl_wrap, (void *)cfresnl_wrap };
static void * shichi_data[] = { (void *)shichi, (void *)shichi, };
static void * sici_data[] = { (void *)sici, (void *)sici, };


static void * itj0y0_data[] = { (void *)it1j0y0_wrap, (void *)it1j0y0_wrap, };
static void * it2j0y0_data[] = { (void *)it2j0y0_wrap, (void *)it2j0y0_wrap, };
static void * iti0k0_data[] = { (void *)it1i0k0_wrap, (void *)it1i0k0_wrap, };
static void * it2i0k0_data[] = { (void *)it2i0k0_wrap, (void *)it2i0k0_wrap, };

/*
static void * stdtr_data[] = { (void *)stdtr, (void *)stdtr, };
static void * stdtri_data[] = { (void *)stdtri, (void *)stdtri, };
*/

static void * yn_data[] = { (void *)yn, (void *)yn, };
static void * smirnov_data[] = { (void *)smirnov, (void *)smirnov, };
static void * smirnovi_data[] = { (void *)smirnovi, (void *)smirnovi, };

static void * bdtrc_data[] = { (void *)bdtrc, (void *)bdtrc, };
static void * bdtr_data[] = { (void *)bdtr, (void *)bdtr, };
static void * bdtri_data[] = { (void *)bdtri, (void *)bdtri, };
static void * btdtr_data[] = { (void *)btdtr, (void *)btdtr, };
static void * btdtri_data[] = { (void *)incbi, (void *)incbi, };

static void * fdtrc_data[] = { (void *)fdtrc, (void *)fdtrc, };
static void * fdtr_data[] = { (void *)fdtr, (void *)fdtr, };
static void * fdtri_data[] = { (void *)fdtri, (void *)fdtri, };

static void * gdtrc_data[] = { (void *)gdtrc, (void *)gdtrc, };
static void * gdtr_data[] = { (void *)gdtr, (void *)gdtr, };
/*
static void * gdtri_data[] = { (void *)gdtri, (void *)gdtri, };
*/
static void * hyp2f1_data[] = { (void *)hyp2f1, (void *)hyp2f1, (void *)chyp2f1_wrap, (void *)chyp2f1_wrap};
static void * hyperg_data[] = { (void *)hyperg, (void *)hyperg, (void *)chyp1f1_wrap, (void *)chyp1f1_wrap};
static void * hypU_data[] = { (void *)hypU_wrap, (void *)hypU_wrap, };
static void * hyp2f0_data[] = { (void *)hyp2f0, (void *)hyp2f0, };
static void * threef0_data[] = { (void *)threef0, (void *)threef0, };
static void * onef2_data[] = { (void *)onef2, (void *)onef2, };

static void * incbet_data[] = { (void *)incbet, (void *)incbet, };
static void * incbi_data[] = { (void *)incbi, (void *)incbi, };

static void * nbdtrc_data[] = { (void *)nbdtrc, (void *)nbdtrc, };
static void * nbdtr_data[] = { (void *)nbdtr, (void *)nbdtr, };
static void * nbdtri_data[] = { (void *)nbdtri, (void *)nbdtri, };

static void * beta_data[] = { (void *)beta, (void *)beta, };
static void * lbeta_data[] = { (void *)lbeta, (void *)lbeta, };
static void * cbrt_data[] = { (void *)cbrt, (void *)cbrt, };
static void * chdtrc_data[] = { (void *)chdtrc, (void *)chdtrc, };
static void * chdtr_data[] = { (void *)chdtr, (void *)chdtr, };
static void * chdtri_data[] = { (void *)chdtri, (void *)chdtri, };
static void * dawsn_data[] = {  (void *)dawsn, (void *)dawsn, };
static void * ellie_data[] = { (void *)ellie, (void *)ellie, };
static void * ellik_data[] = { (void *)ellik, (void *)ellik, };
static void * ellpe_data[] = { (void *)ellpe, (void *)ellpe, };
static void * ellpk_data[] = { (void *)ellpk, (void *)ellpk, };
static void * exp10_data[] = { (void *)exp10, (void *)exp10, };
static void * exp2_data[] = { (void *)exp2, (void *)exp2, };
static void * Gamma_data[] = { (void *)Gamma, (void *)Gamma, (void *)cgamma_wrap, (void *)cgamma_wrap};
static void * lgam_data[] = { (void *)lgam, (void *)lgam, (void *)clngamma_wrap, (void *)clngamma_wrap};
static void * i0_data[] = { (void *)i0, (void *)i0, };
static void * i0e_data[] = { (void *)i0e, (void *)i0e, };
static void * i1_data[] = { (void *)i1, (void *)i1, };
static void * i1e_data[] = { (void *)i1e, (void *)i1e, };
static void * igamc_data[] = { (void *)igamc, (void *)igamc, };
static void * igam_data[] = { (void *)igam, (void *)igam, };
static void * igami_data[] = { (void *)igami, (void *)igami, };

static void * iv_data[] = { (void *)iv, (void *)iv, (void *)cbesi_wrap, (void *)cbesi_wrap,};
static void * ive_data[] = { (void *)cbesi_wrap_e, (void *)cbesi_wrap_e, (void *)cbesi_wrap_e, (void *)cbesi_wrap_e, };
static void * j0_data[] = { (void *)j0,  (void *)j0,  };
static void * y0_data[] = { (void *)y0, (void *)y0, };
static void * j1_data[] = { (void *)j1,  (void *)j1,  };
static void * y1_data[] = { (void *)y1, (void *)y1, };
static void * jv_data[] = { (void *)jv, (void *)jv, (void *)cbesj_wrap, (void *)cbesj_wrap,};
static void * jve_data[] = { (void *)cbesj_wrap_e, (void *)cbesj_wrap_e, (void *)cbesj_wrap_e, (void *)cbesj_wrap_e, };
static void * yv_data[] = { (void *)yv, (void *)yv, (void *)cbesy_wrap, (void *)cbesy_wrap,};
static void * yve_data[] = { (void *)cbesy_wrap_e, (void *)cbesy_wrap_e, (void *)cbesy_wrap_e, (void *)cbesy_wrap_e, };

static void * k0_data[] = { (void *)k0, (void *)k0, };
static void * k0e_data[] = { (void *)k0e, (void *)k0e, };
static void * k1_data[] = { (void *)k1, (void *)k1, };
static void * k1e_data[] = { (void *)k1e, (void *)k1e, };
static void * kv_data[] = { (void *)cbesk_wrap, (void *)cbesk_wrap, (void *)cbesk_wrap, (void *)cbesk_wrap,};
static void * kve_data[] = { (void *)cbesk_wrap_e, (void *)cbesk_wrap_e, (void *)cbesk_wrap_e, (void *)cbesk_wrap_e,};
static void * hankel1_data[] = { (void *)cbesh_wrap1, (void *)cbesh_wrap1,};
static void * hankel1e_data[] = { (void *)cbesh_wrap1_e, (void *)cbesh_wrap1_e,};
static void * hankel2_data[] = { (void *)cbesh_wrap2, (void *)cbesh_wrap2,};
static void * hankel2e_data[] = { (void *)cbesh_wrap2_e, (void *)cbesh_wrap2_e,};

static void * ndtr_data[] = { (void *)ndtr, (void *)ndtr, };
static void * erfc_data[] = { (void *)erfc, (void *)erfc, };
static void * erf_data[] = { (void *)erf, (void *)erf, (void *)cerf_wrap, (void *)cerf_wrap};
static void * ndtri_data[] = { (void *)ndtri, (void *)ndtri, };

static void * psi_data[] = { (void *)psi, (void *)psi, (void *)cpsi_wrap, (void *)cpsi_wrap};
static void * rgamma_data[] = { (void *)rgamma, (void *)rgamma, (void *)crgamma_wrap, (void *)crgamma_wrap};
static void * round_data[] = { (void *)round, (void *)round, };
static void * sindg_data[] = { (void *)sindg, (void *)sindg, };
static void * cosdg_data[] = { (void *)cosdg, (void *)cosdg, };
static void * radian_data[] = { (void *)radian, (void *)radian, };
static void * tandg_data[] = { (void *)tandg, (void *)tandg, };
static void * cotdg_data[] = { (void *)cotdg, (void *)cotdg, };
static void * log1p_data[] = { (void *)log1p, (void *)log1p, };
static void * expm1_data[] = { (void *)expm1, (void *)expm1, };
static void * cosm1_data[] = { (void *)cosm1, (void *)cosm1, };

static void * spence_data[] = { (void *)spence, (void *)spence, };
/* static void * struve_data[] = { (void *)struve, (void *)struve, };*/
static void * struve_data[] = { (void *)struve_wrap, (void *)struve_wrap, };
static void * modstruve_data[] = { (void *)modstruve_wrap, (void *)modstruve_wrap, };
static void * itmodstruve0_data[] = { (void *)itmodstruve0_wrap, (void *)itmodstruve0_wrap, };
static void * itstruve0_data[] = { (void *)itstruve0_wrap, (void *)itstruve0_wrap, };
static void * it2struve0_data[] = { (void *)it2struve0_wrap, (void *)it2struve0_wrap, };


static void * zeta_data[] = { (void *)zeta, (void *)zeta, };
static void * zetac_data[] = { (void *)zetac, (void *)zetac, };

static void * kolmogorov_data[] = { (void *)kolmogorov, (void *)kolmogorov, };
static void * kolmogi_data[] = { (void *)kolmogi, (void *)kolmogi, };

static void * wofz_data[] = { (void *)cwofz_wrap, (void *)cwofz_wrap, };

static void * besselpoly_data[] = {(void *)besselpoly, (void *)besselpoly,};

static void * cdfbet3_data[] = {(void *)cdfbet3_wrap, (void *)cdfbet3_wrap};
static void * cdfbet4_data[] = {(void *)cdfbet4_wrap, (void *)cdfbet4_wrap};
static void * cdfbin2_data[] = {(void *)cdfbin2_wrap, (void *)cdfbin2_wrap};
static void * cdfbin3_data[] = {(void *)cdfbin3_wrap, (void *)cdfbin3_wrap};
static void * cdfchi3_data[] = {(void *)cdfchi3_wrap, (void *)cdfchi3_wrap};
static void * cdfchn1_data[] = {(void *)cdfchn1_wrap, (void *)cdfchn1_wrap};
static void * cdfchn2_data[] = {(void *)cdfchn2_wrap, (void *)cdfchn2_wrap};
static void * cdfchn3_data[] = {(void *)cdfchn3_wrap, (void *)cdfchn3_wrap};
static void * cdfchn4_data[] = {(void *)cdfchn4_wrap, (void *)cdfchn4_wrap};
/*
static void * cdff1_data[] = {(void *)cdff1_wrap, (void *)cdff1_wrap};
static void * cdff2_data[] = {(void *)cdff2_wrap, (void *)cdff2_wrap};
static void * cdff3_data[] = {(void *)cdff3_wrap, (void *)cdff3_wrap};
*/
static void * cdff4_data[] = {(void *)cdff4_wrap, (void *)cdff4_wrap};

static void * cdffnc1_data[] = {(void *)cdffnc1_wrap, (void *)cdffnc1_wrap};
static void * cdffnc2_data[] = {(void *)cdffnc2_wrap, (void *)cdffnc2_wrap};
static void * cdffnc3_data[] = {(void *)cdffnc3_wrap, (void *)cdffnc3_wrap};
static void * cdffnc4_data[] = {(void *)cdffnc4_wrap, (void *)cdffnc4_wrap};
static void * cdffnc5_data[] = {(void *)cdffnc5_wrap, (void *)cdffnc5_wrap};
/*
static void * cdfgam1_data[] = {(void *)cdfgam1_wrap, (void *)cdfgam1_wrap};
*/
static void * cdfgam2_data[] = {(void *)cdfgam2_wrap, (void *)cdfgam2_wrap};
static void * cdfgam3_data[] = {(void *)cdfgam3_wrap, (void *)cdfgam3_wrap};
static void * cdfgam4_data[] = {(void *)cdfgam4_wrap, (void *)cdfgam4_wrap};

static void * cdfnbn2_data[] = {(void *)cdfnbn2_wrap, (void *)cdfnbn2_wrap};
static void * cdfnbn3_data[] = {(void *)cdfnbn3_wrap, (void *)cdfnbn3_wrap};

static void * cdfnor3_data[] = {(void *)cdfnor3_wrap, (void *)cdfnor3_wrap};
static void * cdfnor4_data[] = {(void *)cdfnor4_wrap, (void *)cdfnor4_wrap};

static void * cdfpoi2_data[] = {(void *)cdfpoi2_wrap, (void *)cdfpoi2_wrap};

static void * cdft1_data[] = {(void *)cdft1_wrap, (void *)cdft1_wrap};
static void * cdft2_data[] = {(void *)cdft2_wrap, (void *)cdft2_wrap};
static void * cdft3_data[] = {(void *)cdft3_wrap, (void *)cdft3_wrap};

static void * cdftnc1_data[] = {(void *)cdftnc1_wrap, (void *)cdftnc1_wrap};
static void * cdftnc2_data[] = {(void *)cdftnc2_wrap, (void *)cdftnc2_wrap};
static void * cdftnc3_data[] = {(void *)cdftnc3_wrap, (void *)cdftnc3_wrap};
static void * cdftnc4_data[] = {(void *)cdftnc4_wrap, (void *)cdftnc4_wrap};

static void * tklambda_data[] = {(void *)tukeylambdacdf, (void *)tukeylambdacdf};

static void * mathieu_a_data[] = {(void *)cem_cva_wrap, (void *)cem_cva_wrap};
static void * mathieu_b_data[] = {(void *)sem_cva_wrap, (void *)sem_cva_wrap};
static void * mathieu_cem_data[] = {(void *)cem_wrap, (void *)cem_wrap};
static void * mathieu_sem_data[] = {(void *)sem_wrap, (void *)sem_wrap};
static void * mathieu_mcem1_data[] = {(void *)mcm1_wrap, (void *)mcm1_wrap};
static void * mathieu_mcem2_data[] = {(void *)mcm2_wrap, (void *)mcm2_wrap};
static void * mathieu_msem1_data[] = {(void *)msm1_wrap, (void *)msm1_wrap};
static void * mathieu_msem2_data[] = {(void *)msm2_wrap, (void *)msm2_wrap};

static void * lpmv_data[] = {(void *)pmv_wrap, (void *)pmv_wrap};
static void * pbwa_data[] = {(void *)pbwa_wrap, (void *)pbwa_wrap};
static void * pbdv_data[] = {(void *)pbdv_wrap, (void *)pbdv_wrap};
static void * pbvv_data[] = {(void *)pbvv_wrap, (void *)pbvv_wrap};
static void * prolate_aswfa_data[] = {(void *)prolate_aswfa_wrap, (void *)prolate_aswfa_wrap};
static void * prolate_radial1_data[] = {(void *)prolate_radial1_wrap, (void *)prolate_radial1_wrap};
static void * prolate_radial2_data[] = {(void *)prolate_radial2_wrap, (void *)prolate_radial2_wrap};
static void * oblate_aswfa_data[] = {(void *)oblate_aswfa_wrap, (void *)oblate_aswfa_wrap};
static void * oblate_radial1_data[] = {(void *)oblate_radial1_wrap, (void *)oblate_radial1_wrap};
static void * oblate_radial2_data[] = {(void *)oblate_radial2_wrap, (void *)oblate_radial2_wrap};
static void * prolate_aswfa_nocv_data[] = {(void *)prolate_aswfa_nocv_wrap, (void *)prolate_aswfa_nocv_wrap};
static void * prolate_radial1_nocv_data[] = {(void *)prolate_radial1_nocv_wrap, (void *)prolate_radial1_nocv_wrap};
static void * prolate_radial2_nocv_data[] = {(void *)prolate_radial2_nocv_wrap, (void *)prolate_radial2_nocv_wrap};
static void * oblate_aswfa_nocv_data[] = {(void *)oblate_aswfa_nocv_wrap, (void *)oblate_aswfa_nocv_wrap};
static void * oblate_radial1_nocv_data[] = {(void *)oblate_radial1_nocv_wrap, (void *)oblate_radial1_nocv_wrap};
static void * oblate_radial2_nocv_data[] = {(void *)oblate_radial2_nocv_wrap, (void *)oblate_radial2_nocv_wrap};
static void * prolate_segv_data[] = {(void *)prolate_segv_wrap, (void *)prolate_segv_wrap};
static void * oblate_segv_data[] = {(void *)oblate_segv_wrap, (void *)oblate_segv_wrap};

static void * modfresnelp_data[] = {(void *)modified_fresnel_plus_wrap, (void *)modified_fresnel_plus_wrap};
static void * modfresnelm_data[] = {(void *)modified_fresnel_minus_wrap, (void *)modified_fresnel_minus_wrap};


static char cephes_7_types[] = { PyArray_FLOAT,  PyArray_FLOAT,  PyArray_FLOAT, PyArray_FLOAT, PyArray_FLOAT, PyArray_FLOAT, PyArray_FLOAT, PyArray_DOUBLE,  PyArray_DOUBLE, PyArray_DOUBLE, PyArray_DOUBLE, PyArray_DOUBLE, PyArray_DOUBLE, PyArray_DOUBLE,};
static char cephes_6_types[] = { PyArray_FLOAT,  PyArray_FLOAT,  PyArray_FLOAT, PyArray_FLOAT, PyArray_FLOAT, PyArray_FLOAT, PyArray_DOUBLE,  PyArray_DOUBLE, PyArray_DOUBLE, PyArray_DOUBLE, PyArray_DOUBLE, PyArray_DOUBLE,};
static char cephes_5_types[] = { PyArray_FLOAT,  PyArray_FLOAT,  PyArray_FLOAT, PyArray_FLOAT, PyArray_FLOAT, PyArray_DOUBLE,  PyArray_DOUBLE, PyArray_DOUBLE, PyArray_DOUBLE, PyArray_DOUBLE,};

static char cephes_5b2_types[] = { PyArray_FLOAT,  PyArray_CFLOAT,  PyArray_CFLOAT, PyArray_CFLOAT, PyArray_CFLOAT, PyArray_DOUBLE,  PyArray_CDOUBLE, PyArray_CDOUBLE, PyArray_CDOUBLE, PyArray_CDOUBLE,};

static char cephes_5c_types[] = { PyArray_FLOAT,  PyArray_FLOAT,  PyArray_FLOAT, PyArray_FLOAT, PyArray_FLOAT, PyArray_DOUBLE,  PyArray_DOUBLE, PyArray_DOUBLE, PyArray_DOUBLE, PyArray_DOUBLE, PyArray_CFLOAT, PyArray_CFLOAT, PyArray_CFLOAT, PyArray_CFLOAT, PyArray_CFLOAT, PyArray_CDOUBLE, PyArray_CDOUBLE, PyArray_CDOUBLE, PyArray_CDOUBLE, PyArray_CDOUBLE, };

static char cephes_5c2_types[] = { PyArray_FLOAT,  PyArray_FLOAT,  PyArray_FLOAT, PyArray_FLOAT, PyArray_FLOAT, PyArray_DOUBLE,  PyArray_DOUBLE, PyArray_DOUBLE, PyArray_DOUBLE, PyArray_DOUBLE, PyArray_FLOAT, PyArray_FLOAT, PyArray_FLOAT, PyArray_CFLOAT, PyArray_CFLOAT, PyArray_DOUBLE, PyArray_DOUBLE, PyArray_DOUBLE, PyArray_CDOUBLE, PyArray_CDOUBLE, };

static char cephes_4_types[] = { PyArray_FLOAT,  PyArray_FLOAT,  PyArray_FLOAT, PyArray_FLOAT, PyArray_DOUBLE,  PyArray_DOUBLE, PyArray_DOUBLE, PyArray_DOUBLE,};

static char cephes_4c_types[] = { PyArray_FLOAT,  PyArray_FLOAT,  PyArray_FLOAT, PyArray_FLOAT, PyArray_DOUBLE,  PyArray_DOUBLE, PyArray_DOUBLE, PyArray_DOUBLE, PyArray_FLOAT, PyArray_FLOAT, PyArray_CFLOAT, PyArray_CFLOAT, PyArray_DOUBLE,  PyArray_DOUBLE, PyArray_CDOUBLE, PyArray_CDOUBLE};

static char cephes_3_types[] = { PyArray_FLOAT,  PyArray_FLOAT,  PyArray_FLOAT,   PyArray_DOUBLE,  PyArray_DOUBLE, PyArray_DOUBLE, };
static char cephes_3_cmplx_types[] = { PyArray_FLOAT,  PyArray_FLOAT,  PyArray_FLOAT,   PyArray_DOUBLE,  PyArray_DOUBLE, PyArray_DOUBLE, PyArray_CFLOAT,  PyArray_CFLOAT,  PyArray_CFLOAT,   PyArray_CDOUBLE,  PyArray_CDOUBLE, PyArray_CDOUBLE, };
static char cephes_3c_types[] = { PyArray_FLOAT, PyArray_FLOAT, PyArray_FLOAT, PyArray_DOUBLE, PyArray_DOUBLE, PyArray_DOUBLE, PyArray_FLOAT, PyArray_CFLOAT,  PyArray_CFLOAT, PyArray_DOUBLE, PyArray_CDOUBLE, PyArray_CDOUBLE, };
static char cephes_3cp_types[] = { PyArray_FLOAT, PyArray_CFLOAT,  PyArray_CFLOAT, PyArray_DOUBLE, PyArray_CDOUBLE, PyArray_CDOUBLE, };
static char cephes_2_types[] = { PyArray_FLOAT,  PyArray_FLOAT,  PyArray_DOUBLE,  PyArray_DOUBLE,  };
static char cephes_1rc_types[] = { PyArray_FLOAT,  PyArray_FLOAT,  PyArray_DOUBLE,  PyArray_DOUBLE,  PyArray_CFLOAT, PyArray_CFLOAT, PyArray_CDOUBLE, PyArray_CDOUBLE };
static char cephes_1c_types[] = { PyArray_CFLOAT, PyArray_CFLOAT, PyArray_CDOUBLE, PyArray_CDOUBLE, };


/* Some functions needed from ufunc object, so that Py_complex's aren't being returned 
between code possibly compiled with different compilers.
*/

typedef Py_complex ComplexUnaryFunc(Py_complex x);

static void cephes_F_F_As_D_D(char **args, int *dimensions, int *steps, void *func) {
    int i; Py_complex x;
    char *ip1=args[0], *op=args[1];
    for(i=0; i<*dimensions; i++, ip1+=steps[0], op+=steps[1]) {
	x.real = ((float *)ip1)[0]; x.imag = ((float *)ip1)[1];
	x = ((ComplexUnaryFunc *)func)(x);
	((float *)op)[0] = (float)x.real;
	((float *)op)[1] = (float)x.imag;
    }
}

static void cephes_D_D(char **args, int *dimensions, int *steps, void *func) {
    int i; Py_complex x;
    char *ip1=args[0], *op=args[1];
    for(i=0; i<*dimensions; i++, ip1+=steps[0], op+=steps[1]) {
	x.real = ((double *)ip1)[0]; x.imag = ((double *)ip1)[1];
	x = ((ComplexUnaryFunc *)func)(x);
	((double *)op)[0] = x.real;
	((double *)op)[1] = x.imag;
    }
}

static void Cephes_InitOperators(PyObject *dictionary) {
	PyObject *f;

        cephes1_functions[0] = PyUFunc_f_f_As_d_d;
        cephes1_functions[1] = PyUFunc_d_d;
        cephes1c_functions[0] = cephes_F_F_As_D_D;
	cephes1c_functions[1] = cephes_D_D;
        cephes1rc_functions[0] = PyUFunc_f_f_As_d_d;
        cephes1rc_functions[1] = PyUFunc_d_d;
        cephes1rc_functions[2] = cephes_F_F_As_D_D;
	cephes1rc_functions[3] = cephes_D_D;
        cephes1_2_functions[0] = PyUFunc_f_ff_As_d_dd;
        cephes1_2_functions[1] = PyUFunc_d_dd;
        cephes1_2_functions[2] = PyUFunc_F_FF_As_D_DD;
        cephes1_2_functions[3] = PyUFunc_D_DD;
        cephes1_2c_functions[0] = PyUFunc_f_FF_As_d_DD;
        cephes1_2c_functions[1] = PyUFunc_d_DD;
        cephes1c_4_functions[0] = PyUFunc_f_ffff_As_d_dddd;
        cephes1c_4_functions[1] = PyUFunc_d_dddd;
        cephes1c_4_functions[2] = PyUFunc_F_FFFF_As_D_DDDD;
        cephes1c_4_functions[3] = PyUFunc_D_DDDD;
        cephes1cp_4_functions[0] = PyUFunc_f_ffff_As_D_DDDD;
        cephes1cp_4_functions[1] = PyUFunc_d_dddd_As_D_DDDD;
        cephes1cp_4_functions[2] = PyUFunc_F_FFFF_As_D_DDDD;
        cephes1cp_4_functions[3] = PyUFunc_D_DDDD;
        cephes1cpb_4_functions[0] = PyUFunc_f_FFFF_As_d_DDDD;
        cephes1cpb_4_functions[1] = PyUFunc_d_DDDD;
        cephes2_functions[0] = PyUFunc_ff_f_As_dd_d;
        cephes2_functions[1] = PyUFunc_dd_d;
        cephes2_2_functions[0] = PyUFunc_ff_ff_As_dd_dd;
        cephes2_2_functions[1] = PyUFunc_dd_dd;
        cephes2a_functions[0] = PyUFunc_ff_f_As_id_d;
        cephes2a_functions[1] = PyUFunc_dd_d_As_id_d;
        cephes2c_functions[0] = PyUFunc_ff_f_As_dd_d;
        cephes2c_functions[1] = PyUFunc_dd_d;
        cephes2c_functions[2] = PyUFunc_fF_F_As_dD_D;
        cephes2c_functions[3] = PyUFunc_dD_D;
        cephes2cp_functions[0] = PyUFunc_ff_f_As_dD_D;
        cephes2cp_functions[1] = PyUFunc_dd_d_As_dD_D;
        cephes2cp_functions[2] = PyUFunc_fF_F_As_dD_D;
        cephes2cp_functions[3] = PyUFunc_dD_D;
        cephes2cpp_functions[0] = PyUFunc_fF_F_As_dD_D;
        cephes2cpp_functions[1] = PyUFunc_dD_D;
        cephes2_4_functions[0] = PyUFunc_ff_ffff_As_dd_dddd;
        cephes2_4_functions[1] = PyUFunc_dd_dddd;
        cephes3_functions[0] = PyUFunc_fff_f_As_ddd_d;
        cephes3_functions[1] = PyUFunc_ddd_d;
        cephes3_functions[2] = PyUFunc_ffF_F_As_ddD_D;
        cephes3_functions[3] = PyUFunc_ddD_D;
        cephes3a_functions[0] = PyUFunc_fff_f_As_iid_d;
        cephes3a_functions[1] = PyUFunc_ddd_d_As_iid_d;
        cephes3_2_functions[0] = PyUFunc_fff_ff_As_ddd_dd;
        cephes3_2_functions[1] = PyUFunc_ddd_dd;
        cephes4_functions[0] = PyUFunc_ffff_f_As_dddd_d;
        cephes4_functions[1] = PyUFunc_dddd_d;
        cephes4_functions[2] = PyUFunc_fffF_F_As_dddD_D;
        cephes4_functions[3] = PyUFunc_dddD_D;
        cephes4_2_functions[0] = PyUFunc_ffff_ff_As_dddd_dd;
        cephes4_2_functions[1] = PyUFunc_dddd_dd;
        cephes4a_2_functions[0] = PyUFunc_ffff_ff_As_dddi_dd;
        cephes4a_2_functions[1] = PyUFunc_dddd_dd_As_dddi_dd;
        cephes5_2_functions[0] = PyUFunc_fffff_ff_As_ddddd_dd;
        cephes5_2_functions[1] = PyUFunc_ddddd_dd;
	
	/* Create function objects for each function call and insert
	   them in the dictionary */
	f = PyUFunc_FromFuncAndData(cephes3a_functions, bdtrc_data, cephes_4_types, 2, 3, 1, PyUFunc_None, "bdtrc", bdtrc_doc, 0);
	PyDict_SetItemString(dictionary, "bdtrc", f);
	Py_DECREF(f);
	f = PyUFunc_FromFuncAndData(cephes3a_functions, bdtr_data, cephes_4_types, 2, 3, 1, PyUFunc_None, "bdtr", bdtr_doc, 0);
	PyDict_SetItemString(dictionary, "bdtr", f);
	Py_DECREF(f);
	f = PyUFunc_FromFuncAndData(cephes3a_functions, bdtri_data, cephes_4_types, 2, 3, 1, PyUFunc_None, "bdtri", bdtri_doc, 0);
	PyDict_SetItemString(dictionary, "bdtri", f);
	Py_DECREF(f);
	f = PyUFunc_FromFuncAndData(cephes3_functions, btdtr_data, cephes_4_types, 2, 3, 1, PyUFunc_None, "btdtr", btdtr_doc, 0);
	PyDict_SetItemString(dictionary, "btdtr", f);
	Py_DECREF(f);
	f = PyUFunc_FromFuncAndData(cephes3_functions, btdtri_data, cephes_4_types, 2, 3, 1, PyUFunc_None, "btdtri", btdtri_doc, 0);
	PyDict_SetItemString(dictionary, "btdtri", f);
	Py_DECREF(f);

	f = PyUFunc_FromFuncAndData(cephes3_functions, fdtrc_data, cephes_4_types, 2, 3, 1, PyUFunc_None, "fdtrc", fdtrc_doc, 0);
        PyDict_SetItemString(dictionary, "fdtrc", f);
        Py_DECREF(f);
        f = PyUFunc_FromFuncAndData(cephes3_functions, fdtr_data, cephes_4_types, 2, 3, 1, PyUFunc_None, "fdtr", fdtr_doc, 0);
        PyDict_SetItemString(dictionary, "fdtr", f);
	Py_DECREF(f);
	f = PyUFunc_FromFuncAndData(cephes3_functions, fdtri_data, cephes_4_types, 2, 3, 1, PyUFunc_None, "fdtri", fdtri_doc, 0);
	PyDict_SetItemString(dictionary, "fdtri", f);
 	Py_DECREF(f);

	f = PyUFunc_FromFuncAndData(cephes3_functions, gdtrc_data, cephes_4_types, 2, 3, 1, PyUFunc_None, "gdtrc", gdtrc_doc, 0);
	PyDict_SetItemString(dictionary, "gdtrc", f);
	Py_DECREF(f);
	f = PyUFunc_FromFuncAndData(cephes3_functions, gdtr_data, cephes_4_types, 2, 3, 1, PyUFunc_None, "gdtr", gdtr_doc, 0);
	PyDict_SetItemString(dictionary, "gdtr", f);
	Py_DECREF(f);
        /* Use inverse from cdflib (a little faster)
	f = PyUFunc_FromFuncAndData(cephes3_functions, gdtri_data, cephes_4_types, 2, 3, 1, PyUFunc_None, "gdtri", gdtri_doc, 0);
	PyDict_SetItemString(dictionary, "gdtri", f);
	Py_DECREF(f);
        */

	f = PyUFunc_FromFuncAndData(cephes4_functions, hyp2f1_data, cephes_5c2_types, 4, 4, 1, PyUFunc_None, "hyp2f1", hyp2f1_doc, 0);
	PyDict_SetItemString(dictionary, "hyp2f1", f);
	Py_DECREF(f);
	f = PyUFunc_FromFuncAndData(cephes3_functions, hyperg_data, cephes_4c_types, 4, 3, 1, PyUFunc_None, "hyp1f1", hyp1f1_doc, 0);
	PyDict_SetItemString(dictionary, "hyp1f1", f);
	Py_DECREF(f);

	f = PyUFunc_FromFuncAndData(cephes3_functions, hypU_data, cephes_4_types, 2, 3, 1, PyUFunc_None, "hyperu", hyperu_doc, 0);
	PyDict_SetItemString(dictionary, "hyperu", f);
	Py_DECREF(f);

	f = PyUFunc_FromFuncAndData(cephes4a_2_functions, hyp2f0_data, cephes_6_types, 2, 4, 2, PyUFunc_None, "hyp2f0", hyp2f0_doc, 0);
	PyDict_SetItemString(dictionary, "hyp2f0", f);
	Py_DECREF(f);
	f = PyUFunc_FromFuncAndData(cephes4_2_functions, onef2_data, cephes_6_types, 2, 4, 2, PyUFunc_None, "hyp1f2", hyp1f2_doc, 0);
	PyDict_SetItemString(dictionary, "hyp1f2", f);
	Py_DECREF(f);
	f = PyUFunc_FromFuncAndData(cephes4_2_functions, threef0_data, cephes_6_types, 2, 4, 2, PyUFunc_None, "hyp3f0", hyp3f0_doc, 0);
	PyDict_SetItemString(dictionary, "hyp3f0", f);
	Py_DECREF(f);

	f = PyUFunc_FromFuncAndData(cephes3_functions, incbet_data, cephes_4_types, 2, 3, 1, PyUFunc_None, "betainc", betainc_doc, 0);
	PyDict_SetItemString(dictionary, "betainc", f);
	Py_DECREF(f);
	f = PyUFunc_FromFuncAndData(cephes3_functions, incbi_data, cephes_4_types, 2, 3, 1, PyUFunc_None, "betaincinv", betaincinv_doc, 0);
	PyDict_SetItemString(dictionary, "betaincinv", f);
	Py_DECREF(f);

	f = PyUFunc_FromFuncAndData(cephes3a_functions, nbdtrc_data, cephes_4_types, 2, 3, 1, PyUFunc_None, "nbdtrc", nbdtrc_doc, 0);
	PyDict_SetItemString(dictionary, "nbdtrc", f);
	Py_DECREF(f);
	f = PyUFunc_FromFuncAndData(cephes3a_functions, nbdtr_data, cephes_4_types, 2, 3, 1, PyUFunc_None, "nbdtr", nbdtr_doc, 0);
	PyDict_SetItemString(dictionary, "nbdtr", f);
	Py_DECREF(f);
	f = PyUFunc_FromFuncAndData(cephes3a_functions, nbdtri_data, cephes_4_types, 2, 3, 1, PyUFunc_None, "nbdtri", nbdtri_doc, 0);
	PyDict_SetItemString(dictionary, "nbdtri", f);
	Py_DECREF(f);


	f = PyUFunc_FromFuncAndData(cephes2_functions, beta_data, cephes_3_types, 2, 2, 1, PyUFunc_None, "beta", beta_doc, 0);
	PyDict_SetItemString(dictionary, "beta", f);
	Py_DECREF(f);
	f = PyUFunc_FromFuncAndData(cephes2_functions, lbeta_data, cephes_3_types, 2, 2, 1, PyUFunc_None, "betaln", betaln_doc, 0);
	PyDict_SetItemString(dictionary, "betaln", f);
	Py_DECREF(f);
	f = PyUFunc_FromFuncAndData(cephes1_functions, cbrt_data, cephes_2_types, 2, 1, 1, PyUFunc_None, "cbrt", cbrt_doc, 0);
	PyDict_SetItemString(dictionary, "cbrt", f);
	Py_DECREF(f);
	f = PyUFunc_FromFuncAndData(cephes2_functions, chdtrc_data, cephes_3_types, 2, 2, 1, PyUFunc_None, "chdtrc", chdtrc_doc, 0);
	PyDict_SetItemString(dictionary, "chdtrc", f);
	Py_DECREF(f);
	f = PyUFunc_FromFuncAndData(cephes2_functions, chdtr_data, cephes_3_types, 2, 2, 1, PyUFunc_None, "chdtr", chdtr_doc, 0);
	PyDict_SetItemString(dictionary, "chdtr", f);
	Py_DECREF(f);
	f = PyUFunc_FromFuncAndData(cephes2_functions, chdtri_data, cephes_3_types, 2, 2, 1, PyUFunc_None, "chdtri", chdtri_doc, 0);
	PyDict_SetItemString(dictionary, "chdtri", f);
	Py_DECREF(f);
	f = PyUFunc_FromFuncAndData(cephes1_functions, dawsn_data, cephes_2_types, 2, 1, 1, PyUFunc_None, "dawsn", dawsn_doc, 0);
	PyDict_SetItemString(dictionary, "dawsn", f);
	Py_DECREF(f);
	f = PyUFunc_FromFuncAndData(cephes2_functions, ellie_data, cephes_3_types, 2, 2, 1, PyUFunc_None, "ellipeinc", ellipeinc_doc, 0);
	PyDict_SetItemString(dictionary, "ellipeinc", f);
	Py_DECREF(f);
	f = PyUFunc_FromFuncAndData(cephes2_functions, ellik_data, cephes_3_types, 2, 2, 1, PyUFunc_None, "ellipkinc", ellipkinc_doc, 0);
	PyDict_SetItemString(dictionary, "ellipkinc", f);
	Py_DECREF(f);
	f = PyUFunc_FromFuncAndData(cephes1_functions, ellpe_data, cephes_2_types, 2, 1, 1, PyUFunc_None, "ellipe", ellipe_doc, 0);
	PyDict_SetItemString(dictionary, "ellipe", f);
	Py_DECREF(f);
	f = PyUFunc_FromFuncAndData(cephes1_functions, ellpk_data, cephes_2_types, 2, 1, 1, PyUFunc_None, "ellipk", ellipk_doc, 0);
	PyDict_SetItemString(dictionary, "ellipk", f);
	Py_DECREF(f);
	f = PyUFunc_FromFuncAndData(cephes1_functions, exp10_data, cephes_2_types, 2, 1, 1, PyUFunc_None, "exp10", exp10_doc, 0);
	PyDict_SetItemString(dictionary, "exp10", f);
	Py_DECREF(f);

	f = PyUFunc_FromFuncAndData(cephes1_functions, exp2_data, cephes_2_types, 2, 1, 1, PyUFunc_None, "exp2", exp2_doc, 0);
	PyDict_SetItemString(dictionary, "exp2", f);
	Py_DECREF(f);
	f = PyUFunc_FromFuncAndData(cephes1rc_functions, Gamma_data, cephes_1rc_types, 4, 1, 1, PyUFunc_None, "gamma", gamma_doc, 0);
	PyDict_SetItemString(dictionary, "gamma", f);
	Py_DECREF(f);
	f = PyUFunc_FromFuncAndData(cephes1rc_functions, lgam_data, cephes_1rc_types, 4, 1, 1, PyUFunc_None, "gammaln", gammaln_doc, 0);
	PyDict_SetItemString(dictionary, "gammaln", f);
	Py_DECREF(f);
	f = PyUFunc_FromFuncAndData(cephes1_functions, i0_data, cephes_2_types, 2, 1, 1, PyUFunc_None, "i0", i0_doc, 0);
	PyDict_SetItemString(dictionary, "i0", f);
	Py_DECREF(f);
	f = PyUFunc_FromFuncAndData(cephes1_functions, i0e_data, cephes_2_types, 2, 1, 1, PyUFunc_None, "i0e", i0e_doc, 0);
	PyDict_SetItemString(dictionary, "i0e", f);
	Py_DECREF(f);
	f = PyUFunc_FromFuncAndData(cephes1_functions, i1_data, cephes_2_types, 2, 1, 1, PyUFunc_None, "i1", i1_doc, 0);
	PyDict_SetItemString(dictionary, "i1", f);
	Py_DECREF(f);
	f = PyUFunc_FromFuncAndData(cephes1_functions, i1e_data, cephes_2_types, 2, 1, 1, PyUFunc_None, "i1e", i1e_doc, 0);
	PyDict_SetItemString(dictionary, "i1e", f);
	Py_DECREF(f);

	f = PyUFunc_FromFuncAndData(cephes2_functions, igamc_data, cephes_3_types, 2, 2, 1, PyUFunc_None, "gammaincc", gammaincc_doc, 0);
	PyDict_SetItemString(dictionary, "gammaincc", f);
	Py_DECREF(f);
	f = PyUFunc_FromFuncAndData(cephes2_functions, igam_data, cephes_3_types, 2, 2, 1, PyUFunc_None, "gammainc", gammainc_doc, 0);
	PyDict_SetItemString(dictionary, "gammainc", f);
	Py_DECREF(f);
	f = PyUFunc_FromFuncAndData(cephes2_functions, igami_data, cephes_3_types, 2, 2, 1, PyUFunc_None, "gammainccinv", gammainccinv_doc, 0);
	PyDict_SetItemString(dictionary, "gammainccinv", f);
	Py_DECREF(f);
	f = PyUFunc_FromFuncAndData(cephes2c_functions, iv_data, cephes_3c_types, 4, 2, 1, PyUFunc_None, "iv", iv_doc, 0);
	PyDict_SetItemString(dictionary, "iv", f);
	Py_DECREF(f);
	f = PyUFunc_FromFuncAndData(cephes2cp_functions, ive_data, cephes_3c_types, 4, 2, 1, PyUFunc_None, "ive", ive_doc, 0);
	PyDict_SetItemString(dictionary, "ive", f);
	Py_DECREF(f);

	f = PyUFunc_FromFuncAndData(cephes2_4_functions, ellpj_data, cephes_6_types, 2, 2, 4, PyUFunc_None, "ellipj", ellipj_doc, 0);
	PyDict_SetItemString(dictionary, "ellipj", f);
	Py_DECREF(f);

	f = PyUFunc_FromFuncAndData(cephes2a_functions, expn_data, cephes_3_types, 2, 2, 1, PyUFunc_None, "expn", expn_doc, 0);
	PyDict_SetItemString(dictionary, "expn", f);
	Py_DECREF(f);
	f = PyUFunc_FromFuncAndData(cephes1rc_functions, exp1_data, cephes_1rc_types, 4, 1, 1, PyUFunc_None, "exp1", exp1_doc, 0);
	PyDict_SetItemString(dictionary, "exp1", f);
	Py_DECREF(f);
	f = PyUFunc_FromFuncAndData(cephes1_functions, expi_data, cephes_2_types, 2, 1, 1, PyUFunc_None, "expi", expi_doc, 0);
	PyDict_SetItemString(dictionary, "expi", f);
	Py_DECREF(f);



	f = PyUFunc_FromFuncAndData(cephes2a_functions, jn_data, cephes_3_types, 2, 2, 1, PyUFunc_None, "jn", jn_doc, 0);
	PyDict_SetItemString(dictionary, "jn", f);
	Py_DECREF(f);
	f = PyUFunc_FromFuncAndData(cephes2a_functions, kn_data, cephes_3_types, 2, 2, 1, PyUFunc_None, "kn", kn_doc, 0);
	PyDict_SetItemString(dictionary, "kn", f);
	Py_DECREF(f);
	f = PyUFunc_FromFuncAndData(cephes2a_functions, pdtrc_data, cephes_3_types, 2, 2, 1, PyUFunc_None, "pdtrc", pdtrc_doc, 0);
	PyDict_SetItemString(dictionary, "pdtrc", f);
	Py_DECREF(f);
	f = PyUFunc_FromFuncAndData(cephes2a_functions, pdtr_data, cephes_3_types, 2, 2, 1, PyUFunc_None, "pdtr", pdtr_doc, 0);
	PyDict_SetItemString(dictionary, "pdtr", f);
	Py_DECREF(f);
	f = PyUFunc_FromFuncAndData(cephes2a_functions, pdtri_data, cephes_3_types, 2, 2, 1, PyUFunc_None, "pdtri", pdtri_doc, 0);
	PyDict_SetItemString(dictionary, "pdtri", f);
	Py_DECREF(f);
        /*  Use the student t library from cdflib (it supports doubles for
              degrees of freedom 
	f = PyUFunc_FromFuncAndData(cephes2a_functions, stdtr_data, cephes_3_types, 2, 2, 1, PyUFunc_None, "stdtr", stdtr_doc, 0);
	PyDict_SetItemString(dictionary, "stdtr", f);
	Py_DECREF(f);
	f = PyUFunc_FromFuncAndData(cephes2a_functions, stdtri_data, cephes_3_types, 2, 2, 1, PyUFunc_None, "stdtri", stdtri_doc, 0);
	PyDict_SetItemString(dictionary, "stdtri", f);
	Py_DECREF(f);
        */
	f = PyUFunc_FromFuncAndData(cephes2a_functions, yn_data, cephes_3_types, 2, 2, 1, PyUFunc_None, "yn", yn_doc, 0);
	PyDict_SetItemString(dictionary, "yn", f);
	Py_DECREF(f);
	f = PyUFunc_FromFuncAndData(cephes2a_functions, smirnov_data, cephes_3_types, 2, 2, 1, PyUFunc_None, "smirnov", smirnov_doc, 0);
	PyDict_SetItemString(dictionary, "smirnov", f);
	Py_DECREF(f);
	f = PyUFunc_FromFuncAndData(cephes2a_functions, smirnovi_data, cephes_3_types, 2, 2, 1, PyUFunc_None, "smirnovi", smirnovi_doc, 0);
	PyDict_SetItemString(dictionary, "smirnovi", f);
	Py_DECREF(f);

	f = PyUFunc_FromFuncAndData(cephes1c_4_functions, airy_data, cephes_5c_types, 4, 1, 4, PyUFunc_None, "airy", airy_doc, 0);
	PyDict_SetItemString(dictionary, "airy", f);
	Py_DECREF(f);

	f = PyUFunc_FromFuncAndData(cephes1c_4_functions, itairy_data, cephes_5_types, 2, 1, 4, PyUFunc_None, "itairy", itairy_doc, 0);
	PyDict_SetItemString(dictionary, "itairy", f);
	Py_DECREF(f);


	f = PyUFunc_FromFuncAndData(cephes1cp_4_functions, airye_data, cephes_5c_types, 4, 1, 4, PyUFunc_None, "airye", airye_doc, 0);
	PyDict_SetItemString(dictionary, "airye", f);
	Py_DECREF(f);

	f = PyUFunc_FromFuncAndData(cephes1_2_functions, fresnl_data, cephes_3_cmplx_types, 4, 1, 2, PyUFunc_None, "fresnel", fresnel_doc, 0);
	PyDict_SetItemString(dictionary, "fresnel", f);
	Py_DECREF(f);
	f = PyUFunc_FromFuncAndData(cephes1_2_functions, shichi_data, cephes_3_types, 2, 1, 2, PyUFunc_None, "shichi", shichi_doc, 0);
	PyDict_SetItemString(dictionary, "shichi", f);
	Py_DECREF(f);
	f = PyUFunc_FromFuncAndData(cephes1_2_functions, sici_data, cephes_3_types, 2, 1, 2, PyUFunc_None, "sici", sici_doc, 0);
	PyDict_SetItemString(dictionary, "sici", f);
	Py_DECREF(f);


	f = PyUFunc_FromFuncAndData(cephes1_2_functions, itj0y0_data, cephes_3_types, 2, 1, 2, PyUFunc_None, "itj0y0", itj0y0_doc, 0);
	PyDict_SetItemString(dictionary, "itj0y0", f);
	Py_DECREF(f);
	f = PyUFunc_FromFuncAndData(cephes1_2_functions, it2j0y0_data, cephes_3_types, 2, 1, 2, PyUFunc_None, "it2j0y0", it2j0y0_doc, 0);
	PyDict_SetItemString(dictionary, "it2j0y0", f);
	Py_DECREF(f);
	f = PyUFunc_FromFuncAndData(cephes1_2_functions, iti0k0_data, cephes_3_types, 2, 1, 2, PyUFunc_None, "iti0k0", iti0k0_doc, 0);
	PyDict_SetItemString(dictionary, "iti0k0", f);
	Py_DECREF(f);
	f = PyUFunc_FromFuncAndData(cephes1_2_functions, it2i0k0_data, cephes_3_types, 2, 1, 2, PyUFunc_None, "it2i0k0", it2i0k0_doc, 0);
	PyDict_SetItemString(dictionary, "it2i0k0", f);
	Py_DECREF(f);


	f = PyUFunc_FromFuncAndData(cephes1_functions, j0_data, cephes_2_types, 2, 1, 1, PyUFunc_None, "j0", j0_doc, 0);
	PyDict_SetItemString(dictionary, "j0", f);
	Py_DECREF(f);
	f = PyUFunc_FromFuncAndData(cephes1_functions, y0_data, cephes_2_types, 2, 1, 1, PyUFunc_None, "y0", y0_doc, 0);
	PyDict_SetItemString(dictionary, "y0", f);
	Py_DECREF(f);
	f = PyUFunc_FromFuncAndData(cephes1_functions, j1_data, cephes_2_types, 2, 1, 1, PyUFunc_None, "j1", j1_doc, 0);
	PyDict_SetItemString(dictionary, "j1", f);
	Py_DECREF(f);
	f = PyUFunc_FromFuncAndData(cephes1_functions, y1_data, cephes_2_types, 2, 1, 1, PyUFunc_None, "y1", y1_doc, 0);
	PyDict_SetItemString(dictionary, "y1", f);
	Py_DECREF(f);

	f = PyUFunc_FromFuncAndData(cephes2c_functions, jv_data, cephes_3c_types, 4, 2, 1, PyUFunc_None, "jv", jv_doc, 0);
	PyDict_SetItemString(dictionary, "jv", f);
	Py_DECREF(f);
	f = PyUFunc_FromFuncAndData(cephes2cp_functions, jve_data, cephes_3c_types, 4, 2, 1, PyUFunc_None, "jve", jve_doc, 0);
	PyDict_SetItemString(dictionary, "jve", f);
	Py_DECREF(f);
	f = PyUFunc_FromFuncAndData(cephes2c_functions, yv_data, cephes_3c_types, 4, 2, 1, PyUFunc_None, "yv", yv_doc, 0);
	PyDict_SetItemString(dictionary, "yv", f);
	Py_DECREF(f);
	f = PyUFunc_FromFuncAndData(cephes2cp_functions, yve_data, cephes_3c_types, 4, 2, 1, PyUFunc_None, "yve", yve_doc, 0);
	PyDict_SetItemString(dictionary, "yve", f);
	Py_DECREF(f);


	f = PyUFunc_FromFuncAndData(cephes1_functions, k0_data, cephes_2_types, 2, 1, 1, PyUFunc_None, "k0", k0_doc, 0);
	PyDict_SetItemString(dictionary, "k0", f);
	Py_DECREF(f);
	f = PyUFunc_FromFuncAndData(cephes1_functions, k0e_data, cephes_2_types, 2, 1, 1, PyUFunc_None, "k0e", k0e_doc, 0);
	PyDict_SetItemString(dictionary, "k0e", f);
	Py_DECREF(f);
	f = PyUFunc_FromFuncAndData(cephes1_functions, k1_data, cephes_2_types, 2, 1, 1, PyUFunc_None, "k1", k1_doc, 0);
	PyDict_SetItemString(dictionary, "k1", f);
	Py_DECREF(f);
	f = PyUFunc_FromFuncAndData(cephes1_functions, k1e_data, cephes_2_types, 2, 1, 1, PyUFunc_None, "k1e", k1e_doc, 0);
	PyDict_SetItemString(dictionary, "k1e", f);
	Py_DECREF(f);
	f = PyUFunc_FromFuncAndData(cephes2cp_functions, kv_data, cephes_3c_types, 4, 2, 1, PyUFunc_None, "kv", kv_doc, 0);
	PyDict_SetItemString(dictionary, "kv", f);
	Py_DECREF(f);
	f = PyUFunc_FromFuncAndData(cephes2cp_functions, kve_data, cephes_3c_types, 4, 2, 1, PyUFunc_None, "kve", kve_doc, 0);
	PyDict_SetItemString(dictionary, "kve", f);
	Py_DECREF(f);

	f = PyUFunc_FromFuncAndData(cephes2cpp_functions, hankel1_data, cephes_3cp_types, 2, 2, 1, PyUFunc_None, "hankel1", hankel1_doc, 0);
	PyDict_SetItemString(dictionary, "hankel1", f);
	Py_DECREF(f);
	f = PyUFunc_FromFuncAndData(cephes2cpp_functions, hankel1e_data, cephes_3cp_types, 2, 2, 1, PyUFunc_None, "hankel1e", hankel1e_doc, 0);
	PyDict_SetItemString(dictionary, "hankel1e", f);
	Py_DECREF(f);
	f = PyUFunc_FromFuncAndData(cephes2cpp_functions, hankel2_data, cephes_3cp_types, 2, 2, 1, PyUFunc_None, "hankel2", hankel2_doc, 0);
	PyDict_SetItemString(dictionary, "hankel2", f);
	Py_DECREF(f);
	f = PyUFunc_FromFuncAndData(cephes2cpp_functions, hankel2e_data, cephes_3cp_types, 2, 2, 1, PyUFunc_None, "hankel2e", hankel2e_doc, 0);
	PyDict_SetItemString(dictionary, "hankel2e", f);
	Py_DECREF(f);


	f = PyUFunc_FromFuncAndData(cephes1_functions, ndtr_data, cephes_2_types, 2, 1, 1, PyUFunc_None, "ndtr", ndtr_doc, 0);
	PyDict_SetItemString(dictionary, "ndtr", f);
	Py_DECREF(f);

	f = PyUFunc_FromFuncAndData(cephes1_functions, erfc_data, cephes_2_types, 2, 1, 1, PyUFunc_None, "erfc", erfc_doc, 0);
	PyDict_SetItemString(dictionary, "erfc", f);
	Py_DECREF(f);
	f = PyUFunc_FromFuncAndData(cephes1rc_functions, erf_data, cephes_1rc_types, 4, 1, 1, PyUFunc_None, "erf", erf_doc, 0);
	PyDict_SetItemString(dictionary, "erf", f);
	Py_DECREF(f);

	f = PyUFunc_FromFuncAndData(cephes1_functions, ndtri_data, cephes_2_types, 2, 1, 1, PyUFunc_None, "ndtri", ndtri_doc, 0);
	PyDict_SetItemString(dictionary, "ndtri", f);
	Py_DECREF(f);
	f = PyUFunc_FromFuncAndData(cephes1rc_functions, psi_data, cephes_1rc_types, 4, 1, 1, PyUFunc_None, "psi", psi_doc, 0);
	PyDict_SetItemString(dictionary, "psi", f);
	Py_DECREF(f);
	f = PyUFunc_FromFuncAndData(cephes1rc_functions, rgamma_data, cephes_1rc_types, 4, 1, 1, PyUFunc_None, "rgamma", rgamma_doc, 0);
	PyDict_SetItemString(dictionary, "rgamma", f);
	Py_DECREF(f);
	f = PyUFunc_FromFuncAndData(cephes1_functions, round_data, cephes_2_types, 2, 1, 1, PyUFunc_None, "round", round_doc, 0);
	PyDict_SetItemString(dictionary, "round", f);
	Py_DECREF(f);

	f = PyUFunc_FromFuncAndData(cephes1_functions, sindg_data, cephes_2_types, 2, 1, 1, PyUFunc_None, "sindg", sindg_doc, 0);
	PyDict_SetItemString(dictionary, "sindg", f);
	Py_DECREF(f);
	f = PyUFunc_FromFuncAndData(cephes1_functions, cosdg_data, cephes_2_types, 2, 1, 1, PyUFunc_None, "cosdg", cosdg_doc, 0);
	PyDict_SetItemString(dictionary, "cosdg", f);
	Py_DECREF(f);
	f = PyUFunc_FromFuncAndData(cephes3_functions, radian_data, cephes_4_types, 2, 3, 1, PyUFunc_None, "radian", radian_doc, 0);
	PyDict_SetItemString(dictionary, "radian", f);
	Py_DECREF(f);
	f = PyUFunc_FromFuncAndData(cephes1_functions, tandg_data, cephes_2_types, 2, 1, 1, PyUFunc_None, "tandg", tandg_doc, 0);
	PyDict_SetItemString(dictionary, "tandg", f);
	Py_DECREF(f);
	f = PyUFunc_FromFuncAndData(cephes1_functions, cotdg_data, cephes_2_types, 2, 1, 1, PyUFunc_None, "cotdg", cotdg_doc, 0);
	PyDict_SetItemString(dictionary, "cotdg", f);
	Py_DECREF(f);
	f = PyUFunc_FromFuncAndData(cephes1_functions, log1p_data, cephes_2_types, 2, 1, 1, PyUFunc_None, "log1p", log1p_doc, 0);
	PyDict_SetItemString(dictionary, "log1p", f);
	Py_DECREF(f);
	f = PyUFunc_FromFuncAndData(cephes1_functions, expm1_data, cephes_2_types, 2, 1, 1, PyUFunc_None, "expm1", expm1_doc, 0);
	PyDict_SetItemString(dictionary, "expm1", f);
	Py_DECREF(f);
	f = PyUFunc_FromFuncAndData(cephes1_functions, cosm1_data, cephes_2_types, 2, 1, 1, PyUFunc_None, "cosm1", cosm1_doc, 0);
	PyDict_SetItemString(dictionary, "cosm1", f);
	Py_DECREF(f);

	f = PyUFunc_FromFuncAndData(cephes1_functions, spence_data, cephes_2_types, 2, 1, 1, PyUFunc_None, "spence", spence_doc, 0);
	PyDict_SetItemString(dictionary, "spence", f);
	Py_DECREF(f);
	f = PyUFunc_FromFuncAndData(cephes1_functions, zetac_data, cephes_2_types, 2, 1, 1, PyUFunc_None, "zetac", zetac_doc, 0);
	PyDict_SetItemString(dictionary, "zetac", f);
	Py_DECREF(f);

	f = PyUFunc_FromFuncAndData(cephes2_functions, struve_data, cephes_3_types, 2, 2, 1, PyUFunc_None, "struve", struve_doc, 0);
	PyDict_SetItemString(dictionary, "struve", f);
	Py_DECREF(f);
	f = PyUFunc_FromFuncAndData(cephes2_functions, modstruve_data, cephes_3_types, 2, 2, 1, PyUFunc_None, "modstruve", modstruve_doc, 0);
	PyDict_SetItemString(dictionary, "modstruve", f);
	Py_DECREF(f);
	f = PyUFunc_FromFuncAndData(cephes1_functions, itstruve0_data, cephes_2_types, 2, 1, 1, PyUFunc_None, "itstruve0", itstruve0_doc, 0);
	PyDict_SetItemString(dictionary, "itstruve0", f);
	Py_DECREF(f);
	f = PyUFunc_FromFuncAndData(cephes1_functions, it2struve0_data, cephes_2_types, 2, 1, 1, PyUFunc_None, "it2struve0", it2struve0_doc, 0);
	PyDict_SetItemString(dictionary, "it2struve0", f);
	Py_DECREF(f);
	f = PyUFunc_FromFuncAndData(cephes1_functions, itmodstruve0_data, cephes_2_types, 2, 1, 1, PyUFunc_None, "itmodstruve0", itmodstruve0_doc, 0);
	PyDict_SetItemString(dictionary, "itmodstruve0", f);
	Py_DECREF(f);

	f = PyUFunc_FromFuncAndData(cephes1cpb_4_functions, kelvin_data, cephes_5b2_types, 2, 1, 4, PyUFunc_None, "kelvin", kelvin_doc, 0);
	PyDict_SetItemString(dictionary, "kelvin", f);
	Py_DECREF(f);
	f = PyUFunc_FromFuncAndData(cephes1_functions, ber_data, cephes_2_types, 2, 1, 1, PyUFunc_None, "ber", ber_doc, 0);
	PyDict_SetItemString(dictionary, "ber", f);
	Py_DECREF(f);
	f = PyUFunc_FromFuncAndData(cephes1_functions, bei_data, cephes_2_types, 2, 1, 1, PyUFunc_None, "bei", bei_doc, 0);
	PyDict_SetItemString(dictionary, "bei", f);
	Py_DECREF(f);
	f = PyUFunc_FromFuncAndData(cephes1_functions, ker_data, cephes_2_types, 2, 1, 1, PyUFunc_None, "ker", ker_doc, 0);
	PyDict_SetItemString(dictionary, "ker", f);
	Py_DECREF(f);
	f = PyUFunc_FromFuncAndData(cephes1_functions, kei_data, cephes_2_types, 2, 1, 1, PyUFunc_None, "kei", kei_doc, 0);
	PyDict_SetItemString(dictionary, "kei", f);
	Py_DECREF(f);
	f = PyUFunc_FromFuncAndData(cephes1_functions, berp_data, cephes_2_types, 2, 1, 1, PyUFunc_None, "berp", berp_doc, 0);
	PyDict_SetItemString(dictionary, "berp", f);
	Py_DECREF(f);
	f = PyUFunc_FromFuncAndData(cephes1_functions, beip_data, cephes_2_types, 2, 1, 1, PyUFunc_None, "beip", beip_doc, 0);
	PyDict_SetItemString(dictionary, "beip", f);
	Py_DECREF(f);
	f = PyUFunc_FromFuncAndData(cephes1_functions, kerp_data, cephes_2_types, 2, 1, 1, PyUFunc_None, "kerp", kerp_doc, 0);
	PyDict_SetItemString(dictionary, "kerp", f);
	Py_DECREF(f);
	f = PyUFunc_FromFuncAndData(cephes1_functions, keip_data, cephes_2_types, 2, 1, 1, PyUFunc_None, "keip", keip_doc, 0);
	PyDict_SetItemString(dictionary, "keip", f);
	Py_DECREF(f);


	f = PyUFunc_FromFuncAndData(cephes2_functions, zeta_data, cephes_3_types, 2, 2, 1, PyUFunc_None, "zeta", zeta_doc, 0);
	PyDict_SetItemString(dictionary, "zeta", f);
	Py_DECREF(f);

	f = PyUFunc_FromFuncAndData(cephes1_functions, kolmogorov_data, cephes_2_types, 2, 1, 1, PyUFunc_None, "kolmogorov", kolmogorov_doc, 0);
	PyDict_SetItemString(dictionary, "kolmogorov", f);
	Py_DECREF(f);

	f = PyUFunc_FromFuncAndData(cephes1_functions, kolmogi_data, cephes_2_types, 2, 1, 1, PyUFunc_None, "kolmogi", kolmogi_doc, 0);
	PyDict_SetItemString(dictionary, "kolmogi", f);
	Py_DECREF(f);

	f = PyUFunc_FromFuncAndData(cephes1c_functions, wofz_data, cephes_1c_types, 2, 1, 1, PyUFunc_None, "wofz", wofz_doc, 0); 
	PyDict_SetItemString(dictionary, "wofz", f);
	Py_DECREF(f);

	f = PyUFunc_FromFuncAndData(cephes3_functions, besselpoly_data, cephes_4_types, 2, 3, 1, PyUFunc_None, "besselpoly", besselpoly_doc, 0);
	PyDict_SetItemString(dictionary, "besselpoly", f);
	Py_DECREF(f);

	f = PyUFunc_FromFuncAndData(cephes3_functions, cdfbet3_data, cephes_4_types, 2, 3, 1, PyUFunc_None, "btdtria", "", 0);
	PyDict_SetItemString(dictionary, "btdtria", f);
	Py_DECREF(f);
	f = PyUFunc_FromFuncAndData(cephes3_functions, cdfbet4_data, cephes_4_types, 2, 3, 1, PyUFunc_None, "btdtrib", "", 0);
	PyDict_SetItemString(dictionary, "btdtrib", f);
	Py_DECREF(f);

	f = PyUFunc_FromFuncAndData(cephes3_functions, cdfbin2_data, cephes_4_types, 2, 3, 1, PyUFunc_None, "bdtrik", "", 0);
	PyDict_SetItemString(dictionary, "bdtrik", f);
	Py_DECREF(f);
	f = PyUFunc_FromFuncAndData(cephes3_functions, cdfbin3_data, cephes_4_types, 2, 3, 1, PyUFunc_None, "bdtrin", "", 0);
	PyDict_SetItemString(dictionary, "bdtrin", f);
	Py_DECREF(f);

	f = PyUFunc_FromFuncAndData(cephes2_functions, cdfchi3_data, cephes_3_types, 2, 2, 1, PyUFunc_None, "chdtriv", "", 0);
	PyDict_SetItemString(dictionary, "chdtriv", f);
	Py_DECREF(f);


	f = PyUFunc_FromFuncAndData(cephes3_functions, cdfchn1_data, cephes_4_types, 2, 3, 1, PyUFunc_None, "chndtr", "", 0);
	PyDict_SetItemString(dictionary, "chndtr", f);
	Py_DECREF(f);
	f = PyUFunc_FromFuncAndData(cephes3_functions, cdfchn2_data, cephes_4_types, 2, 3, 1, PyUFunc_None, "chndtrix", "", 0);
	PyDict_SetItemString(dictionary, "chndtrix", f);
	Py_DECREF(f);
	f = PyUFunc_FromFuncAndData(cephes3_functions, cdfchn3_data, cephes_4_types, 2, 3, 1, PyUFunc_None, "chndtridf", "", 0);
	PyDict_SetItemString(dictionary, "chndtridf", f);
	Py_DECREF(f);
	f = PyUFunc_FromFuncAndData(cephes3_functions, cdfchn4_data, cephes_4_types, 2, 3, 1, PyUFunc_None, "chndtrinc", "", 0);
	PyDict_SetItemString(dictionary, "chndtrinc", f);
	Py_DECREF(f);

        /*
        f = PyUFunc_FromFuncAndData(cephes3_functions, cdff1_data, cephes_4_types, 2, 3, 1, PyUFunc_None, "fdtr", fdtr_doc, 0);
        PyDict_SetItemString(dictionary, "fdtr", f);
	Py_DECREF(f);
	f = PyUFunc_FromFuncAndData(cephes3_functions, cdff2_data, cephes_4_types, 2, 3, 1, PyUFunc_None, "fdtrix", fdtri_doc, 0);
	PyDict_SetItemString(dictionary, "fdtrix", f);
 	Py_DECREF(f);
        */
        
        /*  The Fortran code for this one seems not to be working properly.
	f = PyUFunc_FromFuncAndData(cephes3_functions, cdff3_data, cephes_4_types, 2, 3, 1, PyUFunc_None, "fdtridfn", "", 0);
	PyDict_SetItemString(dictionary, "fdtridfn", f);
	Py_DECREF(f);
        */ 
	f = PyUFunc_FromFuncAndData(cephes3_functions, cdff4_data, cephes_4_types, 2, 3, 1, PyUFunc_None, "fdtridfd", "", 0);
	PyDict_SetItemString(dictionary, "fdtridfd", f);
	Py_DECREF(f);

	f = PyUFunc_FromFuncAndData(cephes4_functions, cdffnc1_data, cephes_5_types, 2, 4, 1, PyUFunc_None, "ncfdtr", "", 0);
	PyDict_SetItemString(dictionary, "ncfdtr", f);
	Py_DECREF(f);
	f = PyUFunc_FromFuncAndData(cephes4_functions, cdffnc2_data, cephes_5_types, 2, 4, 1, PyUFunc_None, "ncfdtri", "", 0);
	PyDict_SetItemString(dictionary, "ncfdtri", f);
	Py_DECREF(f);
	f = PyUFunc_FromFuncAndData(cephes4_functions, cdffnc3_data, cephes_5_types, 2, 4, 1, PyUFunc_None, "ncfdtridfn", "", 0);
	PyDict_SetItemString(dictionary, "ncfdtridfn", f);
	Py_DECREF(f);

	f = PyUFunc_FromFuncAndData(cephes4_functions, cdffnc4_data, cephes_5_types, 2, 4, 1, PyUFunc_None, "ncfdtridfd", "", 0);
	PyDict_SetItemString(dictionary, "ncfdtridfd", f);
	Py_DECREF(f);
	f = PyUFunc_FromFuncAndData(cephes4_functions, cdffnc5_data, cephes_5_types, 2, 4, 1, PyUFunc_None, "ncfdtrinc", "", 0);
	PyDict_SetItemString(dictionary, "ncfdtrinc", f);
	Py_DECREF(f);

        /*
	f = PyUFunc_FromFuncAndData(cephes3_functions, cdfgam1_data, cephes_4_types, 2, 3, 1, PyUFunc_None, "gdtr2", "", 0);
	PyDict_SetItemString(dictionary, "gdtr2", f);
	Py_DECREF(f);
        */
	f = PyUFunc_FromFuncAndData(cephes3_functions, cdfgam2_data, cephes_4_types, 2, 3, 1, PyUFunc_None, "gdtrix", "", 0);
	PyDict_SetItemString(dictionary, "gdtrix", f);
	Py_DECREF(f);
	f = PyUFunc_FromFuncAndData(cephes3_functions, cdfgam3_data, cephes_4_types, 2, 3, 1, PyUFunc_None, "gdtrib", "", 0);
	PyDict_SetItemString(dictionary, "gdtrib", f);
	Py_DECREF(f);
	f = PyUFunc_FromFuncAndData(cephes3_functions, cdfgam4_data, cephes_4_types, 2, 3, 1, PyUFunc_None, "gdtria", "", 0);
	PyDict_SetItemString(dictionary, "gdtria", f);
	Py_DECREF(f);

	f = PyUFunc_FromFuncAndData(cephes3_functions, cdfnbn2_data, cephes_4_types, 2, 3, 1, PyUFunc_None, "nbdtrik", "", 0);
	PyDict_SetItemString(dictionary, "nbdtrik", f);
	Py_DECREF(f);
	f = PyUFunc_FromFuncAndData(cephes3_functions, cdfnbn3_data, cephes_4_types, 2, 3, 1, PyUFunc_None, "nbdtrin", "", 0);
	PyDict_SetItemString(dictionary, "nbdtrin", f);
	Py_DECREF(f);

	f = PyUFunc_FromFuncAndData(cephes3_functions, cdfnor3_data, cephes_4_types, 2, 3, 1, PyUFunc_None, "nrdtrimn", "", 0);
	PyDict_SetItemString(dictionary, "nrdtrimn", f);
	Py_DECREF(f);
	f = PyUFunc_FromFuncAndData(cephes3_functions, cdfnor4_data, cephes_4_types, 2, 3, 1, PyUFunc_None, "nrdtrisd", "", 0);
	PyDict_SetItemString(dictionary, "nrdtrisd", f);
	Py_DECREF(f);

	f = PyUFunc_FromFuncAndData(cephes2_functions, cdfpoi2_data, cephes_3_types, 2, 2, 1, PyUFunc_None, "pdtrik", "", 0);
	PyDict_SetItemString(dictionary, "pdtrik", f);
	Py_DECREF(f);

	f = PyUFunc_FromFuncAndData(cephes2_functions, cdft1_data, cephes_3_types, 2, 2, 1, PyUFunc_None, "stdtr", stdtr_doc, 0);
	PyDict_SetItemString(dictionary, "stdtr", f);
	Py_DECREF(f);
	f = PyUFunc_FromFuncAndData(cephes2_functions, cdft2_data, cephes_3_types, 2, 2, 1, PyUFunc_None, "stdtrit", stdtrit_doc, 0);
	PyDict_SetItemString(dictionary, "stdtrit", f);
	Py_DECREF(f);
	f = PyUFunc_FromFuncAndData(cephes2_functions, cdft3_data, cephes_3_types, 2, 2, 1, PyUFunc_None, "stdtridf", stdtridf_doc, 0);
	PyDict_SetItemString(dictionary, "stdtridf", f);
	Py_DECREF(f);

	f = PyUFunc_FromFuncAndData(cephes3_functions, cdftnc1_data, cephes_4_types, 2, 3, 1, PyUFunc_None, "nctdtr", "", 0);
	PyDict_SetItemString(dictionary, "nctdtr", f);
	Py_DECREF(f);
	f = PyUFunc_FromFuncAndData(cephes3_functions, cdftnc2_data, cephes_4_types, 2, 3, 1, PyUFunc_None, "nctdtrit", "", 0);
	PyDict_SetItemString(dictionary, "nctdtrit", f);
	Py_DECREF(f);
	f = PyUFunc_FromFuncAndData(cephes3_functions, cdftnc3_data, cephes_4_types, 2, 3, 1, PyUFunc_None, "nctdtridf", "", 0);
	PyDict_SetItemString(dictionary, "nctdtridf", f);
	Py_DECREF(f);
	f = PyUFunc_FromFuncAndData(cephes3_functions, cdftnc4_data, cephes_4_types, 2, 3, 1, PyUFunc_None, "nctdtrinc", "", 0);
	PyDict_SetItemString(dictionary, "nctdtrinc", f);
	Py_DECREF(f);

	f = PyUFunc_FromFuncAndData(cephes2_functions, tklambda_data, cephes_3_types, 2, 2, 1, PyUFunc_None, "tklmbda", "", 0);
	PyDict_SetItemString(dictionary, "tklmbda", f);
	Py_DECREF(f);

        
	f = PyUFunc_FromFuncAndData(cephes2_functions, mathieu_a_data, cephes_3_types, 2, 2, 1, PyUFunc_None, "mathieu_a", mathieu_a_doc, 0);
	PyDict_SetItemString(dictionary, "mathieu_a", f);
	Py_DECREF(f);
	f = PyUFunc_FromFuncAndData(cephes2_functions, mathieu_b_data, cephes_3_types, 2, 2, 1, PyUFunc_None, "mathieu_b", mathieu_b_doc, 0);
	PyDict_SetItemString(dictionary, "mathieu_b", f);
	Py_DECREF(f);
	f = PyUFunc_FromFuncAndData(cephes3_2_functions, mathieu_cem_data, cephes_5_types, 2, 3, 2, PyUFunc_None, "mathieu_cem", mathieu_cem_doc, 0);
	PyDict_SetItemString(dictionary, "mathieu_cem", f);
	Py_DECREF(f);
	f = PyUFunc_FromFuncAndData(cephes3_2_functions, mathieu_sem_data, cephes_5_types, 2, 3, 2, PyUFunc_None, "mathieu_sem", mathieu_sem_doc, 0);
	PyDict_SetItemString(dictionary, "mathieu_sem", f);
	Py_DECREF(f);
	f = PyUFunc_FromFuncAndData(cephes3_2_functions, mathieu_mcem1_data, cephes_5_types, 2, 3, 2, PyUFunc_None, "mathieu_modcem1", mathieu_modcem1_doc, 0);
	PyDict_SetItemString(dictionary, "mathieu_modcem1", f);
	Py_DECREF(f);
	f = PyUFunc_FromFuncAndData(cephes3_2_functions, mathieu_mcem2_data, cephes_5_types, 2, 3, 2, PyUFunc_None, "mathieu_modcem2", mathieu_modcem2_doc, 0);
	PyDict_SetItemString(dictionary, "mathieu_modcem2", f);
	Py_DECREF(f);
	f = PyUFunc_FromFuncAndData(cephes3_2_functions, mathieu_msem1_data, cephes_5_types, 2, 3, 2, PyUFunc_None, "mathieu_modsem1", mathieu_modsem1_doc, 0);
	PyDict_SetItemString(dictionary, "mathieu_modsem1", f);
	Py_DECREF(f);
	f = PyUFunc_FromFuncAndData(cephes3_2_functions, mathieu_msem2_data, cephes_5_types, 2, 3, 2, PyUFunc_None, "mathieu_modsem2", mathieu_modsem2_doc, 0);
	PyDict_SetItemString(dictionary, "mathieu_modsem2", f);
	Py_DECREF(f);

	f = PyUFunc_FromFuncAndData(cephes3_functions, lpmv_data, cephes_4_types, 2, 3, 1, PyUFunc_None, "lpmv", lpmv_doc, 0);
	PyDict_SetItemString(dictionary, "lpmv", f);
	Py_DECREF(f);

	f = PyUFunc_FromFuncAndData(cephes2_2_functions, pbwa_data, cephes_4_types, 2, 2, 2, PyUFunc_None, "pbwa", pbwa_doc, 0);
	PyDict_SetItemString(dictionary, "pbwa", f);
	Py_DECREF(f);
	f = PyUFunc_FromFuncAndData(cephes2_2_functions, pbdv_data, cephes_4_types, 2, 2, 2, PyUFunc_None, "pbdv", pbdv_doc, 0);
	PyDict_SetItemString(dictionary, "pbdv", f);
	Py_DECREF(f);
	f = PyUFunc_FromFuncAndData(cephes2_2_functions, pbvv_data, cephes_4_types, 2, 2, 2, PyUFunc_None, "pbvv", pbvv_doc, 0);
	PyDict_SetItemString(dictionary, "pbvv", f);
	Py_DECREF(f);

	f = PyUFunc_FromFuncAndData(cephes3_functions, prolate_segv_data, cephes_4_types, 2, 3, 1, PyUFunc_None, "pro_cv", pro_cv_doc, 0);
	PyDict_SetItemString(dictionary, "pro_cv", f);
	Py_DECREF(f);
	f = PyUFunc_FromFuncAndData(cephes3_functions, oblate_segv_data, cephes_4_types, 2, 3, 1, PyUFunc_None, "obl_cv", obl_cv_doc, 0);
	PyDict_SetItemString(dictionary, "obl_cv", f);
	Py_DECREF(f);
	f = PyUFunc_FromFuncAndData(cephes5_2_functions, prolate_aswfa_data, cephes_7_types, 2, 5, 2, PyUFunc_None, "pro_ang1_cv", pro_ang1_cv_doc, 0);
	PyDict_SetItemString(dictionary, "pro_ang1_cv", f);
	Py_DECREF(f);
	f = PyUFunc_FromFuncAndData(cephes5_2_functions, prolate_radial1_data, cephes_7_types, 2, 5, 2, PyUFunc_None, "pro_rad1_cv", pro_rad1_cv_doc, 0);
	PyDict_SetItemString(dictionary, "pro_rad1_cv", f);
	Py_DECREF(f);
	f = PyUFunc_FromFuncAndData(cephes5_2_functions, prolate_radial2_data, cephes_7_types, 2, 5, 2, PyUFunc_None, "pro_rad2_cv", pro_rad2_cv_doc, 0);
	PyDict_SetItemString(dictionary, "pro_rad2_cv", f);
	Py_DECREF(f);
	f = PyUFunc_FromFuncAndData(cephes5_2_functions, oblate_aswfa_data, cephes_7_types, 2, 5, 2, PyUFunc_None, "obl_ang1_cv", obl_ang1_cv_doc, 0);
	PyDict_SetItemString(dictionary, "obl_ang1_cv", f);
	Py_DECREF(f);
	f = PyUFunc_FromFuncAndData(cephes5_2_functions, oblate_radial1_data, cephes_7_types, 2, 5, 2, PyUFunc_None, "obl_rad1_cv", obl_rad1_cv_doc, 0);
	PyDict_SetItemString(dictionary, "obl_rad1_cv", f);
	Py_DECREF(f);
	f = PyUFunc_FromFuncAndData(cephes5_2_functions, oblate_radial2_data, cephes_7_types, 2, 5, 2, PyUFunc_None, "obl_rad2_cv", obl_rad2_cv_doc, 0);
	PyDict_SetItemString(dictionary, "obl_rad2_cv", f);
	Py_DECREF(f);
	f = PyUFunc_FromFuncAndData(cephes4_2_functions, prolate_aswfa_nocv_data, cephes_6_types, 2, 4, 2, PyUFunc_None, "pro_ang1", pro_ang1_doc, 0);
	PyDict_SetItemString(dictionary, "pro_ang1", f);
	Py_DECREF(f);
	f = PyUFunc_FromFuncAndData(cephes4_2_functions, prolate_radial1_nocv_data, cephes_6_types, 2, 4, 2, PyUFunc_None, "pro_rad1", pro_rad1_doc, 0);
	PyDict_SetItemString(dictionary, "pro_rad1", f);
	Py_DECREF(f);
	f = PyUFunc_FromFuncAndData(cephes4_2_functions, prolate_radial2_nocv_data, cephes_6_types, 2, 4, 2, PyUFunc_None, "pro_rad2", pro_rad2_doc, 0);
	PyDict_SetItemString(dictionary, "pro_rad2", f);
	Py_DECREF(f);
	f = PyUFunc_FromFuncAndData(cephes4_2_functions, oblate_aswfa_nocv_data, cephes_6_types, 2, 4, 2, PyUFunc_None, "obl_ang1", obl_ang1_doc, 0);
	PyDict_SetItemString(dictionary, "obl_ang1", f);
	Py_DECREF(f);
	f = PyUFunc_FromFuncAndData(cephes4_2_functions, oblate_radial1_nocv_data, cephes_6_types, 2, 4, 2, PyUFunc_None, "obl_rad1", obl_rad1_doc, 0);
	PyDict_SetItemString(dictionary, "obl_rad1", f);
	Py_DECREF(f);
	f = PyUFunc_FromFuncAndData(cephes4_2_functions, oblate_radial2_nocv_data, cephes_6_types, 2, 4, 2, PyUFunc_None, "obl_rad2", obl_rad2_doc, 0);
	PyDict_SetItemString(dictionary, "obl_rad2", f);
	Py_DECREF(f);



	f = PyUFunc_FromFuncAndData(cephes1_2c_functions, modfresnelp_data, cephes_3cp_types, 2, 1, 2, PyUFunc_None, "modfresnelp", modfresnelp_doc, 0);
	PyDict_SetItemString(dictionary, "modfresnelp", f);
	Py_DECREF(f);

	f = PyUFunc_FromFuncAndData(cephes1_2c_functions, modfresnelm_data, cephes_3cp_types, 2, 1, 2, PyUFunc_None, "modfresnelm", modfresnelm_doc, 0);
	PyDict_SetItemString(dictionary, "modfresnelm", f);
	Py_DECREF(f);



}

static char errprint_doc[] = \
"errprint({flag}) sets the error printing flag for special functions\n" \
"    (from the cephesmodule). The output is the previous state.\n" \
"    With errprint(0) no error messages are shown;\n" \
"    the default is errprint(1).\n" \
"    If no argument is given the current state of\n" \
"    the flag is returned and no change occurs.\n";


static PyObject *errprint_func(PyObject *self, PyObject *args)
{
  int inflag = -37;
  int oldflag = 0;
  if (!PyArg_ParseTuple ( args, "|i;cephes.errprint", &inflag)) return NULL;

  oldflag = scipy_special_print_error_messages;  
  if (inflag != -37) {
    scipy_special_print_error_messages = (inflag != 0);
  }
  return PyInt_FromLong((long) oldflag);
}

  
static struct PyMethodDef methods[] = {
  {"errprint", errprint_func, METH_VARARGS, errprint_doc},
  {NULL,		NULL, 0}		/* sentinel */
};


void initcephes(void) {
  PyObject *m, *d, *s;
  
  /* Create the module and add the functions */
  m = Py_InitModule("cephes", methods); 

  /* Import the ufunc objects */
  import_array();
  import_ufunc();

  /* Add some symbolic constants to the module */
  d = PyModule_GetDict(m);

  s = PyString_FromString("2.0");
  PyDict_SetItemString(d, "__version__", s);
  Py_DECREF(s);

  /* Add scipy_special_print_error_message global variable */
  /*  No, instead acessible through errprint */

  /* Load the cephes operators into the array module's namespace */
  Cephes_InitOperators(d); 
  
  /* Check for errors */
  if (PyErr_Occurred())
    Py_FatalError("can't initialize module cephes");
}

