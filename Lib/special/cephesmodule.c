
/* Cephes module version 1.3
 *  This module defines the functions in the cephes and amos libraries as
 *   Numerical python ufunc objects so that they can operate on arbitrary 
 *   NumPy arrays with broadcasting and typecasting rules implemented.
 *  
 *  Copyright 1999  Travis E. Oliphant
 *  This program may be freely modified and distributed under the conditions of the 
 *  LGPL, provided this notification remain.  No warranty is implied.
 *  USE at your own risk.
 */

#include "Python.h"
#include "Numeric/arrayobject.h"
#include "Numeric/ufuncobject.h" 
#include "ufunc_extras.h"
#include "abstract.h"
#include "cephes.h"
#include "amos_wrappers.h"
#include "toms_wrappers.h"
#include "c_misc/misc.h"
#ifdef macintosh
#include "mymath.h"
#else
#include <math.h>
#endif

#include "cephes_doc.h"

#define ABS(x) (y=(x); y < 0 ? -y : y)

/* isnan and isinf and isfinite functions */
static void FLOAT_isnan(char **args, int *dimensions, int *steps, void *func) {
    int i, is1=steps[0], os=steps[1], n=dimensions[0];
    char *i1=args[0], *op=args[1];
    for (i=0; i < n; i++, i1+=is1, op+=os) {
	*((signed char *)op) = (signed char) ABS(isnan((double)(*((float *)i1))));
    }
}

static void DOUBLE_isnan(char **args, int *dimensions, int *steps, void *func) {
    int i, is1=steps[0], os=steps[1], n=dimensions[0];
    char *i1=args[0], *op=args[1];
    for (i=0; i < n; i++, i1+=is1, op+=os) {
	*((signed char *)op) = (signed char) ABS(isnan((double)(*((double *)i1))));
    }
}

static void CFLOAT_isnan(char **args, int *dimensions, int *steps, void *func) {
    int i, is1=steps[0], os=steps[1], n=dimensions[0];
    char *i1=args[0], *op=args[1];
    for (i=0; i < n; i++, i1+=is1, op+=os) {
	*((signed char *)op) = (signed char) isnan((double)((float *)i1)[0]) || isnan((double)((float *)i1)[1]);
    }
}

static void CDOUBLE_isnan(char **args, int *dimensions, int *steps, void *func) {
    int i, is1=steps[0], os=steps[1], n=dimensions[0];
    char *i1=args[0], *op=args[1];
    for (i=0; i < n; i++, i1+=is1, op+=os) {
	*((signed char *)op) = (signed char) isnan((double)((double *)i1)[0]) || isnan((double)((double *)i1)[1]);
    }
}


static void FLOAT_isinf(char **args, int *dimensions, int *steps, void *func) {
    int i, is1=steps[0], os=steps[1], n=dimensions[0];
    char *i1=args[0], *op=args[1];
    for (i=0; i < n; i++, i1+=is1, op+=os) {
	*((signed char *)op) = (signed char) !(isfinite((double)(*((float *)i1))) || isnan((double)(*((float *)i1))));
    }
}

static void DOUBLE_isinf(char **args, int *dimensions, int *steps, void *func) {
    int i, is1=steps[0], os=steps[1], n=dimensions[0];
    char *i1=args[0], *op=args[1];
    for (i=0; i < n; i++, i1+=is1, op+=os) {
	*((signed char *)op)= (signed char) !(isfinite((double)(*((double *)i1))) || isnan((double)(*((double *)i1))));
    }
}

static void CFLOAT_isinf(char **args, int *dimensions, int *steps, void *func) {
    int i, is1=steps[0], os=steps[1], n=dimensions[0];
    char *i1=args[0], *op=args[1];
    for (i=0; i < n; i++, i1+=is1, op+=os) {
	*((signed char *)op)= (signed char) !((isfinite((double)(((float *)i1)[0])) && isfinite((double)(((float *)i1)[1]))) || isnan((double)(((float *)i1)[0])) || isnan((double)(((float *)i1)[1])));
    }
}

static void CDOUBLE_isinf(char **args, int *dimensions, int *steps, void *func) {
    int i, is1=steps[0], os=steps[1], n=dimensions[0];
    char *i1=args[0], *op=args[1];
    for (i=0; i < n; i++, i1+=is1, op+=os) {
	*((signed char *)op)= (signed char) !((isfinite((double)(((double *)i1)[0])) && isfinite((double)(((double *)i1)[1]))) || isnan((double)(((double *)i1)[0])) || isnan((double)(((double *)i1)[1])));
    }
}


static void FLOAT_isfinite(char **args, int *dimensions, int *steps, void *func) {
    int i, is1=steps[0], os=steps[1], n=dimensions[0];
    char *i1=args[0], *op=args[1];
    for (i=0; i < n; i++, i1+=is1, op+=os) {
	*((signed char *)op) = (signed char) isfinite((double)(*((float *)i1)));
    }
}

static void DOUBLE_isfinite(char **args, int *dimensions, int *steps, void *func) {
    int i, is1=steps[0], os=steps[1], n=dimensions[0];
    char *i1=args[0], *op=args[1];
    for (i=0; i < n; i++, i1+=is1, op+=os) {
	*((signed char *)op) = (signed char) isfinite((double)(*((double *)i1)));
    }
}

static void CFLOAT_isfinite(char **args, int *dimensions, int *steps, void *func) {
    int i, is1=steps[0], os=steps[1], n=dimensions[0];
    char *i1=args[0], *op=args[1];
    for (i=0; i < n; i++, i1+=is1, op+=os) {
	*((signed char *)op) = (signed char) isfinite((double)((float *)i1)[0]) && isfinite((double)((float *)i1)[1]);
    }
}

static void CDOUBLE_isfinite(char **args, int *dimensions, int *steps, void *func) {
    int i, is1=steps[0], os=steps[1], n=dimensions[0];
    char *i1=args[0], *op=args[1];
    for (i=0; i < n; i++, i1+=is1, op+=os) {
	*((signed char *)op) = (signed char) isfinite((double)((double *)i1)[0]) && isfinite((double)((double *)i1)[1]);
    }
}

static PyUFuncGenericFunction isnan_functions[] = {FLOAT_isnan, DOUBLE_isnan, CFLOAT_isnan, CDOUBLE_isnan, NULL};
static PyUFuncGenericFunction isinf_functions[] = {FLOAT_isinf, DOUBLE_isinf, CFLOAT_isinf, CDOUBLE_isinf, NULL};
static PyUFuncGenericFunction isfinite_functions[] = {FLOAT_isfinite, DOUBLE_isfinite, CFLOAT_isfinite, CDOUBLE_isfinite, NULL};

static char isinf_signatures[] = { PyArray_FLOAT, PyArray_SBYTE, PyArray_DOUBLE, PyArray_SBYTE, PyArray_CFLOAT, PyArray_SBYTE,  PyArray_CDOUBLE, PyArray_SBYTE, };

static void * isnan_data[] = {(void *)NULL, (void *)NULL, (void *)NULL, (void *)NULL};
static void * isinf_data[] = {(void *)NULL, (void *)NULL, (void *)NULL, (void *)NULL};
static void * isfinite_data[] = {(void *)NULL, (void *)NULL, (void *)NULL, (void *)NULL};

static PyUFuncGenericFunction cephes1_functions[] = { NULL, NULL, };
static PyUFuncGenericFunction cephes1_2_functions[] = { NULL, NULL, };
static PyUFuncGenericFunction cephes1c_4_functions[] = { NULL, NULL, NULL, NULL };
static PyUFuncGenericFunction cephes1cp_4_functions[] = { NULL, NULL, NULL, NULL};
static PyUFuncGenericFunction cephes2_functions[] = { NULL, NULL, };
static PyUFuncGenericFunction cephes2_2_functions[] = { NULL, NULL, };
static PyUFuncGenericFunction cephes2_4_functions[] = { NULL, NULL, };
static PyUFuncGenericFunction cephes2a_functions[] = { NULL, NULL, };
static PyUFuncGenericFunction cephes2c_functions[] = { NULL, NULL, NULL, NULL };
static PyUFuncGenericFunction cephes2cp_functions[] = { NULL, NULL, NULL, NULL, };
static PyUFuncGenericFunction cephes2cpp_functions[] = { NULL, NULL, };
static PyUFuncGenericFunction cephes3_functions[] = { NULL, NULL, };
static PyUFuncGenericFunction cephes3a_functions[] = { NULL, NULL, };
static PyUFuncGenericFunction cephes4_functions[] = { NULL, NULL, };
static PyUFuncGenericFunction cephes4a_2_functions[] = { NULL, NULL, };
static PyUFuncGenericFunction cephes4_2_functions[] = { NULL, NULL, };

static PyUFuncGenericFunction cephes1c_functions[] = { NULL, };

static void * airy_data[] = { (void *)airy, (void *)airy, (void *)cairy_wrap, (void *)cairy_wrap,};
static void * airye_data[] = { (void *)cairy_wrap_e, (void *)cairy_wrap_e, };

static void * ellpj_data[] = { (void *)ellpj, (void *)ellpj,};

static void * expn_data[] = { (void *)expn, (void *)expn, };
static void * jn_data[] = { (void *)jn, (void *)jn, };
static void * kn_data[] = { (void *)kn, (void *)kn, };

static void * pdtrc_data[] = { (void *)pdtrc, (void *)pdtrc, };
static void * pdtr_data[] = { (void *)pdtr, (void *)pdtr, };
static void * pdtri_data[] = { (void *)pdtri, (void *)pdtri, };

static void * fresnl_data[] = { (void *)fresnl, (void *)fresnl };
static void * shichi_data[] = { (void *)shichi, (void *)shichi, };
static void * sici_data[] = { (void *)sici, (void *)sici, };

static void * stdtr_data[] = { (void *)stdtr, (void *)stdtr, };
static void * stdtri_data[] = { (void *)stdtri, (void *)stdtri, };

static void * yn_data[] = { (void *)yn, (void *)yn, };
static void * smirnov_data[] = { (void *)smirnov, (void *)smirnov, };
static void * smirnovi_data[] = { (void *)smirnovi, (void *)smirnovi, };

static void * bdtrc_data[] = { (void *)bdtrc, (void *)bdtrc, };
static void * bdtr_data[] = { (void *)bdtr, (void *)bdtr, };
static void * bdtri_data[] = { (void *)bdtri, (void *)bdtri, };
static void * btdtr_data[] = { (void *)btdtr, (void *)btdtr, };

static void * fdtrc_data[] = { (void *)fdtrc, (void *)fdtrc, };
static void * fdtr_data[] = { (void *)fdtr, (void *)fdtr, };
static void * fdtri_data[] = { (void *)fdtri, (void *)fdtri, };

static void * gdtrc_data[] = { (void *)gdtrc, (void *)gdtrc, };
static void * gdtr_data[] = { (void *)gdtr, (void *)gdtr, };

static void * hyp2f1_data[] = { (void *)hyp2f1, (void *)hyp2f1, };
static void * hyperg_data[] = { (void *)hyperg, (void *)hyperg, };
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
static void * Gamma_data[] = { (void *)Gamma, (void *)Gamma, };
static void * lgam_data[] = { (void *)lgam, (void *)lgam, };
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
static void * erf_data[] = { (void *)erf, (void *)erf, };
static void * ndtri_data[] = { (void *)ndtri, (void *)ndtri, };

static void * psi_data[] = { (void *)psi, (void *)psi, };
static void * rgamma_data[] = { (void *)rgamma, (void *)rgamma, };
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
static void * struve_data[] = { (void *)struve, (void *)struve, };

static void * zeta_data[] = { (void *)zeta, (void *)zeta, };
static void * zetac_data[] = { (void *)zetac, (void *)zetac, };

static void * kolmogorov_data[] = { (void *)kolmogorov, (void *)kolmogorov, };
static void * kolmogi_data[] = { (void *)kolmogi, (void *)kolmogi, };

static void * wofz_data[] = { (void *)cwofz_wrap, (void *)cwofz_wrap, };

static void * besselpoly_data[] = {(void *)besselpoly, (void *)besselpoly,};

static char cephes_6_types[] = { PyArray_FLOAT,  PyArray_FLOAT,  PyArray_FLOAT, PyArray_FLOAT, PyArray_FLOAT, PyArray_FLOAT, PyArray_DOUBLE,  PyArray_DOUBLE, PyArray_DOUBLE, PyArray_DOUBLE, PyArray_DOUBLE, PyArray_DOUBLE,};
static char cephes_5_types[] = { PyArray_FLOAT,  PyArray_FLOAT,  PyArray_FLOAT, PyArray_FLOAT, PyArray_FLOAT, PyArray_DOUBLE,  PyArray_DOUBLE, PyArray_DOUBLE, PyArray_DOUBLE, PyArray_DOUBLE,};
static char cephes_5c_types[] = { PyArray_FLOAT,  PyArray_FLOAT,  PyArray_FLOAT, PyArray_FLOAT, PyArray_FLOAT, PyArray_DOUBLE,  PyArray_DOUBLE, PyArray_DOUBLE, PyArray_DOUBLE, PyArray_DOUBLE, PyArray_CFLOAT, PyArray_CFLOAT, PyArray_CFLOAT, PyArray_CFLOAT, PyArray_CFLOAT, PyArray_CDOUBLE, PyArray_CDOUBLE, PyArray_CDOUBLE, PyArray_CDOUBLE, PyArray_CDOUBLE, };
static char cephes_4_types[] = { PyArray_FLOAT,  PyArray_FLOAT,  PyArray_FLOAT, PyArray_FLOAT, PyArray_DOUBLE,  PyArray_DOUBLE, PyArray_DOUBLE, PyArray_DOUBLE,};
static char cephes_3_types[] = { PyArray_FLOAT,  PyArray_FLOAT,  PyArray_FLOAT,   PyArray_DOUBLE,  PyArray_DOUBLE, PyArray_DOUBLE, };
static char cephes_3c_types[] = { PyArray_FLOAT, PyArray_FLOAT, PyArray_FLOAT, PyArray_DOUBLE, PyArray_DOUBLE, PyArray_DOUBLE, PyArray_FLOAT, PyArray_CFLOAT,  PyArray_CFLOAT, PyArray_DOUBLE, PyArray_CDOUBLE, PyArray_CDOUBLE, };
static char cephes_3cp_types[] = { PyArray_FLOAT, PyArray_CFLOAT,  PyArray_CFLOAT, PyArray_DOUBLE, PyArray_CDOUBLE, PyArray_CDOUBLE, };
static char cephes_2_types[] = { PyArray_FLOAT,  PyArray_FLOAT,  PyArray_DOUBLE,  PyArray_DOUBLE,  };
static char cephes_1c_types[] = { PyArray_CFLOAT, PyArray_CFLOAT, PyArray_CDOUBLE, PyArray_CDOUBLE, };

static void Cephes_InitOperators(PyObject *dictionary) {
	PyObject *f;

        cephes1_functions[0] = PyUFunc_f_f_As_d_d;
        cephes1_functions[1] = PyUFunc_d_d;
        cephes1c_functions[0] = PyUFunc_F_F_As_D_D;
	cephes1c_functions[1] = PyUFunc_D_D;
        cephes1_2_functions[0] = PyUFunc_f_ff_As_d_dd;
        cephes1_2_functions[1] = PyUFunc_d_dd;
        cephes1c_4_functions[0] = PyUFunc_f_ffff_As_d_dddd;
        cephes1c_4_functions[1] = PyUFunc_d_dddd;
        cephes1c_4_functions[2] = PyUFunc_F_FFFF_As_D_DDDD;
        cephes1c_4_functions[3] = PyUFunc_D_DDDD;
        cephes1cp_4_functions[0] = PyUFunc_f_ffff_As_D_DDDD;
        cephes1cp_4_functions[1] = PyUFunc_d_dddd_As_D_DDDD;
        cephes1cp_4_functions[2] = PyUFunc_F_FFFF_As_D_DDDD;
        cephes1cp_4_functions[3] = PyUFunc_D_DDDD;
        cephes2_functions[0] = PyUFunc_ff_f_As_dd_d;
        cephes2_functions[1] = PyUFunc_dd_d;
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
        cephes3a_functions[0] = PyUFunc_fff_f_As_iid_d;
        cephes3a_functions[1] = PyUFunc_ddd_d_As_iid_d;
        cephes4_functions[0] = PyUFunc_ffff_f_As_dddd_d;
        cephes4_functions[1] = PyUFunc_dddd_d;
        cephes4_2_functions[0] = PyUFunc_ffff_ff_As_dddd_dd;
        cephes4_2_functions[1] = PyUFunc_dddd_dd;
        cephes4a_2_functions[0] = PyUFunc_ffff_ff_As_dddi_dd;
        cephes4a_2_functions[1] = PyUFunc_dddd_dd_As_dddi_dd;
	
	/* Create function objects for each function call and insert
	   them in the dictionary */
	f = PyUFunc_FromFuncAndData(cephes3a_functions, bdtrc_data, cephes_4_types, 2, 3, 1, PyUFunc_None, "bdtrc", bdtrc_doc, 1);
	PyDict_SetItemString(dictionary, "bdtrc", f);
	Py_DECREF(f);
	f = PyUFunc_FromFuncAndData(cephes3a_functions, bdtr_data, cephes_4_types, 2, 3, 1, PyUFunc_None, "bdtr", bdtr_doc, 1);
	PyDict_SetItemString(dictionary, "bdtr", f);
	Py_DECREF(f);
	f = PyUFunc_FromFuncAndData(cephes3a_functions, bdtri_data, cephes_4_types, 2, 3, 1, PyUFunc_None, "bdtri", bdtri_doc, 1);
	PyDict_SetItemString(dictionary, "bdtri", f);
	Py_DECREF(f);
	f = PyUFunc_FromFuncAndData(cephes3_functions, btdtr_data, cephes_4_types, 2, 3, 1, PyUFunc_None, "btdtr", btdtr_doc, 1);
	PyDict_SetItemString(dictionary, "btdtr", f);
	Py_DECREF(f);

	f = PyUFunc_FromFuncAndData(cephes3a_functions, fdtrc_data, cephes_4_types, 2, 3, 1, PyUFunc_None, "fdtrc", fdtrc_doc, 1);
	PyDict_SetItemString(dictionary, "fdtrc", f);
	Py_DECREF(f);
	f = PyUFunc_FromFuncAndData(cephes3a_functions, fdtr_data, cephes_4_types, 2, 3, 1, PyUFunc_None, "fdtr", fdtr_doc, 1);
	PyDict_SetItemString(dictionary, "fdtr", f);
	Py_DECREF(f);
	f = PyUFunc_FromFuncAndData(cephes3a_functions, fdtri_data, cephes_4_types, 2, 3, 1, PyUFunc_None, "fdtri", fdtri_doc, 1);
	PyDict_SetItemString(dictionary, "fdtri", f);
	Py_DECREF(f);

	f = PyUFunc_FromFuncAndData(cephes3_functions, gdtrc_data, cephes_4_types, 2, 3, 1, PyUFunc_None, "gdtrc", gdtrc_doc, 1);
	PyDict_SetItemString(dictionary, "gdtrc", f);
	Py_DECREF(f);
	f = PyUFunc_FromFuncAndData(cephes3_functions, gdtr_data, cephes_4_types, 2, 3, 1, PyUFunc_None, "gdtr", gdtr_doc, 1);
	PyDict_SetItemString(dictionary, "gdtr", f);
	Py_DECREF(f);

	f = PyUFunc_FromFuncAndData(cephes4_functions, hyp2f1_data, cephes_5_types, 2, 4, 1, PyUFunc_None, "hyp2f1", hyp2f1_doc, 1);
	PyDict_SetItemString(dictionary, "hyp2f1", f);
	Py_DECREF(f);
	f = PyUFunc_FromFuncAndData(cephes3_functions, hyperg_data, cephes_4_types, 2, 3, 1, PyUFunc_None, "hyp1f1", hyp1f1_doc, 1);
	PyDict_SetItemString(dictionary, "hyp1f1", f);
	Py_DECREF(f);

	f = PyUFunc_FromFuncAndData(cephes4a_2_functions, hyp2f0_data, cephes_6_types, 2, 4, 2, PyUFunc_None, "hyp2f0", hyp2f0_doc, 1);
	PyDict_SetItemString(dictionary, "hyp2f0", f);
	Py_DECREF(f);
	f = PyUFunc_FromFuncAndData(cephes4_2_functions, onef2_data, cephes_6_types, 2, 4, 2, PyUFunc_None, "hyp1f2", hyp1f2_doc, 1);
	PyDict_SetItemString(dictionary, "hyp1f2", f);
	Py_DECREF(f);
	f = PyUFunc_FromFuncAndData(cephes4_2_functions, threef0_data, cephes_6_types, 2, 4, 2, PyUFunc_None, "hyp3f0", hyp3f0_doc, 1);
	PyDict_SetItemString(dictionary, "hyp3f0", f);
	Py_DECREF(f);

	f = PyUFunc_FromFuncAndData(cephes3_functions, incbet_data, cephes_4_types, 2, 3, 1, PyUFunc_None, "incbet", incbet_doc, 1);
	PyDict_SetItemString(dictionary, "incbet", f);
	Py_DECREF(f);
	f = PyUFunc_FromFuncAndData(cephes3_functions, incbi_data, cephes_4_types, 2, 3, 1, PyUFunc_None, "incbi", incbi_doc, 1);
	PyDict_SetItemString(dictionary, "incbi", f);
	Py_DECREF(f);


	f = PyUFunc_FromFuncAndData(cephes3a_functions, nbdtrc_data, cephes_4_types, 2, 3, 1, PyUFunc_None, "nbdtrc", nbdtrc_doc, 1);
	PyDict_SetItemString(dictionary, "nbdtrc", f);
	Py_DECREF(f);
	f = PyUFunc_FromFuncAndData(cephes3a_functions, nbdtr_data, cephes_4_types, 2, 3, 1, PyUFunc_None, "nbdtr", nbdtr_doc, 1);
	PyDict_SetItemString(dictionary, "nbdtr", f);
	Py_DECREF(f);
	f = PyUFunc_FromFuncAndData(cephes3a_functions, nbdtri_data, cephes_4_types, 2, 3, 1, PyUFunc_None, "nbdtri", nbdtri_doc, 1);
	PyDict_SetItemString(dictionary, "nbdtri", f);
	Py_DECREF(f);


	f = PyUFunc_FromFuncAndData(cephes2_functions, beta_data, cephes_3_types, 2, 2, 1, PyUFunc_None, "beta", beta_doc, 1);
	PyDict_SetItemString(dictionary, "beta", f);
	Py_DECREF(f);
	f = PyUFunc_FromFuncAndData(cephes2_functions, lbeta_data, cephes_3_types, 2, 2, 1, PyUFunc_None, "lbeta", lbeta_doc, 1);
	PyDict_SetItemString(dictionary, "lbeta", f);
	Py_DECREF(f);
	f = PyUFunc_FromFuncAndData(cephes1_functions, cbrt_data, cephes_2_types, 2, 1, 1, PyUFunc_None, "cbrt", cbrt_doc, 1);
	PyDict_SetItemString(dictionary, "cbrt", f);
	Py_DECREF(f);
	f = PyUFunc_FromFuncAndData(cephes2_functions, chdtrc_data, cephes_3_types, 2, 2, 1, PyUFunc_None, "chdtrc", chdtrc_doc, 1);
	PyDict_SetItemString(dictionary, "chdtrc", f);
	Py_DECREF(f);
	f = PyUFunc_FromFuncAndData(cephes2_functions, chdtr_data, cephes_3_types, 2, 2, 1, PyUFunc_None, "chdtr", chdtr_doc, 1);
	PyDict_SetItemString(dictionary, "chdtr", f);
	Py_DECREF(f);
	f = PyUFunc_FromFuncAndData(cephes2_functions, chdtri_data, cephes_3_types, 2, 2, 1, PyUFunc_None, "chdtri", chdtri_doc, 1);
	PyDict_SetItemString(dictionary, "chdtri", f);
	Py_DECREF(f);
	f = PyUFunc_FromFuncAndData(cephes1_functions, dawsn_data, cephes_2_types, 2, 1, 1, PyUFunc_None, "dawsn", dawsn_doc, 1);
	PyDict_SetItemString(dictionary, "dawsn", f);
	Py_DECREF(f);
	f = PyUFunc_FromFuncAndData(cephes2_functions, ellie_data, cephes_3_types, 2, 2, 1, PyUFunc_None, "ellie", ellie_doc, 1);
	PyDict_SetItemString(dictionary, "ellie", f);
	Py_DECREF(f);
	f = PyUFunc_FromFuncAndData(cephes2_functions, ellik_data, cephes_3_types, 2, 2, 1, PyUFunc_None, "ellik", ellik_doc, 1);
	PyDict_SetItemString(dictionary, "ellik", f);
	Py_DECREF(f);
	f = PyUFunc_FromFuncAndData(cephes1_functions, ellpe_data, cephes_2_types, 2, 1, 1, PyUFunc_None, "ellpe", ellpe_doc, 1);
	PyDict_SetItemString(dictionary, "ellpe", f);
	Py_DECREF(f);
	f = PyUFunc_FromFuncAndData(cephes1_functions, ellpk_data, cephes_2_types, 2, 1, 1, PyUFunc_None, "ellpk", ellpk_doc, 1);
	PyDict_SetItemString(dictionary, "ellpk", f);
	Py_DECREF(f);
	f = PyUFunc_FromFuncAndData(cephes1_functions, exp10_data, cephes_2_types, 2, 1, 1, PyUFunc_None, "exp10", exp10_doc, 1);
	PyDict_SetItemString(dictionary, "exp10", f);
	Py_DECREF(f);

	f = PyUFunc_FromFuncAndData(cephes1_functions, exp2_data, cephes_2_types, 2, 1, 1, PyUFunc_None, "exp2", exp2_doc, 1);
	PyDict_SetItemString(dictionary, "exp2", f);
	Py_DECREF(f);
	f = PyUFunc_FromFuncAndData(cephes1_functions, Gamma_data, cephes_2_types, 2, 1, 1, PyUFunc_None, "gamma", gamma_doc, 1);
	PyDict_SetItemString(dictionary, "gamma", f);
	Py_DECREF(f);
	f = PyUFunc_FromFuncAndData(cephes1_functions, lgam_data, cephes_2_types, 2, 1, 1, PyUFunc_None, "lgam", lgam_doc, 1);
	PyDict_SetItemString(dictionary, "lgam", f);
	Py_DECREF(f);
	f = PyUFunc_FromFuncAndData(cephes1_functions, i0_data, cephes_2_types, 2, 1, 1, PyUFunc_None, "i0", i0_doc, 1);
	PyDict_SetItemString(dictionary, "i0", f);
	Py_DECREF(f);
	f = PyUFunc_FromFuncAndData(cephes1_functions, i0e_data, cephes_2_types, 2, 1, 1, PyUFunc_None, "i0e", i0e_doc, 1);
	PyDict_SetItemString(dictionary, "i0e", f);
	Py_DECREF(f);
	f = PyUFunc_FromFuncAndData(cephes1_functions, i1_data, cephes_2_types, 2, 1, 1, PyUFunc_None, "i1", i1_doc, 1);
	PyDict_SetItemString(dictionary, "i1", f);
	Py_DECREF(f);
	f = PyUFunc_FromFuncAndData(cephes1_functions, i1e_data, cephes_2_types, 2, 1, 1, PyUFunc_None, "i1e", i1e_doc, 1);
	PyDict_SetItemString(dictionary, "i1e", f);
	Py_DECREF(f);

	f = PyUFunc_FromFuncAndData(cephes2_functions, igamc_data, cephes_3_types, 2, 2, 1, PyUFunc_None, "igamc", igamc_doc, 1);
	PyDict_SetItemString(dictionary, "igamc", f);
	Py_DECREF(f);
	f = PyUFunc_FromFuncAndData(cephes2_functions, igam_data, cephes_3_types, 2, 2, 1, PyUFunc_None, "igam", igam_doc, 1);
	PyDict_SetItemString(dictionary, "igam", f);
	Py_DECREF(f);
	f = PyUFunc_FromFuncAndData(cephes2_functions, igami_data, cephes_3_types, 2, 2, 1, PyUFunc_None, "igami", igami_doc, 1);
	PyDict_SetItemString(dictionary, "igami", f);
	Py_DECREF(f);
	f = PyUFunc_FromFuncAndData(cephes2c_functions, iv_data, cephes_3c_types, 4, 2, 1, PyUFunc_None, "iv", iv_doc, 1);
	PyDict_SetItemString(dictionary, "iv", f);
	Py_DECREF(f);
	f = PyUFunc_FromFuncAndData(cephes2cp_functions, ive_data, cephes_3c_types, 4, 2, 1, PyUFunc_None, "ive", ive_doc, 1);
	PyDict_SetItemString(dictionary, "ive", f);
	Py_DECREF(f);

	f = PyUFunc_FromFuncAndData(cephes2_4_functions, ellpj_data, cephes_6_types, 2, 2, 4, PyUFunc_None, "ellpj", ellpj_doc, 1);
	PyDict_SetItemString(dictionary, "ellpj", f);
	Py_DECREF(f);

	f = PyUFunc_FromFuncAndData(cephes2a_functions, expn_data, cephes_3_types, 2, 2, 1, PyUFunc_None, "expn", expn_doc, 1);
	PyDict_SetItemString(dictionary, "expn", f);
	Py_DECREF(f);
	f = PyUFunc_FromFuncAndData(cephes2a_functions, jn_data, cephes_3_types, 2, 2, 1, PyUFunc_None, "jn", jn_doc, 1);
	PyDict_SetItemString(dictionary, "jn", f);
	Py_DECREF(f);
	f = PyUFunc_FromFuncAndData(cephes2a_functions, kn_data, cephes_3_types, 2, 2, 1, PyUFunc_None, "kn", kn_doc, 1);
	PyDict_SetItemString(dictionary, "kn", f);
	Py_DECREF(f);
	f = PyUFunc_FromFuncAndData(cephes2a_functions, pdtrc_data, cephes_3_types, 2, 2, 1, PyUFunc_None, "pdtrc", pdtrc_doc, 1);
	PyDict_SetItemString(dictionary, "pdtrc", f);
	Py_DECREF(f);
	f = PyUFunc_FromFuncAndData(cephes2a_functions, pdtr_data, cephes_3_types, 2, 2, 1, PyUFunc_None, "pdtr", pdtr_doc, 1);
	PyDict_SetItemString(dictionary, "pdtr", f);
	Py_DECREF(f);
	f = PyUFunc_FromFuncAndData(cephes2a_functions, pdtri_data, cephes_3_types, 2, 2, 1, PyUFunc_None, "pdtri", pdtri_doc, 1);
	PyDict_SetItemString(dictionary, "pdtri", f);
	Py_DECREF(f);
	f = PyUFunc_FromFuncAndData(cephes2a_functions, stdtr_data, cephes_3_types, 2, 2, 1, PyUFunc_None, "stdtr", stdtr_doc, 1);
	PyDict_SetItemString(dictionary, "stdtr", f);
	Py_DECREF(f);
	f = PyUFunc_FromFuncAndData(cephes2a_functions, stdtri_data, cephes_3_types, 2, 2, 1, PyUFunc_None, "stdtri", stdtri_doc, 1);
	PyDict_SetItemString(dictionary, "stdtri", f);
	Py_DECREF(f);
	f = PyUFunc_FromFuncAndData(cephes2a_functions, yn_data, cephes_3_types, 2, 2, 1, PyUFunc_None, "yn", yn_doc, 1);
	PyDict_SetItemString(dictionary, "yn", f);
	Py_DECREF(f);
	f = PyUFunc_FromFuncAndData(cephes2a_functions, smirnov_data, cephes_3_types, 2, 2, 1, PyUFunc_None, "smirnov", smirnov_doc, 1);
	PyDict_SetItemString(dictionary, "smirnov", f);
	Py_DECREF(f);
	f = PyUFunc_FromFuncAndData(cephes2a_functions, smirnovi_data, cephes_3_types, 2, 2, 1, PyUFunc_None, "smirnovi", smirnovi_doc, 1);
	PyDict_SetItemString(dictionary, "smirnovi", f);
	Py_DECREF(f);

	f = PyUFunc_FromFuncAndData(cephes1c_4_functions, airy_data, cephes_5c_types, 4, 1, 4, PyUFunc_None, "airy", airy_doc, 1);
	PyDict_SetItemString(dictionary, "airy", f);
	Py_DECREF(f);

	f = PyUFunc_FromFuncAndData(cephes1cp_4_functions, airye_data, cephes_5c_types, 4, 1, 4, PyUFunc_None, "airye", airye_doc, 1);
	PyDict_SetItemString(dictionary, "airye", f);
	Py_DECREF(f);

	f = PyUFunc_FromFuncAndData(cephes1_2_functions, fresnl_data, cephes_3_types, 2, 1, 2, PyUFunc_None, "fresnl", fresnl_doc, 1);
	PyDict_SetItemString(dictionary, "fresnl", f);
	Py_DECREF(f);
	f = PyUFunc_FromFuncAndData(cephes1_2_functions, shichi_data, cephes_3_types, 2, 1, 2, PyUFunc_None, "shichi", shichi_doc, 1);
	PyDict_SetItemString(dictionary, "shichi", f);
	Py_DECREF(f);
	f = PyUFunc_FromFuncAndData(cephes1_2_functions, sici_data, cephes_3_types, 2, 1, 2, PyUFunc_None, "sici", sici_doc, 1);
	PyDict_SetItemString(dictionary, "sici", f);
	Py_DECREF(f);


	f = PyUFunc_FromFuncAndData(cephes1_functions, j0_data, cephes_2_types, 2, 1, 1, PyUFunc_None, "j0", j0_doc, 1);
	PyDict_SetItemString(dictionary, "j0", f);
	Py_DECREF(f);
	f = PyUFunc_FromFuncAndData(cephes1_functions, y0_data, cephes_2_types, 2, 1, 1, PyUFunc_None, "y0", y0_doc, 1);
	PyDict_SetItemString(dictionary, "y0", f);
	Py_DECREF(f);
	f = PyUFunc_FromFuncAndData(cephes1_functions, j1_data, cephes_2_types, 2, 1, 1, PyUFunc_None, "j1", j1_doc, 1);
	PyDict_SetItemString(dictionary, "j1", f);
	Py_DECREF(f);
	f = PyUFunc_FromFuncAndData(cephes1_functions, y1_data, cephes_2_types, 2, 1, 1, PyUFunc_None, "y1", y1_doc, 1);
	PyDict_SetItemString(dictionary, "y1", f);
	Py_DECREF(f);

	f = PyUFunc_FromFuncAndData(cephes2c_functions, jv_data, cephes_3c_types, 4, 2, 1, PyUFunc_None, "jv", jv_doc, 1);
	PyDict_SetItemString(dictionary, "jv", f);
	Py_DECREF(f);
	f = PyUFunc_FromFuncAndData(cephes2cp_functions, jve_data, cephes_3c_types, 4, 2, 1, PyUFunc_None, "jve", jve_doc, 1);
	PyDict_SetItemString(dictionary, "jve", f);
	Py_DECREF(f);
	f = PyUFunc_FromFuncAndData(cephes2c_functions, yv_data, cephes_3c_types, 4, 2, 1, PyUFunc_None, "yv", yv_doc, 1);
	PyDict_SetItemString(dictionary, "yv", f);
	Py_DECREF(f);
	f = PyUFunc_FromFuncAndData(cephes2cp_functions, yve_data, cephes_3c_types, 4, 2, 1, PyUFunc_None, "yve", yve_doc, 1);
	PyDict_SetItemString(dictionary, "yve", f);
	Py_DECREF(f);


	f = PyUFunc_FromFuncAndData(cephes1_functions, k0_data, cephes_2_types, 2, 1, 1, PyUFunc_None, "k0", k0_doc, 1);
	PyDict_SetItemString(dictionary, "k0", f);
	Py_DECREF(f);
	f = PyUFunc_FromFuncAndData(cephes1_functions, k0e_data, cephes_2_types, 2, 1, 1, PyUFunc_None, "k0e", k0e_doc, 1);
	PyDict_SetItemString(dictionary, "k0e", f);
	Py_DECREF(f);
	f = PyUFunc_FromFuncAndData(cephes1_functions, k1_data, cephes_2_types, 2, 1, 1, PyUFunc_None, "k1", k1_doc, 1);
	PyDict_SetItemString(dictionary, "k1", f);
	Py_DECREF(f);
	f = PyUFunc_FromFuncAndData(cephes1_functions, k1e_data, cephes_2_types, 2, 1, 1, PyUFunc_None, "k1e", k1e_doc, 1);
	PyDict_SetItemString(dictionary, "k1e", f);
	Py_DECREF(f);
	f = PyUFunc_FromFuncAndData(cephes2cp_functions, kv_data, cephes_3c_types, 4, 2, 1, PyUFunc_None, "kv", kv_doc, 1);
	PyDict_SetItemString(dictionary, "kv", f);
	Py_DECREF(f);
	f = PyUFunc_FromFuncAndData(cephes2cp_functions, kve_data, cephes_3c_types, 4, 2, 1, PyUFunc_None, "kve", kve_doc, 1);
	PyDict_SetItemString(dictionary, "kve", f);
	Py_DECREF(f);

	f = PyUFunc_FromFuncAndData(cephes2cpp_functions, hankel1_data, cephes_3cp_types, 2, 2, 1, PyUFunc_None, "hankel1", hankel1_doc, 1);
	PyDict_SetItemString(dictionary, "hankel1", f);
	Py_DECREF(f);
	f = PyUFunc_FromFuncAndData(cephes2cpp_functions, hankel1e_data, cephes_3cp_types, 2, 2, 1, PyUFunc_None, "hankel1e", hankel1e_doc, 1);
	PyDict_SetItemString(dictionary, "hankel1e", f);
	Py_DECREF(f);
	f = PyUFunc_FromFuncAndData(cephes2cpp_functions, hankel2_data, cephes_3cp_types, 2, 2, 1, PyUFunc_None, "hankel2", hankel2_doc, 1);
	PyDict_SetItemString(dictionary, "hankel2", f);
	Py_DECREF(f);
	f = PyUFunc_FromFuncAndData(cephes2cpp_functions, hankel2e_data, cephes_3cp_types, 2, 2, 1, PyUFunc_None, "hankel2e", hankel2e_doc, 1);
	PyDict_SetItemString(dictionary, "hankel2e", f);
	Py_DECREF(f);


	f = PyUFunc_FromFuncAndData(cephes1_functions, ndtr_data, cephes_2_types, 2, 1, 1, PyUFunc_None, "ndtr", ndtr_doc, 1);
	PyDict_SetItemString(dictionary, "ndtr", f);
	Py_DECREF(f);

	f = PyUFunc_FromFuncAndData(cephes1_functions, erfc_data, cephes_2_types, 2, 1, 1, PyUFunc_None, "erfc", erfc_doc, 1);
	PyDict_SetItemString(dictionary, "erfc", f);
	Py_DECREF(f);
	f = PyUFunc_FromFuncAndData(cephes1_functions, erf_data, cephes_2_types, 2, 1, 1, PyUFunc_None, "erf", erf_doc, 1);
	PyDict_SetItemString(dictionary, "erf", f);
	Py_DECREF(f);
	f = PyUFunc_FromFuncAndData(cephes1_functions, ndtri_data, cephes_2_types, 2, 1, 1, PyUFunc_None, "ndtri", ndtri_doc, 1);
	PyDict_SetItemString(dictionary, "ndtri", f);
	Py_DECREF(f);
	f = PyUFunc_FromFuncAndData(cephes1_functions, psi_data, cephes_2_types, 2, 1, 1, PyUFunc_None, "psi", psi_doc, 1);
	PyDict_SetItemString(dictionary, "psi", f);
	Py_DECREF(f);
	f = PyUFunc_FromFuncAndData(cephes1_functions, rgamma_data, cephes_2_types, 2, 1, 1, PyUFunc_None, "rgamma", rgamma_doc, 1);
	PyDict_SetItemString(dictionary, "rgamma", f);
	Py_DECREF(f);
	f = PyUFunc_FromFuncAndData(cephes1_functions, round_data, cephes_2_types, 2, 1, 1, PyUFunc_None, "round", round_doc, 1);
	PyDict_SetItemString(dictionary, "round", f);
	Py_DECREF(f);

	f = PyUFunc_FromFuncAndData(cephes1_functions, sindg_data, cephes_2_types, 2, 1, 1, PyUFunc_None, "sindg", sindg_doc, 1);
	PyDict_SetItemString(dictionary, "sindg", f);
	Py_DECREF(f);
	f = PyUFunc_FromFuncAndData(cephes1_functions, cosdg_data, cephes_2_types, 2, 1, 1, PyUFunc_None, "cosdg", cosdg_doc, 1);
	PyDict_SetItemString(dictionary, "cosdg", f);
	Py_DECREF(f);
	f = PyUFunc_FromFuncAndData(cephes3_functions, radian_data, cephes_4_types, 2, 3, 1, PyUFunc_None, "radian", radian_doc, 1);
	PyDict_SetItemString(dictionary, "radian", f);
	Py_DECREF(f);
	f = PyUFunc_FromFuncAndData(cephes1_functions, tandg_data, cephes_2_types, 2, 1, 1, PyUFunc_None, "tandg", tandg_doc, 1);
	PyDict_SetItemString(dictionary, "tandg", f);
	Py_DECREF(f);
	f = PyUFunc_FromFuncAndData(cephes1_functions, cotdg_data, cephes_2_types, 2, 1, 1, PyUFunc_None, "cotdg", cotdg_doc, 1);
	PyDict_SetItemString(dictionary, "cotdg", f);
	Py_DECREF(f);
	f = PyUFunc_FromFuncAndData(cephes1_functions, log1p_data, cephes_2_types, 2, 1, 1, PyUFunc_None, "log1p", log1p_doc, 1);
	PyDict_SetItemString(dictionary, "log1p", f);
	Py_DECREF(f);
	f = PyUFunc_FromFuncAndData(cephes1_functions, expm1_data, cephes_2_types, 2, 1, 1, PyUFunc_None, "expm1", expm1_doc, 1);
	PyDict_SetItemString(dictionary, "expm1", f);
	Py_DECREF(f);
	f = PyUFunc_FromFuncAndData(cephes1_functions, cosm1_data, cephes_2_types, 2, 1, 1, PyUFunc_None, "cosm1", cosm1_doc, 1);
	PyDict_SetItemString(dictionary, "cosm1", f);
	Py_DECREF(f);

	f = PyUFunc_FromFuncAndData(cephes1_functions, spence_data, cephes_2_types, 2, 1, 1, PyUFunc_None, "spence", spence_doc, 1);
	PyDict_SetItemString(dictionary, "spence", f);
	Py_DECREF(f);
	f = PyUFunc_FromFuncAndData(cephes1_functions, zetac_data, cephes_2_types, 2, 1, 1, PyUFunc_None, "zetac", zetac_doc, 1);
	PyDict_SetItemString(dictionary, "zetac", f);
	Py_DECREF(f);

	f = PyUFunc_FromFuncAndData(cephes2_functions, struve_data, cephes_3_types, 2, 2, 1, PyUFunc_None, "struve", struve_doc, 1);
	PyDict_SetItemString(dictionary, "struve", f);
	Py_DECREF(f);
	f = PyUFunc_FromFuncAndData(cephes2_functions, zeta_data, cephes_3_types, 2, 2, 1, PyUFunc_None, "zeta", zeta_doc, 1);
	PyDict_SetItemString(dictionary, "zeta", f);
	Py_DECREF(f);

	f = PyUFunc_FromFuncAndData(cephes1_functions, kolmogorov_data, cephes_2_types, 2, 1, 1, PyUFunc_None, "kolmogorov", kolmogorov_doc, 1);
	PyDict_SetItemString(dictionary, "kolmogorov", f);
	Py_DECREF(f);

	f = PyUFunc_FromFuncAndData(cephes1_functions, kolmogi_data, cephes_2_types, 2, 1, 1, PyUFunc_None, "kolmogorovi", kolmogorovi_doc, 1);
	PyDict_SetItemString(dictionary, "kolmogorovi", f);
	Py_DECREF(f);

	f = PyUFunc_FromFuncAndData(cephes1c_functions, wofz_data, cephes_1c_types, 2, 1, 1, PyUFunc_None, "wofz", wofz_doc, 1); 
	PyDict_SetItemString(dictionary, "wofz", f);
	Py_DECREF(f);

	f = PyUFunc_FromFuncAndData(cephes3_functions, besselpoly_data, cephes_4_types, 2, 3, 1, PyUFunc_None, "besselpoly", besselpoly_doc, 1);
	PyDict_SetItemString(dictionary, "besselpoly", f);
	Py_DECREF(f);

	f = PyUFunc_FromFuncAndData(isinf_functions, isinf_data, isinf_signatures, 
				    4, 1, 1, PyUFunc_None, "isinf", 
				    "isinf(x) returns non-zero if x is infinity.", 0);
	PyDict_SetItemString(dictionary, "isinf", f);
	Py_DECREF(f);

	f = PyUFunc_FromFuncAndData(isfinite_functions, isfinite_data, isinf_signatures, 
				    4, 1, 1, PyUFunc_None, "isfinite", 
				    "isfinite(x) returns non-zero if x is not infinity or not a number.", 0);
	PyDict_SetItemString(dictionary, "isfinite", f);
	Py_DECREF(f);

	f = PyUFunc_FromFuncAndData(isnan_functions, isnan_data, isinf_signatures, 
				    4, 1, 1, PyUFunc_None, "isnan", 
				    "isnan(x) returns non-zero if x is not a number.", 0);
	PyDict_SetItemString(dictionary, "isnan", f);
	Py_DECREF(f);

}



/* Decrement the reference count of all objects in **arrays. */
static void cleanup_arrays(PyArrayObject **arrays, int number)
{
  int k;
  for (k=0; k < number; k++)
    Py_XDECREF((PyObject *)arrays[k]);
  return;
}

/* All rank-0 arrays are converted to rank-1 arrays */
/* The number of dimensions of each array with rank less than
    the rank of the array with the most dimensions is increased by 
    prepending with a dimenson length of one so that all arrays have
    the same rank. */
/* Dimensions are checked and unmatched dimensions triggers an error */
/* Strides for dimensions whose real length is one is set to zero but the dimension
   length is set to the maximum dimensions for the collection of inputs  */
static int setup_input_arrays(PyTupleObject *inputs, PyArrayObject **inputarrays, int nin)
{
  int i, k;
  int maxrank=1;
  int *maxdims;
  PyObject *inputobj;
  PyArrayObject *ain, *tmparray;

  /* Convert nested sequences to arrays or just increase reference count
     if already an array */
  for (i=0; i < nin; i++) {
    ain = NULL;
    inputobj = PyTuple_GET_ITEM(inputs,i);
    ain = (PyArrayObject *)PyArray_FromObject(inputobj,PyArray_ObjectType(inputobj,0),0,0);
    if (NULL == ain) {
      cleanup_arrays(inputarrays,i);
      return -1;
    }
    if (ain->nd > maxrank) maxrank = ain->nd;
    inputarrays[i] = ain;
  }

  maxdims = (int*)malloc(sizeof(int)*maxrank);
  if (NULL == maxdims) {
    PyErr_SetString(PyExc_MemoryError, "arraymap: can't allocate memory for input arrays");
    cleanup_arrays(inputarrays,nin);
    return -1;
  }

  /* Reshape all arrays so they have the same rank (pre-pend with length 1 dimensions) */
  /* We want to replace the header information without copying the data. 
     Keeping the reference count correct can be tricky.
     We want to make a new array object with a different header and decrease the 
     reference count of the old one without deallocating the data section */
  for (i=0; i < nin; i++) {
    ain = inputarrays[i];

    /* Initialize all dimensions to 1 */
    /* Change array shape */
    for (k=0; k < maxrank; k++) 
      maxdims[k] = 1; 
    for (k=maxrank-ain->nd; k< maxrank; k++) 
      maxdims[k] = ain->dimensions[k-maxrank+ain->nd];

    tmparray = (PyArrayObject *)PyArray_FromDimsAndData(maxrank,maxdims,ain->descr->type,ain->data);
    if (NULL == tmparray) {
      free(maxdims);
      cleanup_arrays(inputarrays,nin);
      return -1;
    }
    tmparray->base = (PyObject *)ain;  /* When tmparray is deallocated ain will be too */
    inputarrays[i] = tmparray;  /* tmparray is new array */
  }

  /* Find dimension length for the output arrays (maximum length for each
     dimension) */
  for (k=0; k < maxrank; k++) { 
    maxdims[k] = 1;
    for (i=0; i < nin; i++) 
      if (inputarrays[i]->dimensions[k] > maxdims[k])
	maxdims[k] = inputarrays[i]->dimensions[k];
  }

  /* Now set all lengths for input array dimensions to maxdims 
       and make strides equal to zero for arrays whose
       real length is 1 for a particular dimension
  */

  for (i=0; i<nin; i++) {
    ain = inputarrays[i];
    for (k=0; k< maxrank; k++) {
      if (1 == ain->dimensions[k]) {
	ain->strides[k] = 0;
	ain->dimensions[k] = maxdims[k];
      }
      else if (ain->dimensions[k] != maxdims[k]) {
	PyErr_SetString(PyExc_ValueError,"arraymap: Frames are not aligned (mismatched dimensions).");
	cleanup_arrays(inputarrays,nin);
	free(maxdims);
	return -1;
      }
    }
  }

  free(maxdims);
  return 0;

}

static int type_from_object(PyObject *obj)
{
  if (PyArray_Check(obj))
    return ((PyArrayObject *)obj)->descr->type_num;
  if (PyInt_Check(obj) || PyLong_Check(obj)) return PyArray_LONG;
  if (PyFloat_Check(obj)) return PyArray_DOUBLE;
  if (PyComplex_Check(obj)) return PyArray_CDOUBLE;
  PyErr_SetString(PyExc_ValueError, "arraymap: Invalid type for output array.");
  return -1;
}

static int type_from_char(char typechar)
{
  switch(typechar) {
  case 'c': return PyArray_CHAR;
  case 'b': return PyArray_UBYTE;
  case '1': return PyArray_SBYTE;
  case 's': return PyArray_SHORT;
  case 'i': return PyArray_INT;
  case 'l': return PyArray_LONG;
  case 'f': return PyArray_FLOAT;
  case 'd': return PyArray_DOUBLE;
  case 'F': return PyArray_CFLOAT;
  case 'D': return PyArray_CDOUBLE;
  default:
    PyErr_SetString(PyExc_ValueError, "arraymap: Invalid type for array");
    return -1;
  }
}



/* This sets up the output arrays by calling the function with arguments 
     the first element of each input arrays.  If otypes is NULL, the
     returned value type is used to establish the type of the output
     arrays, otherwise the characters in otypes determine the
     output types */
static int setup_output_arrays(PyObject *func, PyArrayObject **inarr, int nin, PyArrayObject ***outarr, char *otypes, int numtypes)
{
  PyObject *arglist, *result;
  PyObject *tmpobject;
  PyArrayObject *tmparr;
  int i, nout;
  int nd, *dimensions, type_num;

  nd = inarr[0]->nd;
  dimensions = inarr[0]->dimensions;

  if ((numtypes == 0) || (otypes == NULL)) { 
    /* Call function to get number of outputs */

    /* Build argument list */
    if ((arglist = PyTuple_New(nin)) == NULL) {
      return -1;
    }
    /* Construct input argument by creating a tuple with an element
     from each input array (cast to an appropriate Python Object) */
    for (i=0; i < nin; i++) {
      tmparr = inarr[i];
      /* Get first data point */
      tmpobject = tmparr->descr->getitem((void *)tmparr->data);
      if (NULL == tmpobject) {
	Py_DECREF(arglist);
	return -1;
      }
      PyTuple_SET_ITEM(arglist, i, tmpobject);  /* arg1 owns reference to tmpobj now */
    }    
    /* Call Python Function */
    if ((result=PyEval_CallObject(func, arglist))==NULL) {
      Py_DECREF(arglist);
      return -1;
    }

    Py_DECREF(arglist);

    /* If result is a tuple, create output_arrays according 
       to output.  */
    if (PyTuple_Check(result)) {
      nout = PyTuple_GET_SIZE(result);
      *outarr = (PyArrayObject **)malloc(nout*sizeof(PyArrayObject *));
      if (NULL == *outarr) {
	PyErr_SetString(PyExc_MemoryError, "arraymap: Cannot allocate memory for output arrays.");
	Py_DECREF(result);
	return -1;
      }
      /* Create nout output arrays */
      for (i=0; i < nout; i++) {
	/* Determine type */
	if ((type_num=type_from_object(PyTuple_GET_ITEM(result, i)))==-1) {
	  cleanup_arrays(*outarr,i);
	  Py_DECREF(result);
	  free(*outarr);
	  return -1;
	}
	/* Create output array */
	(*outarr)[i] = (PyArrayObject *)PyArray_FromDims(nd,dimensions,type_num);
	if (NULL == (*outarr)[i]) {
	  cleanup_arrays(*outarr,i);
	  Py_DECREF(result);
	  free(*outarr);
	  return -1;
	}
      }
    }
    else {           /* Only a single output result */
      nout = 1;
      *outarr = (PyArrayObject **)malloc(nout*sizeof(PyArrayObject *));
      if (NULL==*outarr) {
	PyErr_SetString(PyExc_MemoryError,"arraymap: Cannot allocate memory for output arrays.");
	Py_DECREF(result);
	return -1;
      }
      if ((type_num = type_from_object(result))==-1) {
	Py_DECREF(result);
	free(*outarr);
	return -1;
      }
      (*outarr)[0] = (PyArrayObject *)PyArray_FromDims(nd,dimensions,type_num);
      if (NULL == (*outarr)[0]) {
	Py_DECREF(result);
	free(*outarr);
	return -1;
      }
    }
    Py_DECREF(result);
  }

  else { /* Character output types entered */
    nout = numtypes;
    *outarr = (PyArrayObject **)malloc(nout*sizeof(PyArrayObject *));
    if (NULL==*outarr) {
      PyErr_SetString(PyExc_MemoryError,"arraymap: Cannot allocate memory for output arrays.");
      return -1;
    }
    /* Create Output arrays */
    for (i=0; i < nout; i++) {
      /* Get type */
      if ((type_num = type_from_char(otypes[i]))==-1) {
	cleanup_arrays(*outarr,i);
	free(*outarr);
	return -1;
      }
      /* Create array */
      (*outarr)[i] = (PyArrayObject *)PyArray_FromDims(nd,dimensions,type_num);
      if (NULL == (*outarr)[i]) {
	cleanup_arrays(*outarr,i);
	free(*outarr);
	return -1;
      }
    }     
  } 
  return nout;
}


/* Corresponding dimensions are assumed to match, check before calling. */
/* No rank-0 arrays (make them rank-1 arrays) */

/* This replicates the standard Ufunc broadcasting rule that if the
   dimension length is one, incrementing does not occur for that dimension.  

   This is currently done by setting the stride in that dimension to
   zero during input array setup.

   The purpose of this function is to perform a for loop over arbitrary
   discontiguous N-D arrays, call the Python function for each set of 
   corresponding elements and place the results in the output_array.
*/   
#define INCREMENT(ret_ind, nd, max_ind) \
{ \
  int k; \
  k = (nd) - 1; \
  if (++(ret_ind)[k] >= (max_ind)[k]) { \
    while (k >= 0 && ((ret_ind)[k] >= (max_ind)[k]-1)) \
      (ret_ind)[k--] = 0; \
    if (k >= 0) (ret_ind)[k]++; \
    else (ret_ind)[0] = (max_ind)[0]; \
  }  \
}

#define CALCINDEX(indx, nd_index, strides, ndim) \
{ \
  int i; \
 \
  indx = 0; \
  for (i=0; i < (ndim); i++)  \
    indx += nd_index[i]*strides[i]; \
} 

static int loop_over_arrays(PyObject *func, PyArrayObject **inarr, int nin, PyArrayObject **outarr, int nout)
{
  int i, loop_index;
  int *nd_index, indx_in, indx_out;
  PyArrayObject *in, *out, *tmparr;
  PyObject *result, *tmpobj, *arglist;

  in = inarr[0];     /* For any shape information needed */
  out = outarr[0];
  /* Allocate the N-D index initalized to zero. */
  nd_index = (int *)calloc(in->nd,sizeof(int));
  if (NULL == nd_index) {
    PyErr_SetString(PyExc_MemoryError,"arraymap: Cannot allocate memory for arrays.");
    return -1;
  }
  /* Build argument list */
  if ((arglist = PyTuple_New(nin)) == NULL) {
    free(nd_index);
    return -1;
  }

  loop_index = PyArray_Size((PyObject *)in);  /* Total number of Python function calls */

  while(loop_index--) { 
    /* Create input argument list with current element from the input
       arrays 
    */
    for (i=0; i < nin; i++) {
      tmparr = inarr[i];
      /* Find linear index into this input array */
      CALCINDEX(indx_in,nd_index,tmparr->strides,in->nd);
      /* Get object at this index */
      tmpobj = tmparr->descr->getitem((void *)(tmparr->data+indx_in));
      if (NULL == tmpobj) {
	Py_DECREF(arglist);
	free(nd_index);
	return -1;
      }
      /* This steals reference of tmpobj */
      PyTuple_SET_ITEM(arglist, i, tmpobj);  
    }
    /* Call Python Function for this set of inputs */
    if ((result=PyEval_CallObject(func, arglist))==NULL) {
      Py_DECREF(arglist);
      free(nd_index);
      return -1;
    } 

    /* Find index into (all) output arrays */
    CALCINDEX(indx_out,nd_index,out->strides,out->nd);

    /* Copy the results to the output arrays */
    if (1==nout) {
      if ((outarr[0]->descr->setitem(result,(outarr[0]->data+indx_out)))==-1) {
	free(nd_index);
	Py_DECREF(arglist);
	Py_DECREF(result);
	return -1;
      }
    }
    else if (PyTuple_Check(result)) {
      for (i=0; i<nout; i++) {
	if ((outarr[i]->descr->setitem(PyTuple_GET_ITEM(result,i),(outarr[i]->data+indx_out)))==-1) {
	  free(nd_index);
	  Py_DECREF(arglist);
	  Py_DECREF(result);
	  return -1;
	}
      }
    }
    else { 
      PyErr_SetString(PyExc_ValueError,"arraymap: Function output of incorrect type.");
      free(nd_index);
      Py_DECREF(arglist);
      Py_DECREF(result);
      return -1;
    }

    /* Increment the index counter */
    INCREMENT(nd_index,in->nd,in->dimensions);
    Py_DECREF(result);

  }
  Py_DECREF(arglist);
  free(nd_index);
  return 0;
} 

static PyObject *build_output(PyArrayObject **outarr,int nout)
{
  int i;
  PyObject *out;

  if (1==nout) return PyArray_Return(outarr[0]);
  if ((out=PyTuple_New(nout))==NULL) return NULL;
  for (i=0; i<nout; i++) PyTuple_SET_ITEM(out,i,(PyObject *)(outarr[i]));
  return out;
}

static char arraymap_doc[] = "c1,..,cn = arraymap(pyfunc,inputs{,outputtypes})\n\n  Loop over the elements of the inputs tuple, applying pyfunc to the set\n  formed from each element of inputs.  Place the output in arrays c1,...,cn.\n  This function can make any pyfunc with scalar inputs and scalar outputs\n  emulate a ufunc.\n";

static PyObject *map_PyFunc(PyObject *self, PyObject *args)
{
  PyObject *Pyfunc, *out;
  PyTupleObject *inputs;
  PyArrayObject **inputarrays, **outputarrays;
  char *otypes=NULL;
  int nin, nout, numtypes = 0;
  if (!PyArg_ParseTuple ( args, "OO!|s#;cephes.arraymap", &Pyfunc, &PyTuple_Type, (PyObject **)&inputs, &otypes, &numtypes )) return NULL;

  if (!PyCallable_Check(Pyfunc)) {
    PyErr_SetString(PyExc_TypeError, "arraymap: First argument is not a callable object.");
    return NULL;  
  }
    
  nin = PyTuple_GET_SIZE(inputs);
  inputarrays = calloc(nin,sizeof(PyArrayObject *));
  if (NULL == inputarrays) {
     PyErr_SetString(PyExc_MemoryError,"arraymap: Cannot allocate memory for input arrays.");
     return NULL;
  }
  if (setup_input_arrays(inputs,inputarrays,nin) == -1) {
    free(inputarrays);
    return NULL;
  }

  /* Construct output arrays */
  if (-1 == (nout=setup_output_arrays(Pyfunc,inputarrays,nin,&outputarrays,otypes,numtypes))) {
    cleanup_arrays(inputarrays,nin);
    free(inputarrays);
    return NULL;
  }

  /* Loop over the input arrays and place in output-arrays */
  if (-1 == loop_over_arrays(Pyfunc,inputarrays,nin,outputarrays,nout)) {
    cleanup_arrays(inputarrays,nin);
    free(inputarrays);
    cleanup_arrays(outputarrays,nout);
    free(outputarrays);
    return NULL;
  }

  cleanup_arrays(inputarrays,nin);
  free(inputarrays);
  if ((out = build_output(outputarrays,nout))==NULL) {
    cleanup_arrays(outputarrays,nout);
    free(outputarrays);
    return NULL;
  }
  free(outputarrays);
  return out;
}
  
static struct PyMethodDef methods[] = {
  {"arraymap", map_PyFunc, METH_VARARGS, arraymap_doc},
  {NULL,		NULL, 0}		/* sentinel */
};

void initcephes() {
  PyObject *m, *d, *s;
  
  /* Create the module and add the functions */
  m = Py_InitModule("cephes", methods); 

  /* Import the ufunc objects */
  import_array();
  import_ufunc();

  /* Add some symbolic constants to the module */
  d = PyModule_GetDict(m);

  s = PyString_FromString("1.3");
  PyDict_SetItemString(d, "__version__", s);
  Py_DECREF(s);

  /* Load the cephes operators into the array module's namespace */
  Cephes_InitOperators(d); 
  
  /* Check for errors */
  if (PyErr_Occurred())
    Py_FatalError("can't initialize module cephes");
}

