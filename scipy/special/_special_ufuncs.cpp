#include "ufunc.h"

#include <cmath>
#include <complex>

#include "sf_error.h"
#include "xsf_special.h"
#include "xsf/airy.h"
#include "xsf/amos.h"
#include "xsf/bessel.h"
#include "xsf/binom.h"
#include "xsf/digamma.h"
#include "xsf/expint.h"
#include "xsf/fresnel.h"
#include "xsf/gamma.h"
#include "xsf/hyp2f1.h"
#include "xsf/iv_ratio.h"
#include "xsf/kelvin.h"
#include "xsf/lambertw.h"
#include "xsf/legendre.h"
#include "xsf/log_exp.h"
#include "xsf/mathieu.h"
#include "xsf/par_cyl.h"
#include "xsf/specfun.h"
#include "xsf/sph_bessel.h"
#include "xsf/sph_harm.h"
#include "xsf/sphd_wave.h"
#include "xsf/struve.h"
#include "xsf/trig.h"
#include "xsf/wright_bessel.h"
#include "xsf/zeta.h"

// This is the extension module for the NumPy ufuncs in SciPy's special module. To create such a ufunc, call
// "SpecFun_NewUFunc" with a braced list of kernel functions that will become the ufunc overloads. There are
// many examples in the code below. The documentation of each ufunc is kept in a companion file called
// _special_ufuncs_docs.cpp.
//
// If you are adding a ufunc, you will also need to add the appropriate entry to scipy/special/functions.json.
// This allows the build process to generate a corresponding entry for scipy.special.cython_special.

using namespace std;

// The following are based off NumPy's dtype type codes and functions like PyUFunc_dd_d

using cfloat = complex<float>;
using cdouble = complex<double>;

// 1 input, 1 output
using func_f_f_t = float (*)(float);
using func_d_d_t = double (*)(double);
using func_F_F_t = cfloat (*)(cfloat);
using func_D_D_t = cdouble (*)(cdouble);

// 1 input, 2 outputs
using func_f_ff_t = void (*)(float, float &, float &);
using func_d_dd_t = void (*)(double, double &, double &);
using func_f_FF_t = void (*)(float, cfloat &, cfloat &);
using func_d_DD_t = void (*)(double, cdouble &, cdouble &);

// 1 input, 4 outputs
using func_f_ffff_t = void (*)(float, float &, float &, float &, float &);
using func_d_dddd_t = void (*)(double, double &, double &, double &, double &);
using func_f_FFFF_t = void (*)(float, cfloat &, cfloat &, cfloat &, cfloat &);
using func_d_DDDD_t = void (*)(double, cdouble &, cdouble &, cdouble &, cdouble &);
using func_F_FFFF_t = void (*)(cfloat, cfloat &, cfloat &, cfloat &, cfloat &);
using func_D_DDDD_t = void (*)(cdouble, cdouble &, cdouble &, cdouble &, cdouble &);

// 2 inputs, 1 output
using func_qf_f_t = float (*)(long long int, float);
using func_qd_d_t = double (*)(long long int, double);
using func_ff_f_t = float (*)(float, float);
using func_dd_d_t = double (*)(double, double);
using func_FF_F_t = cfloat (*)(cfloat, cfloat);
using func_DD_D_t = cdouble (*)(cdouble, cdouble);
using func_fF_F_t = cfloat (*)(float, cfloat);
using func_dD_D_t = cdouble (*)(double, cdouble);
using func_lf_f_t = float (*)(long int, float);
using func_ld_d_t = double (*)(long int, double);
using func_lF_F_t = cfloat (*)(long int, cfloat);
using func_lD_D_t = cdouble (*)(long int, cdouble);

// 2 inputs, 2 outputs
using func_qf_ff_t = void (*)(long long int, float, float &, float &);
using func_qd_dd_t = void (*)(long long int, double, double &, double &);

// 2 inputs, 3 outputs
using func_qf_fff_t = void (*)(long long int, float, float &, float &, float &);
using func_qd_ddd_t = void (*)(long long int, double, double &, double &, double &);

// 2 inputs, 2 outputs
using func_ff_ff_t = void (*)(float, float, float &, float &);
using func_dd_dd_t = void (*)(double, double, double &, double &);
using func_lf_ff_t = void (*)(long int, float, float &, float &);
using func_ld_dd_t = void (*)(long int, double, double &, double &);

// 2 inputs, 3 outputs
using func_lf_fff_t = void (*)(long int, float, float &, float &, float &);
using func_ld_ddd_t = void (*)(long int, double, double &, double &, double &);

// 3 inputs, 1 output
using func_fff_f_t = float (*)(float, float, float);
using func_ddd_d_t = double (*)(double, double, double);
using func_Flf_F_t = cfloat (*)(cfloat, long int, float);
using func_Dld_D_t = cdouble (*)(cdouble, long int, double);

// 3 inputs, 2 outputs
using func_fff_ff_t = void (*)(float, float, float, float &, float &);
using func_ddd_dd_t = void (*)(double, double, double, double &, double &);

// 3 inputs, 1 output
using func_qqf_f_t = float (*)(long long int, long long int, float);
using func_qqd_d_t = double (*)(long long int, long long int, double);

// 3 inputs, 2 outputs
using func_qqf_ff_t = void (*)(long long int, long long int, float, float &, float &);
using func_qqd_dd_t = void (*)(long long int, long long int, double, double &, double &);

// 3 inputs, 3 outputs
using func_qqf_fff_t = void (*)(long long int, long long int, float, float &, float &, float &);
using func_qqd_ddd_t = void (*)(long long int, long long int, double, double &, double &, double &);

// 4 inputs, 1 outputs
using func_qqqF_F_t = cfloat (*)(long long int, long long int, long long int, cfloat);
using func_qqqD_D_t = cdouble (*)(long long int, long long int, long long int, cdouble);
using func_qqff_F_t = cfloat (*)(long long int, long long int, float, float);
using func_qqdd_D_t = cdouble (*)(long long int, long long int, double, double);
using func_ffff_f_t = float (*)(float, float, float, float);
using func_dddd_d_t = double (*)(double, double, double, double);
using func_fffF_F_t = cfloat (*)(float, float, float, cfloat);
using func_dddD_D_t = cdouble (*)(double, double, double, cdouble);
using func_ffff_F_t = cfloat (*)(float, float, float, float);
using func_dddd_D_t = cdouble (*)(double, double, double, double);

// 4 inputs, 2 outputs
using func_qqqf_ff_t = void (*)(long long int, long long int, long long int, float, float &, float &);
using func_ffff_ff_t = void (*)(float, float, float, float, float &, float &);
using func_qqqd_dd_t = void (*)(long long int, long long int, long long int, double, double &, double &);
using func_dddd_dd_t = void (*)(double, double, double, double, double &, double &);
using func_qqqF_FF_t = void (*)(long long int, long long int, long long int, cfloat, cfloat &, cfloat &);
using func_qqqD_DD_t = void (*)(long long int, long long int, long long int, cdouble, cdouble &, cdouble &);
using func_qqff_FF_t = void (*)(long long int, long long int, float, float, cfloat &, cfloat &);
using func_qqdd_DD_t = void (*)(long long int, long long int, double, double, cdouble &, cdouble &);
using func_qqff_FF2_t = void (*)(long long int, long long int, float, float, cfloat &, cfloat (&)[2]);
using func_qqdd_DD2_t = void (*)(long long int, long long int, double, double, cdouble &, cdouble (&)[2]);

// 4 inputs, 3 outputs
using func_qqqf_fff_t = void (*)(long long int, long long int, long long int, float, float &, float &, float &);
using func_qqqd_ddd_t = void (*)(long long int, long long int, long long int, double, double &, double &, double &);
using func_qqqF_FFF_t = void (*)(long long int, long long int, long long int, cfloat, cfloat &, cfloat &, cfloat &);
using func_qqqD_DDD_t = void (*)(long long int, long long int, long long int, cdouble, cdouble &, cdouble &, cdouble &);
using func_qqff_FFF_t = void (*)(long long int, long long int, float, float, cfloat &, cfloat &, cfloat &);
using func_qqdd_DDD_t = void (*)(long long int, long long int, double, double, cdouble &, cdouble &, cdouble &);
using func_qqff_FF2F22_t = void (*)(long long int, long long int, float, float, cfloat &, cfloat (&)[2],
                                    cfloat (&)[2][2]);
using func_qqdd_DD2D22_t = void (*)(long long int, long long int, double, double, cdouble &, cdouble (&)[2],
                                    cdouble (&)[2][2]);

// 5 inputs, 2 outputs
using func_fffff_ff_t = void (*)(float, float, float, float, float, float &, float &);
using func_ddddd_dd_t = void (*)(double, double, double, double, double, double &, double &);

#if (NPY_SIZEOF_LONGDOUBLE == NPY_SIZEOF_DOUBLE)
using func_g_g_t = double (*)(double);
using func_gg_g_t = double (*)(double);
#else
using func_g_g_t = long double (*)(long double);
using func_gg_g_t = long double (*)(long double);
#endif

extern const char *_cospi_doc;
extern const char *_sinpi_doc;
extern const char *airy_doc;
extern const char *airye_doc;
extern const char *bei_doc;
extern const char *beip_doc;
extern const char *ber_doc;
extern const char *berp_doc;
extern const char *binom_doc;
extern const char *exp1_doc;
extern const char *expi_doc;
extern const char *expit_doc;
extern const char *exprel_doc;
extern const char *gamma_doc;
extern const char *gammaln_doc;
extern const char *it2i0k0_doc;
extern const char *it2j0y0_doc;
extern const char *it2struve0_doc;
extern const char *itairy_doc;
extern const char *iti0k0_doc;
extern const char *itj0y0_doc;
extern const char *itmodstruve0_doc;
extern const char *itstruve0_doc;
extern const char *hankel1_doc;
extern const char *hankel1e_doc;
extern const char *hankel2_doc;
extern const char *hankel2e_doc;
extern const char *hyp2f1_doc;
extern const char *iv_doc;
extern const char *iv_ratio_doc;
extern const char *ive_doc;
extern const char *jv_doc;
extern const char *jve_doc;
extern const char *kei_doc;
extern const char *keip_doc;
extern const char *kelvin_doc;
extern const char *ker_doc;
extern const char *kerp_doc;
extern const char *kv_doc;
extern const char *kve_doc;
extern const char *lambertw_doc;
extern const char *logit_doc;
extern const char *loggamma_doc;
extern const char *log_expit_doc;
extern const char *log_wright_bessel_doc;
extern const char *mathieu_a_doc;
extern const char *mathieu_b_doc;
extern const char *mathieu_cem_doc;
extern const char *mathieu_modcem1_doc;
extern const char *mathieu_modcem2_doc;
extern const char *mathieu_modsem1_doc;
extern const char *mathieu_modsem2_doc;
extern const char *mathieu_sem_doc;
extern const char *modfresnelm_doc;
extern const char *modfresnelp_doc;
extern const char *obl_ang1_doc;
extern const char *obl_ang1_cv_doc;
extern const char *obl_cv_doc;
extern const char *obl_rad1_doc;
extern const char *obl_rad1_cv_doc;
extern const char *obl_rad2_doc;
extern const char *obl_rad2_cv_doc;
extern const char *_zeta_doc;
extern const char *pbdv_doc;
extern const char *pbvv_doc;
extern const char *pbwa_doc;
extern const char *pro_ang1_doc;
extern const char *pro_ang1_cv_doc;
extern const char *pro_cv_doc;
extern const char *pro_rad1_doc;
extern const char *pro_rad1_cv_doc;
extern const char *pro_rad2_doc;
extern const char *pro_rad2_cv_doc;
extern const char *psi_doc;
extern const char *rgamma_doc;
extern const char *scaled_exp1_doc;
extern const char *spherical_jn_doc;
extern const char *spherical_jn_d_doc;
extern const char *spherical_yn_doc;
extern const char *spherical_yn_d_doc;
extern const char *spherical_in_doc;
extern const char *spherical_in_d_doc;
extern const char *spherical_kn_doc;
extern const char *spherical_kn_d_doc;
extern const char *sph_harm_doc;
extern const char *wright_bessel_doc;
extern const char *yv_doc;
extern const char *yve_doc;

// This is needed by sf_error, it is defined in the Cython "_ufuncs_extra_code_common.pxi" for "_generate_pyx.py".
// It exists to "call PyUFunc_getfperr in a context where PyUFunc_API array is initialized", but here we are
// already in such a context.
extern "C" int wrap_PyUFunc_getfperr() { return PyUFunc_getfperr(); }

static PyModuleDef _special_ufuncs_def = {
    PyModuleDef_HEAD_INIT, "_special_ufuncs", NULL, -1, NULL, NULL, NULL, NULL, NULL};

PyMODINIT_FUNC PyInit__special_ufuncs() {
    import_array();
    import_umath();
    if (PyErr_Occurred()) {
        return NULL;
    }

    PyObject *_special_ufuncs = PyModule_Create(&_special_ufuncs_def);
    if (_special_ufuncs == nullptr) {
        return NULL;
    }

    PyObject *_cospi = SpecFun_NewUFunc({static_cast<func_f_f_t>(xsf::cospi), static_cast<func_d_d_t>(xsf::cospi),
                                         static_cast<func_F_F_t>(xsf::cospi), static_cast<func_D_D_t>(xsf::cospi)},
                                        "_cospi", _cospi_doc);
    PyModule_AddObjectRef(_special_ufuncs, "_cospi", _cospi);

    PyObject *_lambertw =
        SpecFun_NewUFunc({static_cast<func_Dld_D_t>(xsf::lambertw), static_cast<func_Flf_F_t>(xsf::lambertw)},
                         "_lambertw", lambertw_doc);
    PyModule_AddObjectRef(_special_ufuncs, "_lambertw", _lambertw);

    PyObject *_scaled_exp1 =
        SpecFun_NewUFunc({static_cast<func_d_d_t>(xsf::scaled_exp1), static_cast<func_f_f_t>(xsf::scaled_exp1)},
                         "_scaled_exp1", scaled_exp1_doc);
    PyModule_AddObjectRef(_special_ufuncs, "_scaled_exp1", _scaled_exp1);

    PyObject *_sinpi = SpecFun_NewUFunc({static_cast<func_f_f_t>(xsf::sinpi), static_cast<func_d_d_t>(xsf::sinpi),
                                         static_cast<func_F_F_t>(xsf::sinpi), static_cast<func_D_D_t>(xsf::sinpi)},
                                        "_sinpi", _sinpi_doc);
    PyModule_AddObjectRef(_special_ufuncs, "_sinpi", _sinpi);

    PyObject *_zeta = SpecFun_NewUFunc({static_cast<func_ff_f_t>(xsf::zeta), static_cast<func_dd_d_t>(xsf::zeta)},
                                       "_zeta", _zeta_doc);
    PyModule_AddObjectRef(_special_ufuncs, "_zeta", _zeta);

    PyObject *airy = SpecFun_NewUFunc({static_cast<func_f_ffff_t>(xsf::airy), static_cast<func_d_dddd_t>(xsf::airy),
                                       static_cast<func_F_FFFF_t>(xsf::airy), static_cast<func_D_DDDD_t>(xsf::airy)},
                                      4, "airy", airy_doc);
    PyModule_AddObjectRef(_special_ufuncs, "airy", airy);

    PyObject *airye = SpecFun_NewUFunc({static_cast<func_f_ffff_t>(xsf::airye), static_cast<func_d_dddd_t>(xsf::airye),
                                        static_cast<func_F_FFFF_t>(xsf::airye), static_cast<func_D_DDDD_t>(xsf::airye)},
                                       4, "airye", airye_doc);
    PyModule_AddObjectRef(_special_ufuncs, "airye", airye);

    PyObject *bei =
        SpecFun_NewUFunc({static_cast<func_f_f_t>(xsf::bei), static_cast<func_d_d_t>(xsf::bei)}, "bei", bei_doc);
    PyModule_AddObjectRef(_special_ufuncs, "bei", bei);

    PyObject *beip =
        SpecFun_NewUFunc({static_cast<func_f_f_t>(xsf::beip), static_cast<func_d_d_t>(xsf::beip)}, "beip", beip_doc);
    PyModule_AddObjectRef(_special_ufuncs, "beip", beip);

    PyObject *ber =
        SpecFun_NewUFunc({static_cast<func_f_f_t>(xsf::ber), static_cast<func_d_d_t>(xsf::ber)}, "ber", ber_doc);
    PyModule_AddObjectRef(_special_ufuncs, "ber", ber);

    PyObject *berp =
        SpecFun_NewUFunc({static_cast<func_f_f_t>(xsf::berp), static_cast<func_d_d_t>(xsf::berp)}, "berp", berp_doc);
    PyModule_AddObjectRef(_special_ufuncs, "berp", berp);

    PyObject *binom = SpecFun_NewUFunc({static_cast<func_ff_f_t>(xsf::binom), static_cast<func_dd_d_t>(xsf::binom)},
                                       "binom", binom_doc);
    PyModule_AddObjectRef(_special_ufuncs, "binom", binom);

    PyObject *exp1 = SpecFun_NewUFunc({static_cast<func_f_f_t>(xsf::exp1), static_cast<func_d_d_t>(xsf::exp1),
                                       static_cast<func_F_F_t>(xsf::exp1), static_cast<func_D_D_t>(xsf::exp1)},
                                      "exp1", exp1_doc);
    PyModule_AddObjectRef(_special_ufuncs, "exp1", exp1);

    PyObject *expi = SpecFun_NewUFunc({static_cast<func_f_f_t>(xsf::expi), static_cast<func_d_d_t>(xsf::expi),
                                       static_cast<func_F_F_t>(xsf::expi), static_cast<func_D_D_t>(xsf::expi)},
                                      "expi", expi_doc);
    PyModule_AddObjectRef(_special_ufuncs, "expi", expi);

    PyObject *expit = SpecFun_NewUFunc(
        {static_cast<func_d_d_t>(xsf::expit), static_cast<func_f_f_t>(xsf::expit), static_cast<func_g_g_t>(xsf::expit)},
        "expit", expit_doc);
    PyModule_AddObjectRef(_special_ufuncs, "expit", expit);

    PyObject *exprel = SpecFun_NewUFunc({static_cast<func_d_d_t>(xsf::exprel), static_cast<func_f_f_t>(xsf::exprel)},
                                        "exprel", exprel_doc);
    PyModule_AddObjectRef(_special_ufuncs, "exprel", exprel);

    PyObject *gamma = SpecFun_NewUFunc({static_cast<func_d_d_t>(xsf::gamma), static_cast<func_D_D_t>(xsf::gamma),
                                        static_cast<func_f_f_t>(xsf::gamma), static_cast<func_F_F_t>(xsf::gamma)},
                                       "gamma", gamma_doc);
    PyModule_AddObjectRef(_special_ufuncs, "gamma", gamma);

    PyObject *gammaln = SpecFun_NewUFunc({static_cast<func_f_f_t>(xsf::gammaln), static_cast<func_d_d_t>(xsf::gammaln)},
                                         "gammaln", gammaln_doc);
    PyModule_AddObjectRef(_special_ufuncs, "gammaln", gammaln);

    PyObject *hyp2f1 =
        SpecFun_NewUFunc({static_cast<func_dddd_d_t>(xsf::hyp2f1), static_cast<func_dddD_D_t>(xsf::hyp2f1),
                          static_cast<func_ffff_f_t>(xsf::hyp2f1), static_cast<func_fffF_F_t>(xsf::hyp2f1)},
                         "hyp2f1", hyp2f1_doc);
    PyModule_AddObjectRef(_special_ufuncs, "hyp2f1", hyp2f1);

    PyObject *hankel1 =
        SpecFun_NewUFunc({static_cast<func_fF_F_t>(xsf::cyl_hankel_1), static_cast<func_dD_D_t>(xsf::cyl_hankel_1)},
                         "hankel1", hankel1_doc);
    PyModule_AddObjectRef(_special_ufuncs, "hankel1", hankel1);

    PyObject *hankel1e =
        SpecFun_NewUFunc({static_cast<func_fF_F_t>(xsf::cyl_hankel_1e), static_cast<func_dD_D_t>(xsf::cyl_hankel_1e)},
                         "hankel1e", hankel1e_doc);
    PyModule_AddObjectRef(_special_ufuncs, "hankel1e", hankel1e);

    PyObject *hankel2 =
        SpecFun_NewUFunc({static_cast<func_fF_F_t>(xsf::cyl_hankel_2), static_cast<func_dD_D_t>(xsf::cyl_hankel_2)},
                         "hankel2", hankel2_doc);
    PyModule_AddObjectRef(_special_ufuncs, "hankel2", hankel2);

    PyObject *hankel2e =
        SpecFun_NewUFunc({static_cast<func_fF_F_t>(xsf::cyl_hankel_2e), static_cast<func_dD_D_t>(xsf::cyl_hankel_2e)},
                         "hankel2e", hankel2e_doc);
    PyModule_AddObjectRef(_special_ufuncs, "hankel2e", hankel2e);

    PyObject *it2i0k0 = SpecFun_NewUFunc(
        {static_cast<func_f_ff_t>(xsf::it2i0k0), static_cast<func_d_dd_t>(xsf::it2i0k0)}, 2, "it2i0k0", it2i0k0_doc);
    PyModule_AddObjectRef(_special_ufuncs, "it2i0k0", it2i0k0);

    PyObject *it2j0y0 = SpecFun_NewUFunc(
        {static_cast<func_f_ff_t>(xsf::it2j0y0), static_cast<func_d_dd_t>(xsf::it2j0y0)}, 2, "it2j0y0", it2j0y0_doc);
    PyModule_AddObjectRef(_special_ufuncs, "it2j0y0", it2j0y0);

    PyObject *it2struve0 =
        SpecFun_NewUFunc({static_cast<func_f_f_t>(xsf::it2struve0), static_cast<func_d_d_t>(xsf::it2struve0)},
                         "it2struve0", it2struve0_doc);
    PyModule_AddObjectRef(_special_ufuncs, "it2struve0", it2struve0);

    PyObject *itairy = SpecFun_NewUFunc(
        {static_cast<func_f_ffff_t>(xsf::itairy), static_cast<func_d_dddd_t>(xsf::itairy)}, 4, "itairy", itairy_doc);
    PyModule_AddObjectRef(_special_ufuncs, "itairy", itairy);

    PyObject *iti0k0 = SpecFun_NewUFunc(
        {static_cast<func_f_ff_t>(xsf::it1i0k0), static_cast<func_d_dd_t>(xsf::it1i0k0)}, 2, "iti0k0", iti0k0_doc);
    PyModule_AddObjectRef(_special_ufuncs, "iti0k0", iti0k0);

    PyObject *itj0y0 = SpecFun_NewUFunc(
        {static_cast<func_f_ff_t>(xsf::it1j0y0), static_cast<func_d_dd_t>(xsf::it1j0y0)}, 2, "itj0y0", itj0y0_doc);
    PyModule_AddObjectRef(_special_ufuncs, "itj0y0", itj0y0);

    PyObject *itmodstruve0 =
        SpecFun_NewUFunc({static_cast<func_f_f_t>(xsf::itmodstruve0), static_cast<func_d_d_t>(xsf::itmodstruve0)},
                         "itmodstruve0", itmodstruve0_doc);
    PyModule_AddObjectRef(_special_ufuncs, "itmodstruve0", itmodstruve0);

    PyObject *itstruve0 = SpecFun_NewUFunc(
        {static_cast<func_f_f_t>(xsf::itstruve0), static_cast<func_d_d_t>(xsf::itstruve0)}, "itstruve0", itstruve0_doc);
    PyModule_AddObjectRef(_special_ufuncs, "itstruve0", itstruve0);

    PyObject *iv =
        SpecFun_NewUFunc({static_cast<func_ff_f_t>(xsf::cyl_bessel_i), static_cast<func_dd_d_t>(xsf::cyl_bessel_i),
                          static_cast<func_fF_F_t>(xsf::cyl_bessel_i), static_cast<func_dD_D_t>(xsf::cyl_bessel_i)},
                         "iv", iv_doc);
    PyModule_AddObjectRef(_special_ufuncs, "iv", iv);

    PyObject *iv_ratio = SpecFun_NewUFunc(
        {static_cast<func_dd_d_t>(xsf::iv_ratio), static_cast<func_ff_f_t>(xsf::iv_ratio)}, "_iv_ratio", iv_ratio_doc);
    PyModule_AddObjectRef(_special_ufuncs, "_iv_ratio", iv_ratio);

    PyObject *ive =
        SpecFun_NewUFunc({static_cast<func_ff_f_t>(xsf::cyl_bessel_ie), static_cast<func_dd_d_t>(xsf::cyl_bessel_ie),
                          static_cast<func_fF_F_t>(xsf::cyl_bessel_ie), static_cast<func_dD_D_t>(xsf::cyl_bessel_ie)},
                         "ive", ive_doc);
    PyModule_AddObjectRef(_special_ufuncs, "ive", ive);

    PyObject *jv =
        SpecFun_NewUFunc({static_cast<func_ff_f_t>(xsf::cyl_bessel_j), static_cast<func_dd_d_t>(xsf::cyl_bessel_j),
                          static_cast<func_fF_F_t>(xsf::cyl_bessel_j), static_cast<func_dD_D_t>(xsf::cyl_bessel_j)},
                         "jv", jv_doc);
    PyModule_AddObjectRef(_special_ufuncs, "jv", jv);

    PyObject *jve =
        SpecFun_NewUFunc({static_cast<func_ff_f_t>(xsf::cyl_bessel_je), static_cast<func_dd_d_t>(xsf::cyl_bessel_je),
                          static_cast<func_fF_F_t>(xsf::cyl_bessel_je), static_cast<func_dD_D_t>(xsf::cyl_bessel_je)},
                         "jve", jve_doc);
    PyModule_AddObjectRef(_special_ufuncs, "jve", jve);

    PyObject *kei =
        SpecFun_NewUFunc({static_cast<func_f_f_t>(xsf::kei), static_cast<func_d_d_t>(xsf::kei)}, "kei", kei_doc);
    PyModule_AddObjectRef(_special_ufuncs, "kei", kei);

    PyObject *keip =
        SpecFun_NewUFunc({static_cast<func_f_f_t>(xsf::keip), static_cast<func_d_d_t>(xsf::keip)}, "keip", keip_doc);
    PyModule_AddObjectRef(_special_ufuncs, "keip", keip);

    PyObject *kelvin = SpecFun_NewUFunc(
        {static_cast<func_f_FFFF_t>(xsf::kelvin), static_cast<func_d_DDDD_t>(xsf::kelvin)}, 4, "kelvin", kelvin_doc);
    PyModule_AddObjectRef(_special_ufuncs, "kelvin", kelvin);

    PyObject *ker =
        SpecFun_NewUFunc({static_cast<func_f_f_t>(xsf::ker), static_cast<func_d_d_t>(xsf::ker)}, "ker", ker_doc);
    PyModule_AddObjectRef(_special_ufuncs, "ker", ker);

    PyObject *kerp =
        SpecFun_NewUFunc({static_cast<func_f_f_t>(xsf::kerp), static_cast<func_d_d_t>(xsf::kerp)}, "kerp", kerp_doc);
    PyModule_AddObjectRef(_special_ufuncs, "kerp", kerp);

    PyObject *kv =
        SpecFun_NewUFunc({static_cast<func_ff_f_t>(xsf::cyl_bessel_k), static_cast<func_dd_d_t>(xsf::cyl_bessel_k),
                          static_cast<func_fF_F_t>(xsf::cyl_bessel_k), static_cast<func_dD_D_t>(xsf::cyl_bessel_k)},
                         "kv", kv_doc);
    PyModule_AddObjectRef(_special_ufuncs, "kv", kv);

    PyObject *kve =
        SpecFun_NewUFunc({static_cast<func_ff_f_t>(xsf::cyl_bessel_ke), static_cast<func_dd_d_t>(xsf::cyl_bessel_ke),
                          static_cast<func_fF_F_t>(xsf::cyl_bessel_ke), static_cast<func_dD_D_t>(xsf::cyl_bessel_ke)},
                         "kve", kve_doc);
    PyModule_AddObjectRef(_special_ufuncs, "kve", kve);

    PyObject *log_expit =
        SpecFun_NewUFunc({static_cast<func_d_d_t>(xsf::log_expit), static_cast<func_f_f_t>(xsf::log_expit),
                          static_cast<func_g_g_t>(xsf::log_expit)},
                         "log_expit", log_expit_doc);
    PyModule_AddObjectRef(_special_ufuncs, "log_expit", log_expit);

    PyObject *log_wright_bessel = SpecFun_NewUFunc(
        {static_cast<func_ddd_d_t>(xsf::log_wright_bessel), static_cast<func_fff_f_t>(xsf::log_wright_bessel)},
        "log_wright_bessel", log_wright_bessel_doc);
    PyModule_AddObjectRef(_special_ufuncs, "log_wright_bessel", log_wright_bessel);

    PyObject *logit = SpecFun_NewUFunc(
        {static_cast<func_d_d_t>(xsf::logit), static_cast<func_f_f_t>(xsf::logit), static_cast<func_g_g_t>(xsf::logit)},
        "logit", logit_doc);
    PyModule_AddObjectRef(_special_ufuncs, "logit", logit);

    PyObject *loggamma =
        SpecFun_NewUFunc({static_cast<func_d_d_t>(xsf::loggamma), static_cast<func_D_D_t>(xsf::loggamma),
                          static_cast<func_f_f_t>(xsf::loggamma), static_cast<func_F_F_t>(xsf::loggamma)},
                         "loggamma", loggamma_doc);
    PyModule_AddObjectRef(_special_ufuncs, "loggamma", loggamma);

    PyObject *legendre_p = Py_BuildValue(
        "(N,N,N)",
        SpecFun_NewUFunc({static_cast<func_qd_d_t>(::legendre_p), static_cast<func_qf_f_t>(::legendre_p)}, "legendre_p",
                         nullptr),
        SpecFun_NewUFunc({static_cast<func_qd_dd_t>(::legendre_p), static_cast<func_qf_ff_t>(::legendre_p)}, 2,
                         "legendre_p", nullptr),
        SpecFun_NewUFunc({static_cast<func_qd_ddd_t>(::legendre_p), static_cast<func_qf_fff_t>(::legendre_p)}, 3,
                         "legendre_p", nullptr));
    PyModule_AddObjectRef(_special_ufuncs, "legendre_p", legendre_p);

    PyObject *assoc_legendre_p = Py_BuildValue(
        "{(O, i): N, (O, i): N, (O, i): N, (O, i): N, (O, i): N,(O, i): N}", Py_True, 0,
        SpecFun_NewUFunc({[](long long int n, long long int m, double z, long long int branch_cut) {
                              return ::assoc_legendre_p(assoc_legendre_norm, n, m, z, branch_cut);
                          },
                          [](long long int n, long long int m, float z, long long int branch_cut) {
                              return ::assoc_legendre_p(assoc_legendre_norm, n, m, z, branch_cut);
                          },
                          [](long long int n, long long int m, cdouble z, long long int branch_cut) {
                              return ::assoc_legendre_p(assoc_legendre_norm, n, m, z, branch_cut);
                          },
                          [](long long int n, long long int m, cfloat z, long long int branch_cut) {
                              return ::assoc_legendre_p(assoc_legendre_norm, n, m, z, branch_cut);
                          }},
                         "assoc_legendre_p", nullptr),
        Py_True, 1,
        SpecFun_NewUFunc(
            {[](long long int n, long long int m, double z, long long int branch_cut, double &res, double &res_jac) {
                 ::assoc_legendre_p(assoc_legendre_norm, n, m, z, branch_cut, res, res_jac);
             },
             [](long long int n, long long int m, float z, long long int branch_cut, float &res, float &res_jac) {
                 ::assoc_legendre_p(assoc_legendre_norm, n, m, z, branch_cut, res, res_jac);
             },
             [](long long int n, long long int m, cdouble z, long long int branch_cut, cdouble &res, cdouble &res_jac) {
                 ::assoc_legendre_p(assoc_legendre_norm, n, m, z, branch_cut, res, res_jac);
             },
             [](long long int n, long long int m, cfloat z, long long int branch_cut, cfloat &res, cfloat &res_jac) {
                 ::assoc_legendre_p(assoc_legendre_norm, n, m, z, branch_cut, res, res_jac);
             }},
            2, "assoc_legendre_p", nullptr),
        Py_True, 2,
        SpecFun_NewUFunc({[](long long int n, long long int m, double z, long long int branch_cut, double &res,
                             double &res_jac, double &res_hess) {
                              ::assoc_legendre_p(assoc_legendre_norm, n, m, z, branch_cut, res, res_jac, res_hess);
                          },
                          [](long long int n, long long int m, float z, long long int branch_cut, float &res,
                             float &res_jac, float &res_hess) {
                              ::assoc_legendre_p(assoc_legendre_norm, n, m, z, branch_cut, res, res_jac, res_hess);
                          },
                          [](long long int n, long long int m, cdouble z, long long int branch_cut, cdouble &res,
                             cdouble &res_jac, cdouble &res_hess) {
                              ::assoc_legendre_p(assoc_legendre_norm, n, m, z, branch_cut, res, res_jac, res_hess);
                          },
                          [](long long int n, long long int m, cfloat z, long long int branch_cut, cfloat &res,
                             cfloat &res_jac, cfloat &res_hess) {
                              ::assoc_legendre_p(assoc_legendre_norm, n, m, z, branch_cut, res, res_jac, res_hess);
                          }},
                         3, "assoc_legendre_p", nullptr),
        Py_False, 0,
        SpecFun_NewUFunc({[](long long int n, long long int m, double z, long long int branch_cut) {
                              return ::assoc_legendre_p(assoc_legendre_unnorm, n, m, z, branch_cut);
                          },
                          [](long long int n, long long int m, float z, long long int branch_cut) {
                              return ::assoc_legendre_p(assoc_legendre_unnorm, n, m, z, branch_cut);
                          },
                          [](long long int n, long long int m, cdouble z, long long int branch_cut) {
                              return ::assoc_legendre_p(assoc_legendre_unnorm, n, m, z, branch_cut);
                          },
                          [](long long int n, long long int m, cfloat z, long long int branch_cut) {
                              return ::assoc_legendre_p(assoc_legendre_unnorm, n, m, z, branch_cut);
                          }},
                         "assoc_legendre_p", nullptr),
        Py_False, 1,
        SpecFun_NewUFunc(
            {[](long long int n, long long int m, double z, long long int branch_cut, double &res, double &res_jac) {
                 ::assoc_legendre_p(assoc_legendre_unnorm, n, m, z, branch_cut, res, res_jac);
             },
             [](long long int n, long long int m, float z, long long int branch_cut, float &res, float &res_jac) {
                 ::assoc_legendre_p(assoc_legendre_unnorm, n, m, z, branch_cut, res, res_jac);
             },
             [](long long int n, long long int m, cdouble z, long long int branch_cut, cdouble &res, cdouble &res_jac) {
                 ::assoc_legendre_p(assoc_legendre_unnorm, n, m, z, branch_cut, res, res_jac);
             },
             [](long long int n, long long int m, cfloat z, long long int branch_cut, cfloat &res, cfloat &res_jac) {
                 ::assoc_legendre_p(assoc_legendre_unnorm, n, m, z, branch_cut, res, res_jac);
             }},
            2, "assoc_legendre_p", nullptr),
        Py_False, 2,
        SpecFun_NewUFunc({[](long long int n, long long int m, double z, long long int branch_cut, double &res,
                             double &res_jac, double &res_hess) {
                              ::assoc_legendre_p(assoc_legendre_unnorm, n, m, z, branch_cut, res, res_jac, res_hess);
                          },
                          [](long long int n, long long int m, float z, long long int branch_cut, float &res,
                             float &res_jac, float &res_hess) {
                              ::assoc_legendre_p(assoc_legendre_unnorm, n, m, z, branch_cut, res, res_jac, res_hess);
                          },
                          [](long long int n, long long int m, cdouble z, long long int branch_cut, cdouble &res,
                             cdouble &res_jac, cdouble &res_hess) {
                              ::assoc_legendre_p(assoc_legendre_unnorm, n, m, z, branch_cut, res, res_jac, res_hess);
                          },
                          [](long long int n, long long int m, cfloat z, long long int branch_cut, cfloat &res,
                             cfloat &res_jac, cfloat &res_hess) {
                              ::assoc_legendre_p(assoc_legendre_unnorm, n, m, z, branch_cut, res, res_jac, res_hess);
                          }},
                         3, "assoc_legendre_p", nullptr));
    PyModule_AddObjectRef(_special_ufuncs, "assoc_legendre_p", assoc_legendre_p);

    PyObject *mathieu_a = SpecFun_NewUFunc(
        {static_cast<func_ff_f_t>(xsf::cem_cva), static_cast<func_dd_d_t>(xsf::cem_cva)}, "mathieu_a", mathieu_a_doc);
    PyModule_AddObjectRef(_special_ufuncs, "mathieu_a", mathieu_a);

    PyObject *mathieu_b = SpecFun_NewUFunc(
        {static_cast<func_ff_f_t>(xsf::sem_cva), static_cast<func_dd_d_t>(xsf::sem_cva)}, "mathieu_b", mathieu_b_doc);
    PyModule_AddObjectRef(_special_ufuncs, "mathieu_b", mathieu_b);

    PyObject *mathieu_cem =
        SpecFun_NewUFunc({static_cast<func_fff_ff_t>(xsf::cem), static_cast<func_ddd_dd_t>(xsf::cem)}, 2, "mathieu_cem",
                         mathieu_cem_doc);
    PyModule_AddObjectRef(_special_ufuncs, "mathieu_cem", mathieu_cem);

    PyObject *mathieu_modcem1 =
        SpecFun_NewUFunc({static_cast<func_fff_ff_t>(xsf::mcm1), static_cast<func_ddd_dd_t>(xsf::mcm1)}, 2,
                         "mathieu_modcem1", mathieu_modcem1_doc);
    PyModule_AddObjectRef(_special_ufuncs, "mathieu_modcem1", mathieu_modcem1);

    PyObject *mathieu_modcem2 =
        SpecFun_NewUFunc({static_cast<func_fff_ff_t>(xsf::mcm2), static_cast<func_ddd_dd_t>(xsf::mcm2)}, 2,
                         "mathieu_modcem2", mathieu_modcem2_doc);
    PyModule_AddObjectRef(_special_ufuncs, "mathieu_modcem2", mathieu_modcem2);

    PyObject *mathieu_modsem1 =
        SpecFun_NewUFunc({static_cast<func_fff_ff_t>(xsf::msm1), static_cast<func_ddd_dd_t>(xsf::msm1)}, 2,
                         "mathieu_modsem1", mathieu_modsem1_doc);
    PyModule_AddObjectRef(_special_ufuncs, "mathieu_modsem1", mathieu_modsem1);

    PyObject *mathieu_modsem2 =
        SpecFun_NewUFunc({static_cast<func_fff_ff_t>(xsf::msm2), static_cast<func_ddd_dd_t>(xsf::msm2)}, 2,
                         "mathieu_modsem2", mathieu_modsem2_doc);
    PyModule_AddObjectRef(_special_ufuncs, "mathieu_modsem2", mathieu_modsem2);

    PyObject *mathieu_sem =
        SpecFun_NewUFunc({static_cast<func_fff_ff_t>(xsf::sem), static_cast<func_ddd_dd_t>(xsf::sem)}, 2, "mathieu_sem",
                         mathieu_sem_doc);
    PyModule_AddObjectRef(_special_ufuncs, "mathieu_sem", mathieu_sem);

    PyObject *modfresnelm = SpecFun_NewUFunc(
        {static_cast<func_f_FF_t>(xsf::modified_fresnel_minus), static_cast<func_d_DD_t>(xsf::modified_fresnel_minus)},
        2, "modfresnelm", modfresnelm_doc);
    PyModule_AddObjectRef(_special_ufuncs, "modfresnelm", modfresnelm);

    PyObject *modfresnelp = SpecFun_NewUFunc(
        {static_cast<func_f_FF_t>(xsf::modified_fresnel_plus), static_cast<func_d_DD_t>(xsf::modified_fresnel_plus)}, 2,
        "modfresnelp", modfresnelp_doc);
    PyModule_AddObjectRef(_special_ufuncs, "modfresnelp", modfresnelp);

    PyObject *obl_ang1 = SpecFun_NewUFunc(
        {static_cast<func_ffff_ff_t>(xsf::oblate_aswfa_nocv), static_cast<func_dddd_dd_t>(xsf::oblate_aswfa_nocv)}, 2,
        "obl_ang1", obl_ang1_doc);
    PyModule_AddObjectRef(_special_ufuncs, "obl_ang1", obl_ang1);

    PyObject *obl_ang1_cv = SpecFun_NewUFunc(
        {static_cast<func_fffff_ff_t>(xsf::oblate_aswfa), static_cast<func_ddddd_dd_t>(xsf::oblate_aswfa)}, 2,
        "obl_ang1_cv", obl_ang1_cv_doc);
    PyModule_AddObjectRef(_special_ufuncs, "obl_ang1_cv", obl_ang1_cv);

    PyObject *obl_cv =
        SpecFun_NewUFunc({static_cast<func_fff_f_t>(xsf::oblate_segv), static_cast<func_ddd_d_t>(xsf::oblate_segv)},
                         "obl_cv", obl_cv_doc);
    PyModule_AddObjectRef(_special_ufuncs, "obl_cv", obl_cv);

    PyObject *obl_rad1 = SpecFun_NewUFunc(
        {static_cast<func_ffff_ff_t>(xsf::oblate_radial1_nocv), static_cast<func_dddd_dd_t>(xsf::oblate_radial1_nocv)},
        2, "obl_rad1", obl_rad1_doc);
    PyModule_AddObjectRef(_special_ufuncs, "obl_rad1", obl_rad1);

    PyObject *obl_rad1_cv = SpecFun_NewUFunc(
        {static_cast<func_fffff_ff_t>(xsf::oblate_radial1), static_cast<func_ddddd_dd_t>(xsf::oblate_radial1)}, 2,
        "obl_rad1_cv", obl_rad1_cv_doc);
    PyModule_AddObjectRef(_special_ufuncs, "obl_rad1_cv", obl_rad1_cv);

    PyObject *obl_rad2 = SpecFun_NewUFunc(
        {static_cast<func_ffff_ff_t>(xsf::oblate_radial2_nocv), static_cast<func_dddd_dd_t>(xsf::oblate_radial2_nocv)},
        2, "obl_rad2", obl_rad2_doc);
    PyModule_AddObjectRef(_special_ufuncs, "obl_rad2", obl_rad2);

    PyObject *obl_rad2_cv = SpecFun_NewUFunc(
        {static_cast<func_fffff_ff_t>(xsf::oblate_radial2), static_cast<func_ddddd_dd_t>(xsf::oblate_radial2)}, 2,
        "obl_rad2_cv", obl_rad2_cv_doc);
    PyModule_AddObjectRef(_special_ufuncs, "obl_rad2_cv", obl_rad2_cv);

    PyObject *pbdv = SpecFun_NewUFunc({static_cast<func_ff_ff_t>(xsf::pbdv), static_cast<func_dd_dd_t>(xsf::pbdv)}, 2,
                                      "pbdv", pbdv_doc);
    PyModule_AddObjectRef(_special_ufuncs, "pbdv", pbdv);

    PyObject *pbvv = SpecFun_NewUFunc({static_cast<func_ff_ff_t>(xsf::pbvv), static_cast<func_dd_dd_t>(xsf::pbvv)}, 2,
                                      "pbvv", pbvv_doc);
    PyModule_AddObjectRef(_special_ufuncs, "pbvv", pbvv);

    PyObject *pbwa = SpecFun_NewUFunc({static_cast<func_ff_ff_t>(xsf::pbwa), static_cast<func_dd_dd_t>(xsf::pbwa)}, 2,
                                      "pbwa", pbwa_doc);
    PyModule_AddObjectRef(_special_ufuncs, "pbwa", pbwa);

    PyObject *pro_ang1 = SpecFun_NewUFunc(
        {static_cast<func_ffff_ff_t>(xsf::prolate_aswfa_nocv), static_cast<func_dddd_dd_t>(xsf::prolate_aswfa_nocv)}, 2,
        "pro_ang1", pro_ang1_doc);
    PyModule_AddObjectRef(_special_ufuncs, "pro_ang1", pro_ang1);

    PyObject *pro_ang1_cv = SpecFun_NewUFunc(
        {static_cast<func_fffff_ff_t>(xsf::prolate_aswfa), static_cast<func_ddddd_dd_t>(xsf::prolate_aswfa)}, 2,
        "pro_ang1_cv", pro_ang1_cv_doc);
    PyModule_AddObjectRef(_special_ufuncs, "pro_ang1_cv", pro_ang1_cv);

    PyObject *pro_cv =
        SpecFun_NewUFunc({static_cast<func_fff_f_t>(xsf::prolate_segv), static_cast<func_ddd_d_t>(xsf::prolate_segv)},
                         "obl_cv", pro_cv_doc);
    PyModule_AddObjectRef(_special_ufuncs, "pro_cv", pro_cv);

    PyObject *pro_rad1 = SpecFun_NewUFunc({static_cast<func_ffff_ff_t>(xsf::prolate_radial1_nocv),
                                           static_cast<func_dddd_dd_t>(xsf::prolate_radial1_nocv)},
                                          2, "pro_rad1", pro_rad1_doc);
    PyModule_AddObjectRef(_special_ufuncs, "pro_rad1", pro_rad1);

    PyObject *pro_rad1_cv = SpecFun_NewUFunc(
        {static_cast<func_fffff_ff_t>(xsf::prolate_radial1), static_cast<func_ddddd_dd_t>(xsf::prolate_radial1)}, 2,
        "pro_rad1_cv", pro_rad1_cv_doc);
    PyModule_AddObjectRef(_special_ufuncs, "pro_rad1_cv", pro_rad1_cv);

    PyObject *pro_rad2 = SpecFun_NewUFunc({static_cast<func_ffff_ff_t>(xsf::prolate_radial2_nocv),
                                           static_cast<func_dddd_dd_t>(xsf::prolate_radial2_nocv)},
                                          2, "pro_rad2", pro_rad2_doc);
    PyModule_AddObjectRef(_special_ufuncs, "pro_rad2", pro_rad2);

    PyObject *pro_rad2_cv = SpecFun_NewUFunc(
        {static_cast<func_fffff_ff_t>(xsf::prolate_radial2), static_cast<func_ddddd_dd_t>(xsf::prolate_radial2)}, 2,
        "pro_rad2_cv", pro_rad2_cv_doc);
    PyModule_AddObjectRef(_special_ufuncs, "pro_rad2_cv", pro_rad2_cv);

    PyObject *psi = SpecFun_NewUFunc({static_cast<func_d_d_t>(xsf::digamma), static_cast<func_D_D_t>(xsf::digamma),
                                      static_cast<func_f_f_t>(xsf::digamma), static_cast<func_F_F_t>(xsf::digamma)},
                                     "psi", psi_doc);
    PyModule_AddObjectRef(_special_ufuncs, "psi", psi);

    PyObject *rgamma = SpecFun_NewUFunc({static_cast<func_d_d_t>(xsf::rgamma), static_cast<func_D_D_t>(xsf::rgamma),
                                         static_cast<func_f_f_t>(xsf::rgamma), static_cast<func_F_F_t>(xsf::rgamma)},
                                        "rgamma", rgamma_doc);
    PyModule_AddObjectRef(_special_ufuncs, "rgamma", rgamma);

    PyObject *_spherical_jn =
        SpecFun_NewUFunc({static_cast<func_ld_d_t>(xsf::sph_bessel_j), static_cast<func_lD_D_t>(xsf::sph_bessel_j),
                          static_cast<func_lf_f_t>(xsf::sph_bessel_j), static_cast<func_lF_F_t>(xsf::sph_bessel_j)},
                         "_spherical_jn", spherical_jn_doc);
    PyModule_AddObjectRef(_special_ufuncs, "_spherical_jn", _spherical_jn);

    PyObject *_spherical_jn_d = SpecFun_NewUFunc(
        {static_cast<func_ld_d_t>(xsf::sph_bessel_j_jac), static_cast<func_lD_D_t>(xsf::sph_bessel_j_jac),
         static_cast<func_lf_f_t>(xsf::sph_bessel_j_jac), static_cast<func_lF_F_t>(xsf::sph_bessel_j_jac)},
        "_spherical_jn_d", spherical_jn_d_doc);
    PyModule_AddObjectRef(_special_ufuncs, "_spherical_jn_d", _spherical_jn_d);

    PyObject *_spherical_yn =
        SpecFun_NewUFunc({static_cast<func_ld_d_t>(xsf::sph_bessel_y), static_cast<func_lD_D_t>(xsf::sph_bessel_y),
                          static_cast<func_lf_f_t>(xsf::sph_bessel_y), static_cast<func_lF_F_t>(xsf::sph_bessel_y)},
                         "_spherical_yn", spherical_yn_doc);
    PyModule_AddObjectRef(_special_ufuncs, "_spherical_yn", _spherical_yn);

    PyObject *_spherical_yn_d = SpecFun_NewUFunc(
        {static_cast<func_ld_d_t>(xsf::sph_bessel_y_jac), static_cast<func_lD_D_t>(xsf::sph_bessel_y_jac),
         static_cast<func_lf_f_t>(xsf::sph_bessel_y_jac), static_cast<func_lF_F_t>(xsf::sph_bessel_y_jac)},
        "_spherical_yn_d", spherical_yn_d_doc);
    PyModule_AddObjectRef(_special_ufuncs, "_spherical_yn_d", _spherical_yn_d);

    PyObject *_spherical_in =
        SpecFun_NewUFunc({static_cast<func_ld_d_t>(xsf::sph_bessel_i), static_cast<func_lD_D_t>(xsf::sph_bessel_i),
                          static_cast<func_lf_f_t>(xsf::sph_bessel_i), static_cast<func_lF_F_t>(xsf::sph_bessel_i)},
                         "_spherical_in", spherical_in_doc);
    PyModule_AddObjectRef(_special_ufuncs, "_spherical_in", _spherical_in);

    PyObject *_spherical_in_d = SpecFun_NewUFunc(
        {static_cast<func_ld_d_t>(xsf::sph_bessel_i_jac), static_cast<func_lD_D_t>(xsf::sph_bessel_i_jac),
         static_cast<func_lf_f_t>(xsf::sph_bessel_i_jac), static_cast<func_lF_F_t>(xsf::sph_bessel_i_jac)},
        "_spherical_in_d", spherical_in_d_doc);
    PyModule_AddObjectRef(_special_ufuncs, "_spherical_in_d", _spherical_in_d);

    PyObject *_spherical_kn =
        SpecFun_NewUFunc({static_cast<func_ld_d_t>(xsf::sph_bessel_k), static_cast<func_lD_D_t>(xsf::sph_bessel_k),
                          static_cast<func_lf_f_t>(xsf::sph_bessel_k), static_cast<func_lF_F_t>(xsf::sph_bessel_k)},
                         "_spherical_kn", spherical_kn_doc);
    PyModule_AddObjectRef(_special_ufuncs, "_spherical_kn", _spherical_kn);

    PyObject *_spherical_kn_d = SpecFun_NewUFunc(
        {static_cast<func_ld_d_t>(xsf::sph_bessel_k_jac), static_cast<func_lD_D_t>(xsf::sph_bessel_k_jac),
         static_cast<func_lf_f_t>(xsf::sph_bessel_k_jac), static_cast<func_lF_F_t>(xsf::sph_bessel_k_jac)},
        "_spherical_kn_d", spherical_kn_d_doc);
    PyModule_AddObjectRef(_special_ufuncs, "_spherical_kn_d", _spherical_kn_d);

    PyObject *sph_legendre_p = Py_BuildValue(
        "(N,N,N)",
        SpecFun_NewUFunc({static_cast<func_qqd_d_t>(::sph_legendre_p), static_cast<func_qqf_f_t>(::sph_legendre_p)},
                         "sph_legendre_p", nullptr),
        SpecFun_NewUFunc({static_cast<func_qqd_dd_t>(::sph_legendre_p), static_cast<func_qqf_ff_t>(::sph_legendre_p)},
                         2, "sph_legendre_p", nullptr),
        SpecFun_NewUFunc({static_cast<func_qqd_ddd_t>(::sph_legendre_p), static_cast<func_qqf_fff_t>(::sph_legendre_p)},
                         3, "sph_legendre_p", nullptr));
    PyModule_AddObjectRef(_special_ufuncs, "sph_legendre_p", sph_legendre_p);

    PyObject *sph_harm =
        SpecFun_NewUFunc({static_cast<func_qqdd_D_t>(::sph_harm), static_cast<func_dddd_D_t>(::sph_harm),
                          static_cast<func_qqff_F_t>(::sph_harm), static_cast<func_ffff_F_t>(::sph_harm)},
                         "sph_harm", sph_harm_doc);
    PyModule_AddObjectRef(_special_ufuncs, "sph_harm", sph_harm);

    PyObject *sph_harm_y = Py_BuildValue(
        "(N,N,N)",
        SpecFun_NewUFunc({static_cast<func_qqdd_D_t>(::sph_harm_y), static_cast<func_qqff_F_t>(::sph_harm_y)},
                         "sph_harm_y", nullptr),
        SpecFun_NewGUFunc({static_cast<func_qqdd_DD2_t>(::sph_harm_y), static_cast<func_qqff_FF2_t>(::sph_harm_y)}, 2,
                          "sph_harm_y", nullptr, "(),(),(),()->(),(2)",
                          [](const npy_intp *dims, npy_intp *new_dims) { new_dims[0] = 2; }),
        SpecFun_NewGUFunc(
            {static_cast<func_qqdd_DD2D22_t>(::sph_harm_y), static_cast<func_qqff_FF2F22_t>(::sph_harm_y)}, 3,
            "sph_harm_y", nullptr, "(),(),(),()->(),(2),(2,2)", [](const npy_intp *dims, npy_intp *new_dims) {
                new_dims[0] = 2;

                new_dims[1] = 2;
                new_dims[2] = 2;
            }));
    PyModule_AddObjectRef(_special_ufuncs, "sph_harm_y", sph_harm_y);

    PyObject *wright_bessel =
        SpecFun_NewUFunc({static_cast<func_ddd_d_t>(xsf::wright_bessel), static_cast<func_fff_f_t>(xsf::wright_bessel)},
                         "wright_bessel", wright_bessel_doc);
    PyModule_AddObjectRef(_special_ufuncs, "wright_bessel", wright_bessel);

    PyObject *yv =
        SpecFun_NewUFunc({static_cast<func_ff_f_t>(xsf::cyl_bessel_y), static_cast<func_dd_d_t>(xsf::cyl_bessel_y),
                          static_cast<func_fF_F_t>(xsf::cyl_bessel_y), static_cast<func_dD_D_t>(xsf::cyl_bessel_y)},
                         "yv", yv_doc);
    PyModule_AddObjectRef(_special_ufuncs, "yv", yv);

    PyObject *yve =
        SpecFun_NewUFunc({static_cast<func_ff_f_t>(xsf::cyl_bessel_ye), static_cast<func_dd_d_t>(xsf::cyl_bessel_ye),
                          static_cast<func_fF_F_t>(xsf::cyl_bessel_ye), static_cast<func_dD_D_t>(xsf::cyl_bessel_ye)},
                         "yve", yve_doc);
    PyModule_AddObjectRef(_special_ufuncs, "yve", yve);

    return _special_ufuncs;
}
