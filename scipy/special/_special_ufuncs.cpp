#include <cmath>

#include "special/airy.h"
#include "special/bessel.h"
#include "special/fresnel.h"
#include "special/gamma.h"
#include "special/kelvin.h"
#include "special/mathieu.h"
#include "special/par_cyl.h"
#include "special/specfun.h"
#include "special/sphd_wave.h"
#include "special/struve.h"
#include "special/trig.h"
#include "special/zeta.h"
#include "ufunc.h"

// This is the extension module for the NumPy ufuncs in SciPy's special module. To create such a ufunc, call
// "SpecFun_NewUFunc" with a braced list of kernel functions that will become the ufunc overloads. There are
// many examples in the code below. The documentation of each ufunc is kept in a companion file called
// _special_ufuncs_docs.cpp.
//
// If you are adding a ufunc, you will also need to add the appropriate entry to scipy/special/functions.json.
// This allows the build process to generate a corresponding entry for scipy.special.cython_special.

using namespace std;

// This is needed by sf_error, it is defined in the Cython "_ufuncs_extra_code_common.pxi" for "_generate_pyx.py".
// It exists to "call PyUFunc_getfperr in a context where PyUFunc_API array is initialized", but here we are
// already in such a context.
extern "C" int wrap_PyUFunc_getfperr() { return PyUFunc_getfperr(); }

// The following are based off NumPy's dtype type codes and functions like PyUFunc_dd_d
// And also the following modifiers
// p for pointer
// r for reference
// c for const
// v for volative

using func_f_f_t = float (*)(float);
using func_d_d_t = double (*)(double);
using func_F_F_t = complex<float> (*)(complex<float>);
using func_D_D_t = complex<double> (*)(complex<double>);

using func_ff_f_t = float (*)(float, float);
using func_dd_d_t = double (*)(double, double);
using func_FF_F_t = complex<float> (*)(complex<float>, complex<float>);
using func_DD_D_t = complex<double> (*)(complex<double>, complex<double>);

using func_fff_f_t = float (*)(float, float, float);
using func_ddd_d_t = double (*)(double, double, double);

using func_f_FpFp_t = void (*)(float, complex<float> *, complex<float> *);
using func_d_DpDp_t = void (*)(double, complex<double> *, complex<double> *);

using func_f_fpfp_t = void (*)(float, float *, float *);
using func_d_dpdp_t = void (*)(double, double *, double *);

using func_ff_fpfp_t = void (*)(float, float, float *, float *);
using func_dd_dpdp_t = void (*)(double, double, double *, double *);

using func_fff_fpfp_t = void (*)(float, float, float, float *, float *);
using func_ddd_dpdp_t = void (*)(double, double, double, double *, double *);

using func_f_fpfpfpfp_t = void (*)(float, float *, float *, float *, float *);
using func_d_dpdpdpdp_t = void (*)(double, double *, double *, double *, double *);

using func_f_FpFpFpFp_t = void (*)(float, complex<float> *, complex<float> *, complex<float> *, complex<float> *);
using func_d_DpDpDpDp_t = void (*)(double, complex<double> *, complex<double> *, complex<double> *, complex<double> *);

using func_ffff_fpfp_t = void (*)(float, float, float, float, float *, float *);
using func_dddd_dpdp_t = void (*)(double, double, double, double, double *, double *);

using func_fffff_fpfp_t = void (*)(float, float, float, float, float, float *, float *);
using func_ddddd_dpdp_t = void (*)(double, double, double, double, double, double *, double *);

extern const char *_cospi_doc;
extern const char *_sinpi_doc;
extern const char *bei_doc;
extern const char *beip_doc;
extern const char *ber_doc;
extern const char *berp_doc;
extern const char *exp1_doc;
extern const char *expi_doc;
extern const char *gammaln_doc;
extern const char *it2i0k0_doc;
extern const char *it2j0y0_doc;
extern const char *it2struve0_doc;
extern const char *itairy_doc;
extern const char *iti0k0_doc;
extern const char *itj0y0_doc;
extern const char *itmodstruve0_doc;
extern const char *itstruve0_doc;
extern const char *kei_doc;
extern const char *keip_doc;
extern const char *kelvin_doc;
extern const char *ker_doc;
extern const char *kerp_doc;
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

static PyModuleDef _special_ufuncs_def = {
    PyModuleDef_HEAD_INIT,
    .m_name = "_special_ufuncs",
    .m_size = -1,
};

PyMODINIT_FUNC PyInit__special_ufuncs() {
    if (!SpecFun_Initialize()) {
        return nullptr;
    }

    PyObject *_special_ufuncs = PyModule_Create(&_special_ufuncs_def);
    if (_special_ufuncs == nullptr) {
        return nullptr;
    }

    PyObject *_cospi =
        SpecFun_NewUFunc({static_cast<func_f_f_t>(special::cospi), static_cast<func_d_d_t>(special::cospi),
                          static_cast<func_F_F_t>(special::cospi), static_cast<func_D_D_t>(special::cospi)},
                         "_cospi", _cospi_doc);
    PyModule_AddObjectRef(_special_ufuncs, "_cospi", _cospi);

    PyObject *_sinpi =
        SpecFun_NewUFunc({static_cast<func_f_f_t>(special::sinpi), static_cast<func_d_d_t>(special::sinpi),
                          static_cast<func_F_F_t>(special::sinpi), static_cast<func_D_D_t>(special::sinpi)},
                         "_sinpi", _sinpi_doc);
    PyModule_AddObjectRef(_special_ufuncs, "_sinpi", _sinpi);

    PyObject *_zeta = SpecFun_NewUFunc(
        {static_cast<func_ff_f_t>(special::zeta), static_cast<func_dd_d_t>(special::zeta)}, "_zeta", _zeta_doc);
    PyModule_AddObjectRef(_special_ufuncs, "_zeta", _zeta);

    PyObject *bei = SpecFun_NewUFunc({static_cast<func_f_f_t>(special::bei), static_cast<func_d_d_t>(special::bei)},
                                     "bei", bei_doc);
    PyModule_AddObjectRef(_special_ufuncs, "bei", bei);

    PyObject *beip = SpecFun_NewUFunc({static_cast<func_f_f_t>(special::beip), static_cast<func_d_d_t>(special::beip)},
                                      "beip", beip_doc);
    PyModule_AddObjectRef(_special_ufuncs, "beip", beip);

    PyObject *ber = SpecFun_NewUFunc({static_cast<func_f_f_t>(special::ber), static_cast<func_d_d_t>(special::ber)},
                                     "ber", ber_doc);
    PyModule_AddObjectRef(_special_ufuncs, "ber", ber);

    PyObject *berp = SpecFun_NewUFunc({static_cast<func_f_f_t>(special::berp), static_cast<func_d_d_t>(special::berp)},
                                      "berp", berp_doc);
    PyModule_AddObjectRef(_special_ufuncs, "berp", berp);

    PyObject *exp1 = SpecFun_NewUFunc({static_cast<func_f_f_t>(special::exp1), static_cast<func_d_d_t>(special::exp1),
                                       static_cast<func_F_F_t>(special::exp1), static_cast<func_D_D_t>(special::exp1)},
                                      "exp1", exp1_doc);
    PyModule_AddObjectRef(_special_ufuncs, "exp1", exp1);

    PyObject *expi = SpecFun_NewUFunc({static_cast<func_f_f_t>(special::expi), static_cast<func_d_d_t>(special::expi),
                                       static_cast<func_F_F_t>(special::expi), static_cast<func_D_D_t>(special::expi)},
                                      "expi", expi_doc);
    PyModule_AddObjectRef(_special_ufuncs, "expi", expi);

    PyObject *gammaln = SpecFun_NewUFunc(
        {static_cast<func_f_f_t>(special::gammaln), static_cast<func_d_d_t>(special::gammaln)}, "gammaln", gammaln_doc);
    PyModule_AddObjectRef(_special_ufuncs, "gammaln", gammaln);

    PyObject *it2i0k0 =
        SpecFun_NewUFunc({static_cast<func_f_fpfp_t>(special::it2i0k0), static_cast<func_d_dpdp_t>(special::it2i0k0)},
                         2, "it2i0k0", it2i0k0_doc);
    PyModule_AddObjectRef(_special_ufuncs, "it2i0k0", it2i0k0);

    PyObject *it2j0y0 =
        SpecFun_NewUFunc({static_cast<func_f_fpfp_t>(special::it2j0y0), static_cast<func_d_dpdp_t>(special::it2j0y0)},
                         2, "it2j0y0", it2j0y0_doc);
    PyModule_AddObjectRef(_special_ufuncs, "it2j0y0", it2j0y0);

    PyObject *it2struve0 =
        SpecFun_NewUFunc({static_cast<func_f_f_t>(special::it2struve0), static_cast<func_d_d_t>(special::it2struve0)},
                         "it2struve0", it2struve0_doc);
    PyModule_AddObjectRef(_special_ufuncs, "it2struve0", it2struve0);

    PyObject *itairy = SpecFun_NewUFunc(
        {static_cast<func_f_fpfpfpfp_t>(special::itairy), static_cast<func_d_dpdpdpdp_t>(special::itairy)}, 4, "itairy",
        itairy_doc);
    PyModule_AddObjectRef(_special_ufuncs, "itairy", itairy);

    PyObject *iti0k0 =
        SpecFun_NewUFunc({static_cast<func_f_fpfp_t>(special::it1i0k0), static_cast<func_d_dpdp_t>(special::it1i0k0)},
                         2, "iti0k0", iti0k0_doc);
    PyModule_AddObjectRef(_special_ufuncs, "iti0k0", iti0k0);

    PyObject *itj0y0 =
        SpecFun_NewUFunc({static_cast<func_f_fpfp_t>(special::it1j0y0), static_cast<func_d_dpdp_t>(special::it1j0y0)},
                         2, "itj0y0", itj0y0_doc);
    PyModule_AddObjectRef(_special_ufuncs, "itj0y0", itj0y0);

    PyObject *itmodstruve0 = SpecFun_NewUFunc(
        {static_cast<func_f_f_t>(special::itmodstruve0), static_cast<func_d_d_t>(special::itmodstruve0)},
        "itmodstruve0", itmodstruve0_doc);
    PyModule_AddObjectRef(_special_ufuncs, "itmodstruve0", itmodstruve0);

    PyObject *itstruve0 =
        SpecFun_NewUFunc({static_cast<func_f_f_t>(special::itstruve0), static_cast<func_d_d_t>(special::itstruve0)},
                         "itstruve0", itstruve0_doc);
    PyModule_AddObjectRef(_special_ufuncs, "itstruve0", itstruve0);

    PyObject *kei = SpecFun_NewUFunc({static_cast<func_f_f_t>(special::kei), static_cast<func_d_d_t>(special::kei)},
                                     "kei", kei_doc);
    PyModule_AddObjectRef(_special_ufuncs, "kei", kei);

    PyObject *keip = SpecFun_NewUFunc({static_cast<func_f_f_t>(special::keip), static_cast<func_d_d_t>(special::keip)},
                                      "keip", keip_doc);
    PyModule_AddObjectRef(_special_ufuncs, "keip", keip);

    PyObject *kelvin = SpecFun_NewUFunc(
        {static_cast<func_f_FpFpFpFp_t>(special::kelvin), static_cast<func_d_DpDpDpDp_t>(special::kelvin)}, 4, "kelvin",
        kelvin_doc);
    PyModule_AddObjectRef(_special_ufuncs, "kelvin", kelvin);

    PyObject *ker = SpecFun_NewUFunc({static_cast<func_f_f_t>(special::ker), static_cast<func_d_d_t>(special::ker)},
                                     "ker", ker_doc);
    PyModule_AddObjectRef(_special_ufuncs, "ker", ker);

    PyObject *kerp = SpecFun_NewUFunc({static_cast<func_f_f_t>(special::kerp), static_cast<func_d_d_t>(special::kerp)},
                                      "kerp", kerp_doc);
    PyModule_AddObjectRef(_special_ufuncs, "kerp", kerp);

    PyObject *mathieu_a =
        SpecFun_NewUFunc({static_cast<func_ff_f_t>(special::cem_cva), static_cast<func_dd_d_t>(special::cem_cva)},
                         "mathieu_a", mathieu_a_doc);
    PyModule_AddObjectRef(_special_ufuncs, "mathieu_a", mathieu_a);

    PyObject *mathieu_b =
        SpecFun_NewUFunc({static_cast<func_ff_f_t>(special::sem_cva), static_cast<func_dd_d_t>(special::sem_cva)},
                         "mathieu_b", mathieu_b_doc);
    PyModule_AddObjectRef(_special_ufuncs, "mathieu_b", mathieu_b);

    PyObject *mathieu_cem =
        SpecFun_NewUFunc({static_cast<func_fff_fpfp_t>(special::cem), static_cast<func_ddd_dpdp_t>(special::cem)}, 2,
                         "mathieu_cem", mathieu_cem_doc);
    PyModule_AddObjectRef(_special_ufuncs, "mathieu_cem", mathieu_cem);

    PyObject *mathieu_modcem1 =
        SpecFun_NewUFunc({static_cast<func_fff_fpfp_t>(special::mcm1), static_cast<func_ddd_dpdp_t>(special::mcm1)}, 2,
                         "mathieu_modcem1", mathieu_modcem1_doc);
    PyModule_AddObjectRef(_special_ufuncs, "mathieu_modcem1", mathieu_modcem1);

    PyObject *mathieu_modcem2 =
        SpecFun_NewUFunc({static_cast<func_fff_fpfp_t>(special::mcm2), static_cast<func_ddd_dpdp_t>(special::mcm2)}, 2,
                         "mathieu_modcem2", mathieu_modcem2_doc);
    PyModule_AddObjectRef(_special_ufuncs, "mathieu_modcem2", mathieu_modcem2);

    PyObject *mathieu_modsem1 =
        SpecFun_NewUFunc({static_cast<func_fff_fpfp_t>(special::msm1), static_cast<func_ddd_dpdp_t>(special::msm1)}, 2,
                         "mathieu_modsem1", mathieu_modsem1_doc);
    PyModule_AddObjectRef(_special_ufuncs, "mathieu_modsem1", mathieu_modsem1);

    PyObject *mathieu_modsem2 =
        SpecFun_NewUFunc({static_cast<func_fff_fpfp_t>(special::msm2), static_cast<func_ddd_dpdp_t>(special::msm2)}, 2,
                         "mathieu_modsem2", mathieu_modsem2_doc);
    PyModule_AddObjectRef(_special_ufuncs, "mathieu_modsem2", mathieu_modsem2);

    PyObject *mathieu_sem =
        SpecFun_NewUFunc({static_cast<func_fff_fpfp_t>(special::sem), static_cast<func_ddd_dpdp_t>(special::sem)}, 2,
                         "mathieu_sem", mathieu_sem_doc);
    PyModule_AddObjectRef(_special_ufuncs, "mathieu_sem", mathieu_sem);

    PyObject *modfresnelm = SpecFun_NewUFunc({static_cast<func_f_FpFp_t>(special::modified_fresnel_minus),
                                              static_cast<func_d_DpDp_t>(special::modified_fresnel_minus)},
                                             2, "modfresnelm", modfresnelm_doc);
    PyModule_AddObjectRef(_special_ufuncs, "modfresnelm", modfresnelm);

    PyObject *modfresnelp = SpecFun_NewUFunc({static_cast<func_f_FpFp_t>(special::modified_fresnel_plus),
                                              static_cast<func_d_DpDp_t>(special::modified_fresnel_plus)},
                                             2, "modfresnelp", modfresnelp_doc);
    PyModule_AddObjectRef(_special_ufuncs, "modfresnelp", modfresnelp);

    PyObject *obl_ang1 = SpecFun_NewUFunc({static_cast<func_ffff_fpfp_t>(special::oblate_aswfa_nocv),
                                           static_cast<func_dddd_dpdp_t>(special::oblate_aswfa_nocv)},
                                          2, "obl_ang1", obl_ang1_doc);
    PyModule_AddObjectRef(_special_ufuncs, "obl_ang1", obl_ang1);

    PyObject *obl_ang1_cv = SpecFun_NewUFunc(
        {static_cast<func_fffff_fpfp_t>(special::oblate_aswfa), static_cast<func_ddddd_dpdp_t>(special::oblate_aswfa)},
        2, "obl_ang1_cv", obl_ang1_cv_doc);
    PyModule_AddObjectRef(_special_ufuncs, "obl_ang1_cv", obl_ang1_cv);

    PyObject *obl_cv = SpecFun_NewUFunc(
        {static_cast<func_fff_f_t>(special::oblate_segv), static_cast<func_ddd_d_t>(special::oblate_segv)}, "obl_cv",
        obl_cv_doc);
    PyModule_AddObjectRef(_special_ufuncs, "obl_cv", obl_cv);

    PyObject *obl_rad1 = SpecFun_NewUFunc({static_cast<func_ffff_fpfp_t>(special::oblate_radial1_nocv),
                                           static_cast<func_dddd_dpdp_t>(special::oblate_radial1_nocv)},
                                          2, "obl_rad1", obl_rad1_doc);
    PyModule_AddObjectRef(_special_ufuncs, "obl_rad1", obl_rad1);

    PyObject *obl_rad1_cv = SpecFun_NewUFunc({static_cast<func_fffff_fpfp_t>(special::oblate_radial1),
                                              static_cast<func_ddddd_dpdp_t>(special::oblate_radial1)},
                                             2, "obl_rad1_cv", obl_rad1_cv_doc);
    PyModule_AddObjectRef(_special_ufuncs, "obl_rad1_cv", obl_rad1_cv);

    PyObject *obl_rad2 = SpecFun_NewUFunc({static_cast<func_ffff_fpfp_t>(special::oblate_radial2_nocv),
                                           static_cast<func_dddd_dpdp_t>(special::oblate_radial2_nocv)},
                                          2, "obl_rad2", obl_rad2_doc);
    PyModule_AddObjectRef(_special_ufuncs, "obl_rad2", obl_rad2);

    PyObject *obl_rad2_cv = SpecFun_NewUFunc({static_cast<func_fffff_fpfp_t>(special::oblate_radial2),
                                              static_cast<func_ddddd_dpdp_t>(special::oblate_radial2)},
                                             2, "obl_rad2_cv", obl_rad2_cv_doc);
    PyModule_AddObjectRef(_special_ufuncs, "obl_rad2_cv", obl_rad2_cv);

    PyObject *pbdv = SpecFun_NewUFunc(
        {static_cast<func_ff_fpfp_t>(special::pbdv), static_cast<func_dd_dpdp_t>(special::pbdv)}, 2, "pbdv", pbdv_doc);
    PyModule_AddObjectRef(_special_ufuncs, "pbdv", pbdv);

    PyObject *pbvv = SpecFun_NewUFunc(
        {static_cast<func_ff_fpfp_t>(special::pbvv), static_cast<func_dd_dpdp_t>(special::pbvv)}, 2, "pbvv", pbvv_doc);
    PyModule_AddObjectRef(_special_ufuncs, "pbvv", pbvv);

    PyObject *pbwa = SpecFun_NewUFunc(
        {static_cast<func_ff_fpfp_t>(special::pbwa), static_cast<func_dd_dpdp_t>(special::pbwa)}, 2, "pbwa", pbwa_doc);
    PyModule_AddObjectRef(_special_ufuncs, "pbwa", pbwa);

    PyObject *pro_ang1 = SpecFun_NewUFunc({static_cast<func_ffff_fpfp_t>(special::prolate_aswfa_nocv),
                                           static_cast<func_dddd_dpdp_t>(special::prolate_aswfa_nocv)},
                                          2, "pro_ang1", pro_ang1_doc);
    PyModule_AddObjectRef(_special_ufuncs, "pro_ang1", pro_ang1);

    PyObject *pro_ang1_cv = SpecFun_NewUFunc({static_cast<func_fffff_fpfp_t>(special::prolate_aswfa),
                                              static_cast<func_ddddd_dpdp_t>(special::prolate_aswfa)},
                                             2, "pro_ang1_cv", pro_ang1_cv_doc);
    PyModule_AddObjectRef(_special_ufuncs, "pro_ang1_cv", pro_ang1_cv);

    PyObject *pro_cv = SpecFun_NewUFunc(
        {static_cast<func_fff_f_t>(special::prolate_segv), static_cast<func_ddd_d_t>(special::prolate_segv)}, "obl_cv",
        pro_cv_doc);
    PyModule_AddObjectRef(_special_ufuncs, "pro_cv", pro_cv);

    PyObject *pro_rad1 = SpecFun_NewUFunc({static_cast<func_ffff_fpfp_t>(special::prolate_radial1_nocv),
                                           static_cast<func_dddd_dpdp_t>(special::prolate_radial1_nocv)},
                                          2, "pro_rad1", pro_rad1_doc);
    PyModule_AddObjectRef(_special_ufuncs, "pro_rad1", pro_rad1);

    PyObject *pro_rad1_cv = SpecFun_NewUFunc({static_cast<func_fffff_fpfp_t>(special::prolate_radial1),
                                              static_cast<func_ddddd_dpdp_t>(special::prolate_radial1)},
                                             2, "pro_rad1_cv", pro_rad1_cv_doc);
    PyModule_AddObjectRef(_special_ufuncs, "pro_rad1_cv", pro_rad1_cv);

    PyObject *pro_rad2 = SpecFun_NewUFunc({static_cast<func_ffff_fpfp_t>(special::prolate_radial2_nocv),
                                           static_cast<func_dddd_dpdp_t>(special::prolate_radial2_nocv)},
                                          2, "pro_rad2", pro_rad2_doc);
    PyModule_AddObjectRef(_special_ufuncs, "pro_rad2", pro_rad2);

    PyObject *pro_rad2_cv = SpecFun_NewUFunc({static_cast<func_fffff_fpfp_t>(special::prolate_radial2),
                                              static_cast<func_ddddd_dpdp_t>(special::prolate_radial2)},
                                             2, "pro_rad2_cv", pro_rad2_cv_doc);
    PyModule_AddObjectRef(_special_ufuncs, "pro_rad2_cv", pro_rad2_cv);

    return _special_ufuncs;
}
