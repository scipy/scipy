#include <xsf/numpy.h>

#include <cmath>
#include <complex>

#include "sf_error.h"
#include <xsf/airy.h>
#include <xsf/alg.h>
#include <xsf/amos.h>
#include <xsf/bessel.h>
#include <xsf/beta.h>
#include <xsf/binom.h>
#include <xsf/digamma.h>
#include <xsf/ellip.h>
#include <xsf/erf.h>
#include <xsf/exp.h>
#include <xsf/expint.h>
#include <xsf/fresnel.h>
#include <xsf/gamma.h>
#include <xsf/hyp2f1.h>
#include <xsf/iv_ratio.h>
#include <xsf/kelvin.h>
#include <xsf/lambertw.h>
#include <xsf/legendre.h>
#include <xsf/log.h>
#include <xsf/log_exp.h>
#include <xsf/mathieu.h>
#include <xsf/par_cyl.h>
#include <xsf/specfun.h>
#include <xsf/sph_bessel.h>
#include <xsf/sph_harm.h>
#include <xsf/sphd_wave.h>
#include <xsf/stats.h>
#include <xsf/struve.h>
#include <xsf/trig.h>
#include <xsf/wright_bessel.h>
#include <xsf/zeta.h>
#include "gen_harmonic.h"

// This is the extension module for the NumPy ufuncs in SciPy's special module. To create such a ufunc, call
// "xsf::numpy::ufunc" with a braced list of kernel functions that will become the ufunc overloads. There are
// many examples in the code below. The documentation of each ufunc is kept in a companion file called
// _special_ufuncs_docs.cpp.
//
// If you are adding a ufunc, you will also need to add the appropriate entry to scipy/special/functions.json.
// This allows the build process to generate a corresponding entry for scipy.special.cython_special.

extern const char *_cospi_doc;
extern const char *_sinpi_doc;
extern const char *_gen_harmonic_doc;
extern const char *_log1mexp_doc;
extern const char *_log1pmx_doc;
extern const char *_normalized_gen_harmonic_doc;
extern const char *airy_doc;
extern const char *airye_doc;
extern const char *bei_doc;
extern const char *beip_doc;
extern const char *ber_doc;
extern const char *berp_doc;
extern const char *besselpoly_doc;
extern const char *beta_doc;
extern const char *betaln_doc;
extern const char *binom_doc;
extern const char *cbrt_doc;
extern const char *cosdg_doc;
extern const char *cosm1_doc;
extern const char *cotdg_doc;
extern const char *dawsn_doc;
extern const char *ellipe_doc;
extern const char *ellipeinc_doc;
extern const char *ellipj_doc;
extern const char *ellipk_doc;
extern const char *ellipkm1_doc;
extern const char *ellipkinc_doc;
extern const char *erf_doc;
extern const char *erfc_doc;
extern const char *erfcx_doc;
extern const char *erfi_doc;
extern const char *exp1_doc;
extern const char *expm1_doc;
extern const char *exp2_doc;
extern const char *exp10_doc;
extern const char *expi_doc;
extern const char *expit_doc;
extern const char *exprel_doc;
extern const char *fresnel_doc;
extern const char *gamma_doc;
extern const char *gammainc_doc;
extern const char *gammaincinv_doc;
extern const char *gammaincc_doc;
extern const char *gammainccinv_doc;
extern const char *gammaln_doc;
extern const char *gammasgn_doc;
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
extern const char *i0_doc;
extern const char *i0e_doc;
extern const char *i1_doc;
extern const char *i1e_doc;
extern const char *iv_doc;
extern const char *iv_ratio_doc;
extern const char *iv_ratio_c_doc;
extern const char *ive_doc;
extern const char *j0_doc;
extern const char *j1_doc;
extern const char *jv_doc;
extern const char *jve_doc;
extern const char *kei_doc;
extern const char *keip_doc;
extern const char *kelvin_doc;
extern const char *ker_doc;
extern const char *kerp_doc;
extern const char *k0_doc;
extern const char *k0e_doc;
extern const char *k1_doc;
extern const char *k1e_doc;
extern const char *kv_doc;
extern const char *kve_doc;
extern const char *lambertw_doc;
extern const char *log1p_doc;
extern const char *logit_doc;
extern const char *loggamma_doc;
extern const char *log_expit_doc;
extern const char *log_ndtr_doc;
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
extern const char *ndtr_doc;
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
extern const char *radian_doc;
extern const char *rgamma_doc;
extern const char *_riemann_zeta_doc;
extern const char *scaled_exp1_doc;
extern const char *sindg_doc;
extern const char *spherical_jn_doc;
extern const char *spherical_jn_d_doc;
extern const char *spherical_yn_doc;
extern const char *spherical_yn_d_doc;
extern const char *spherical_in_doc;
extern const char *spherical_in_d_doc;
extern const char *spherical_kn_doc;
extern const char *spherical_kn_d_doc;
extern const char *struve_h_doc;
extern const char *struve_l_doc;
extern const char *tandg_doc;
extern const char *voigt_profile_doc;
extern const char *wofz_doc;
extern const char *wright_bessel_doc;
extern const char *xlogy_doc;
extern const char *xlog1py_doc;
extern const char *y0_doc;
extern const char *y1_doc;
extern const char *yv_doc;
extern const char *yve_doc;
extern const char *zetac_doc;


// Control error handling policy state
static PyObject* _set_action(PyObject* self, PyObject* args) {
    sf_error_t code;
    sf_action_t action;

    if (!PyArg_ParseTuple(args, "ii", &code, &action)) {
	return NULL;
    }

    sf_error_set_action(code, action);
    Py_RETURN_NONE;
}

static PyMethodDef _methods[] = {
    {"_set_action", _set_action, METH_VARARGS, NULL},
    {NULL, NULL, 0, NULL}
};


static int
_special_ufuncs_module_exec(PyObject *module)
{
    if (_import_array() < 0) { return -1; }
    if (_import_umath() < 0) { return -1; }

    PyObject *_cospi =
        xsf::numpy::ufunc({static_cast<xsf::numpy::f_f>(xsf::cospi), static_cast<xsf::numpy::d_d>(xsf::cospi),
                           static_cast<xsf::numpy::F_F>(xsf::cospi), static_cast<xsf::numpy::D_D>(xsf::cospi)},
                          "_cospi", _cospi_doc);
    PyModule_AddObjectRef(module, "_cospi", _cospi);

    PyObject *_gen_harmonic =
        xsf::numpy::ufunc({static_cast<xsf::numpy::ld_d>(gen_harmonic),
                           static_cast<xsf::numpy::qd_d>(gen_harmonic),
                           static_cast<xsf::numpy::dd_d>(gen_harmonic)},
                           "_gen_harmonic", _gen_harmonic_doc);
    PyModule_AddObjectRef(module, "_gen_harmonic", _gen_harmonic);

    // Define llld_d and qqqd_d here until they are added to xsf::numpy.
    using llld_d = double (*)(long int, long int, long int, double);
    using qqqd_d = double (*)(long long int, long long int, long long int, double);

    PyObject *_normalized_gen_harmonic =
        xsf::numpy::ufunc({static_cast<llld_d>(normalized_gen_harmonic),
                           static_cast<qqqd_d>(normalized_gen_harmonic),
                           static_cast<xsf::numpy::dddd_d>(normalized_gen_harmonic)},
                          "_normalized_gen_harmonic", _normalized_gen_harmonic_doc);
    PyModule_AddObjectRef(module, "_normalized_gen_harmonic", _normalized_gen_harmonic);

    PyObject *_lambertw = xsf::numpy::ufunc(
        {static_cast<xsf::numpy::Dld_D>(xsf::lambertw), static_cast<xsf::numpy::Flf_F>(xsf::lambertw)}, "_lambertw",
        lambertw_doc);
    PyModule_AddObjectRef(module, "_lambertw", _lambertw);

    PyObject *_scaled_exp1 = xsf::numpy::ufunc(
        {static_cast<xsf::numpy::d_d>(xsf::scaled_exp1), static_cast<xsf::numpy::f_f>(xsf::scaled_exp1)},
        "_scaled_exp1", scaled_exp1_doc);
    PyModule_AddObjectRef(module, "_scaled_exp1", _scaled_exp1);

    PyObject *_sinpi =
        xsf::numpy::ufunc({static_cast<xsf::numpy::f_f>(xsf::sinpi), static_cast<xsf::numpy::d_d>(xsf::sinpi),
                           static_cast<xsf::numpy::F_F>(xsf::sinpi), static_cast<xsf::numpy::D_D>(xsf::sinpi)},
                          "_sinpi", _sinpi_doc);
    PyModule_AddObjectRef(module, "_sinpi", _sinpi);

    PyObject *_zeta =
        xsf::numpy::ufunc({static_cast<xsf::numpy::ff_f>(xsf::zeta), static_cast<xsf::numpy::Ff_F>(xsf::zeta),
                           static_cast<xsf::numpy::dd_d>(xsf::zeta), static_cast<xsf::numpy::Dd_D>(xsf::zeta)},
                          "_zeta", _zeta_doc);
    PyModule_AddObjectRef(module, "_zeta", _zeta);

    PyObject *airy =
        xsf::numpy::ufunc({static_cast<xsf::numpy::f_ffff>(xsf::airy), static_cast<xsf::numpy::d_dddd>(xsf::airy),
                           static_cast<xsf::numpy::F_FFFF>(xsf::airy), static_cast<xsf::numpy::D_DDDD>(xsf::airy)},
                          4, "airy", airy_doc);
    PyModule_AddObjectRef(module, "airy", airy);

    PyObject *airye =
        xsf::numpy::ufunc({static_cast<xsf::numpy::f_ffff>(xsf::airye), static_cast<xsf::numpy::d_dddd>(xsf::airye),
                           static_cast<xsf::numpy::F_FFFF>(xsf::airye), static_cast<xsf::numpy::D_DDDD>(xsf::airye)},
                          4, "airye", airye_doc);
    PyModule_AddObjectRef(module, "airye", airye);

    PyObject *bei = xsf::numpy::ufunc({static_cast<xsf::numpy::f_f>(xsf::bei), static_cast<xsf::numpy::d_d>(xsf::bei)},
                                      "bei", bei_doc);
    PyModule_AddObjectRef(module, "bei", bei);

    PyObject *beip = xsf::numpy::ufunc(
        {static_cast<xsf::numpy::f_f>(xsf::beip), static_cast<xsf::numpy::d_d>(xsf::beip)}, "beip", beip_doc);
    PyModule_AddObjectRef(module, "beip", beip);

    PyObject *ber = xsf::numpy::ufunc({static_cast<xsf::numpy::f_f>(xsf::ber), static_cast<xsf::numpy::d_d>(xsf::ber)},
                                      "ber", ber_doc);
    PyModule_AddObjectRef(module, "ber", ber);

    PyObject *berp = xsf::numpy::ufunc(
        {static_cast<xsf::numpy::f_f>(xsf::berp), static_cast<xsf::numpy::d_d>(xsf::berp)}, "berp", berp_doc);
    PyModule_AddObjectRef(module, "berp", berp);

    PyObject *besselpoly = xsf::numpy::ufunc(
        {static_cast<xsf::numpy::ddd_d>(xsf::besselpoly), static_cast<xsf::numpy::fff_f>(xsf::besselpoly)},
        "besselpoly", besselpoly_doc);
    PyModule_AddObjectRef(module, "besselpoly", besselpoly);

    PyObject *beta = xsf::numpy::ufunc(
        {static_cast<xsf::numpy::dd_d>(xsf::beta), static_cast<xsf::numpy::ff_f>(xsf::beta)}, "beta", beta_doc);
    PyModule_AddObjectRef(module, "beta", beta);

    PyObject *betaln = xsf::numpy::ufunc(
        {static_cast<xsf::numpy::dd_d>(xsf::betaln), static_cast<xsf::numpy::ff_f>(xsf::betaln)}, "betaln", betaln_doc);
    PyModule_AddObjectRef(module, "betaln", betaln);

    PyObject *binom = xsf::numpy::ufunc(
        {static_cast<xsf::numpy::ff_f>(xsf::binom), static_cast<xsf::numpy::dd_d>(xsf::binom)}, "binom", binom_doc);
    PyModule_AddObjectRef(module, "binom", binom);

    PyObject *cbrt = xsf::numpy::ufunc(
        {static_cast<xsf::numpy::d_d>(xsf::cbrt), static_cast<xsf::numpy::f_f>(xsf::cbrt)}, "cbrt", cbrt_doc);
    PyModule_AddObjectRef(module, "cbrt", cbrt);

    PyObject *cosdg = xsf::numpy::ufunc(
        {static_cast<xsf::numpy::d_d>(xsf::cosdg), static_cast<xsf::numpy::f_f>(xsf::cosdg)}, "cosdg", cosdg_doc);
    PyModule_AddObjectRef(module, "cosdg", cosdg);

    PyObject *cosm1 = xsf::numpy::ufunc(
        {static_cast<xsf::numpy::d_d>(xsf::cosm1), static_cast<xsf::numpy::f_f>(xsf::cosm1)}, "cosm1", cosm1_doc);
    PyModule_AddObjectRef(module, "cosm1", cosm1);

    PyObject *cotdg = xsf::numpy::ufunc(
        {static_cast<xsf::numpy::d_d>(xsf::cotdg), static_cast<xsf::numpy::f_f>(xsf::cotdg)}, "cotdg", cotdg_doc);
    PyModule_AddObjectRef(module, "cotdg", cotdg);

    PyObject *ellipj = xsf::numpy::ufunc(
        {static_cast<xsf::numpy::dd_dddd>(xsf::ellipj), static_cast<xsf::numpy::ff_ffff>(xsf::ellipj)}, 4, "ellipj",
        ellipj_doc);
    PyModule_AddObjectRef(module, "ellipj", ellipj);

    PyObject *ellipe = xsf::numpy::ufunc(
        {static_cast<xsf::numpy::d_d>(xsf::ellipe), static_cast<xsf::numpy::f_f>(xsf::ellipe)}, "ellipe", ellipe_doc);
    PyModule_AddObjectRef(module, "ellipe", ellipe);

    PyObject *ellipeinc = xsf::numpy::ufunc(
        {static_cast<xsf::numpy::dd_d>(xsf::ellipeinc), static_cast<xsf::numpy::ff_f>(xsf::ellipeinc)}, "ellipeinc",
        ellipeinc_doc);
    PyModule_AddObjectRef(module, "ellipeinc", ellipeinc);

    PyObject *ellipk = xsf::numpy::ufunc(
        {static_cast<xsf::numpy::d_d>(xsf::ellipk), static_cast<xsf::numpy::f_f>(xsf::ellipk)}, "ellipk", ellipk_doc);
    PyModule_AddObjectRef(module, "ellipk", ellipk);

    PyObject *ellipkinc = xsf::numpy::ufunc(
        {static_cast<xsf::numpy::dd_d>(xsf::ellipkinc), static_cast<xsf::numpy::ff_f>(xsf::ellipkinc)}, "ellipkinc",
        ellipkinc_doc);
    PyModule_AddObjectRef(module, "ellipkinc", ellipkinc);

    PyObject *ellipkm1 =
        xsf::numpy::ufunc({static_cast<xsf::numpy::d_d>(xsf::ellipkm1), static_cast<xsf::numpy::f_f>(xsf::ellipkm1)},
                          "ellipkm1", ellipkm1_doc);
    PyModule_AddObjectRef(module, "ellipkm1", ellipkm1);

    PyObject *exp1 =
        xsf::numpy::ufunc({static_cast<xsf::numpy::f_f>(xsf::exp1), static_cast<xsf::numpy::d_d>(xsf::exp1),
                           static_cast<xsf::numpy::F_F>(xsf::exp1), static_cast<xsf::numpy::D_D>(xsf::exp1)},
                          "exp1", exp1_doc);
    PyModule_AddObjectRef(module, "exp1", exp1);

    PyObject *expi =
        xsf::numpy::ufunc({static_cast<xsf::numpy::f_f>(xsf::expi), static_cast<xsf::numpy::d_d>(xsf::expi),
                           static_cast<xsf::numpy::F_F>(xsf::expi), static_cast<xsf::numpy::D_D>(xsf::expi)},
                          "expi", expi_doc);
    PyModule_AddObjectRef(module, "expi", expi);

    PyObject *expit =
        xsf::numpy::ufunc({static_cast<xsf::numpy::d_d>(xsf::expit), static_cast<xsf::numpy::f_f>(xsf::expit),
                           static_cast<xsf::numpy::g_g>(xsf::expit)},
                          "expit", expit_doc);
    PyModule_AddObjectRef(module, "expit", expit);

    PyObject *exprel = xsf::numpy::ufunc(
        {static_cast<xsf::numpy::d_d>(xsf::exprel), static_cast<xsf::numpy::f_f>(xsf::exprel)}, "exprel", exprel_doc);
    PyModule_AddObjectRef(module, "exprel", exprel);

    PyObject *expm1 =
        xsf::numpy::ufunc({static_cast<xsf::numpy::d_d>(xsf::expm1), static_cast<xsf::numpy::f_f>(xsf::expm1),
                           static_cast<xsf::numpy::D_D>(xsf::expm1), static_cast<xsf::numpy::F_F>(xsf::expm1)},
                          "expm1", expm1_doc);
    PyModule_AddObjectRef(module, "expm1", expm1);

    PyObject *exp2 = xsf::numpy::ufunc(
        {static_cast<xsf::numpy::d_d>(xsf::exp2), static_cast<xsf::numpy::f_f>(xsf::exp2)}, "exp2", exp2_doc);
    PyModule_AddObjectRef(module, "exp2", exp2);

    PyObject *exp10 = xsf::numpy::ufunc(
        {static_cast<xsf::numpy::d_d>(xsf::exp10), static_cast<xsf::numpy::f_f>(xsf::exp10)}, "exp10", exp10_doc);
    PyModule_AddObjectRef(module, "exp10", exp10);

    PyObject *erf = xsf::numpy::ufunc({static_cast<xsf::numpy::d_d>(xsf::erf), static_cast<xsf::numpy::f_f>(xsf::erf),
                                       static_cast<xsf::numpy::D_D>(xsf::erf), static_cast<xsf::numpy::F_F>(xsf::erf)},
                                      "erf", erf_doc);
    PyModule_AddObjectRef(module, "erf", erf);

    PyObject *erfc =
        xsf::numpy::ufunc({static_cast<xsf::numpy::d_d>(xsf::erfc), static_cast<xsf::numpy::f_f>(xsf::erfc),
                           static_cast<xsf::numpy::D_D>(xsf::erfc), static_cast<xsf::numpy::F_F>(xsf::erfc)},
                          "erfc", erfc_doc);
    PyModule_AddObjectRef(module, "erfc", erfc);

    PyObject *erfcx =
        xsf::numpy::ufunc({static_cast<xsf::numpy::d_d>(xsf::erfcx), static_cast<xsf::numpy::f_f>(xsf::erfcx),
                           static_cast<xsf::numpy::D_D>(xsf::erfcx), static_cast<xsf::numpy::F_F>(xsf::erfcx)},
                          "erfcx", erfcx_doc);
    PyModule_AddObjectRef(module, "erfcx", erfcx);

    PyObject *erfi =
        xsf::numpy::ufunc({static_cast<xsf::numpy::d_d>(xsf::erfi), static_cast<xsf::numpy::f_f>(xsf::erfi),
                           static_cast<xsf::numpy::D_D>(xsf::erfi), static_cast<xsf::numpy::F_F>(xsf::erfi)},
                          "erfi", erfi_doc);
    PyModule_AddObjectRef(module, "erfi", erfi);

    PyObject *voigt_profile = xsf::numpy::ufunc(
        {static_cast<xsf::numpy::ddd_d>(xsf::voigt_profile), static_cast<xsf::numpy::fff_f>(xsf::voigt_profile)},
        "voigt_profile", voigt_profile_doc);
    PyModule_AddObjectRef(module, "voigt_profile", voigt_profile);

    PyObject *wofz = xsf::numpy::ufunc(
        {static_cast<xsf::numpy::D_D>(xsf::wofz), static_cast<xsf::numpy::F_F>(xsf::wofz)}, "wofz", wofz_doc);
    PyModule_AddObjectRef(module, "wofz", wofz);

    PyObject *dawsn =
        xsf::numpy::ufunc({static_cast<xsf::numpy::d_d>(xsf::dawsn), static_cast<xsf::numpy::f_f>(xsf::dawsn),
                           static_cast<xsf::numpy::D_D>(xsf::dawsn), static_cast<xsf::numpy::F_F>(xsf::dawsn)},
                          "dawsn", dawsn_doc);
    PyModule_AddObjectRef(module, "dawsn", dawsn);

    PyObject *ndtr =
        xsf::numpy::ufunc({static_cast<xsf::numpy::d_d>(xsf::ndtr), static_cast<xsf::numpy::f_f>(xsf::ndtr),
                           static_cast<xsf::numpy::D_D>(xsf::ndtr), static_cast<xsf::numpy::F_F>(xsf::ndtr)},
                          "ndtr", ndtr_doc);
    PyModule_AddObjectRef(module, "ndtr", ndtr);

    PyObject *log_ndtr =
        xsf::numpy::ufunc({static_cast<xsf::numpy::d_d>(xsf::log_ndtr), static_cast<xsf::numpy::f_f>(xsf::log_ndtr),
                           static_cast<xsf::numpy::D_D>(xsf::log_ndtr), static_cast<xsf::numpy::F_F>(xsf::log_ndtr)},
                          "log_ndtr", log_ndtr_doc);
    PyModule_AddObjectRef(module, "log_ndtr", log_ndtr);

    PyObject *fresnel =
        xsf::numpy::ufunc({static_cast<xsf::numpy::d_dd>(xsf::fresnel), static_cast<xsf::numpy::f_ff>(xsf::fresnel),
                           static_cast<xsf::numpy::D_DD>(xsf::fresnel), static_cast<xsf::numpy::F_FF>(xsf::fresnel)},
                          2, "fresnel", fresnel_doc);
    PyModule_AddObjectRef(module, "fresnel", fresnel);

    PyObject *gamma =
        xsf::numpy::ufunc({static_cast<xsf::numpy::d_d>(xsf::gamma), static_cast<xsf::numpy::D_D>(xsf::gamma),
                           static_cast<xsf::numpy::f_f>(xsf::gamma), static_cast<xsf::numpy::F_F>(xsf::gamma)},
                          "gamma", gamma_doc);
    PyModule_AddObjectRef(module, "gamma", gamma);

    PyObject *gammainc =
        xsf::numpy::ufunc({static_cast<xsf::numpy::dd_d>(xsf::gammainc), static_cast<xsf::numpy::ff_f>(xsf::gammainc)},
                          "gammainc", gammainc_doc);
    PyModule_AddObjectRef(module, "gammainc", gammainc);

    PyObject *gammaincinv = xsf::numpy::ufunc(
        {static_cast<xsf::numpy::dd_d>(xsf::gammaincinv), static_cast<xsf::numpy::ff_f>(xsf::gammaincinv)},
        "gammaincinv", gammaincinv_doc);
    PyModule_AddObjectRef(module, "gammaincinv", gammaincinv);

    PyObject *gammaincc = xsf::numpy::ufunc(
        {static_cast<xsf::numpy::dd_d>(xsf::gammaincc), static_cast<xsf::numpy::ff_f>(xsf::gammaincc)}, "gammaincc",
        gammaincc_doc);
    PyModule_AddObjectRef(module, "gammaincc", gammaincc);

    PyObject *gammainccinv = xsf::numpy::ufunc(
        {static_cast<xsf::numpy::dd_d>(xsf::gammainccinv), static_cast<xsf::numpy::ff_f>(xsf::gammainccinv)},
        "gammainccinv", gammainccinv_doc);
    PyModule_AddObjectRef(module, "gammainccinv", gammainccinv);

    PyObject *gammaln =
        xsf::numpy::ufunc({static_cast<xsf::numpy::d_d>(xsf::gammaln), static_cast<xsf::numpy::f_f>(xsf::gammaln)},
                          "gammaln", gammaln_doc);
    PyModule_AddObjectRef(module, "gammaln", gammaln);

    PyObject *gammasgn =
        xsf::numpy::ufunc({static_cast<xsf::numpy::d_d>(xsf::gammasgn), static_cast<xsf::numpy::f_f>(xsf::gammasgn)},
                          "gammasgn", gammasgn_doc);
    PyModule_AddObjectRef(module, "gammasgn", gammasgn);

    PyObject *hyp2f1 =
        xsf::numpy::ufunc({static_cast<xsf::numpy::dddd_d>(xsf::hyp2f1), static_cast<xsf::numpy::dddD_D>(xsf::hyp2f1),
                           static_cast<xsf::numpy::ffff_f>(xsf::hyp2f1), static_cast<xsf::numpy::fffF_F>(xsf::hyp2f1)},
                          "hyp2f1", hyp2f1_doc);
    PyModule_AddObjectRef(module, "hyp2f1", hyp2f1);

    PyObject *hankel1 = xsf::numpy::ufunc(
        {static_cast<xsf::numpy::fF_F>(xsf::cyl_hankel_1), static_cast<xsf::numpy::dD_D>(xsf::cyl_hankel_1)}, "hankel1",
        hankel1_doc);
    PyModule_AddObjectRef(module, "hankel1", hankel1);

    PyObject *hankel1e = xsf::numpy::ufunc(
        {static_cast<xsf::numpy::fF_F>(xsf::cyl_hankel_1e), static_cast<xsf::numpy::dD_D>(xsf::cyl_hankel_1e)},
        "hankel1e", hankel1e_doc);
    PyModule_AddObjectRef(module, "hankel1e", hankel1e);

    PyObject *hankel2 = xsf::numpy::ufunc(
        {static_cast<xsf::numpy::fF_F>(xsf::cyl_hankel_2), static_cast<xsf::numpy::dD_D>(xsf::cyl_hankel_2)}, "hankel2",
        hankel2_doc);
    PyModule_AddObjectRef(module, "hankel2", hankel2);

    PyObject *hankel2e = xsf::numpy::ufunc(
        {static_cast<xsf::numpy::fF_F>(xsf::cyl_hankel_2e), static_cast<xsf::numpy::dD_D>(xsf::cyl_hankel_2e)},
        "hankel2e", hankel2e_doc);
    PyModule_AddObjectRef(module, "hankel2e", hankel2e);

    PyObject *it2i0k0 =
        xsf::numpy::ufunc({static_cast<xsf::numpy::f_ff>(xsf::it2i0k0), static_cast<xsf::numpy::d_dd>(xsf::it2i0k0)}, 2,
                          "it2i0k0", it2i0k0_doc);
    PyModule_AddObjectRef(module, "it2i0k0", it2i0k0);

    PyObject *it2j0y0 =
        xsf::numpy::ufunc({static_cast<xsf::numpy::f_ff>(xsf::it2j0y0), static_cast<xsf::numpy::d_dd>(xsf::it2j0y0)}, 2,
                          "it2j0y0", it2j0y0_doc);
    PyModule_AddObjectRef(module, "it2j0y0", it2j0y0);

    PyObject *it2struve0 = xsf::numpy::ufunc(
        {static_cast<xsf::numpy::f_f>(xsf::it2struve0), static_cast<xsf::numpy::d_d>(xsf::it2struve0)}, "it2struve0",
        it2struve0_doc);
    PyModule_AddObjectRef(module, "it2struve0", it2struve0);

    PyObject *itairy =
        xsf::numpy::ufunc({static_cast<xsf::numpy::f_ffff>(xsf::itairy), static_cast<xsf::numpy::d_dddd>(xsf::itairy)},
                          4, "itairy", itairy_doc);
    PyModule_AddObjectRef(module, "itairy", itairy);

    PyObject *iti0k0 =
        xsf::numpy::ufunc({static_cast<xsf::numpy::f_ff>(xsf::it1i0k0), static_cast<xsf::numpy::d_dd>(xsf::it1i0k0)}, 2,
                          "iti0k0", iti0k0_doc);
    PyModule_AddObjectRef(module, "iti0k0", iti0k0);

    PyObject *itj0y0 =
        xsf::numpy::ufunc({static_cast<xsf::numpy::f_ff>(xsf::it1j0y0), static_cast<xsf::numpy::d_dd>(xsf::it1j0y0)}, 2,
                          "itj0y0", itj0y0_doc);
    PyModule_AddObjectRef(module, "itj0y0", itj0y0);

    PyObject *itmodstruve0 = xsf::numpy::ufunc(
        {static_cast<xsf::numpy::f_f>(xsf::itmodstruve0), static_cast<xsf::numpy::d_d>(xsf::itmodstruve0)},
        "itmodstruve0", itmodstruve0_doc);
    PyModule_AddObjectRef(module, "itmodstruve0", itmodstruve0);

    PyObject *itstruve0 =
        xsf::numpy::ufunc({static_cast<xsf::numpy::f_f>(xsf::itstruve0), static_cast<xsf::numpy::d_d>(xsf::itstruve0)},
                          "itstruve0", itstruve0_doc);
    PyModule_AddObjectRef(module, "itstruve0", itstruve0);

    PyObject *i0 = xsf::numpy::ufunc(
        {static_cast<xsf::numpy::d_d>(xsf::cyl_bessel_i0), static_cast<xsf::numpy::f_f>(xsf::cyl_bessel_i0)}, "i0",
        i0_doc);
    PyModule_AddObjectRef(module, "i0", i0);

    PyObject *i0e = xsf::numpy::ufunc(
        {static_cast<xsf::numpy::d_d>(xsf::cyl_bessel_i0e), static_cast<xsf::numpy::f_f>(xsf::cyl_bessel_i0e)}, "i0e",
        i0e_doc);
    PyModule_AddObjectRef(module, "i0e", i0e);

    PyObject *i1 = xsf::numpy::ufunc(
        {static_cast<xsf::numpy::d_d>(xsf::cyl_bessel_i1), static_cast<xsf::numpy::f_f>(xsf::cyl_bessel_i1)}, "i1",
        i1_doc);
    PyModule_AddObjectRef(module, "i1", i1);

    PyObject *i1e = xsf::numpy::ufunc(
        {static_cast<xsf::numpy::d_d>(xsf::cyl_bessel_i1e), static_cast<xsf::numpy::f_f>(xsf::cyl_bessel_i1e)}, "i1e",
        i1e_doc);
    PyModule_AddObjectRef(module, "i1e", i1e);

    PyObject *iv = xsf::numpy::ufunc(
        {static_cast<xsf::numpy::ff_f>(xsf::cyl_bessel_i), static_cast<xsf::numpy::dd_d>(xsf::cyl_bessel_i),
         static_cast<xsf::numpy::fF_F>(xsf::cyl_bessel_i), static_cast<xsf::numpy::dD_D>(xsf::cyl_bessel_i)},
        "iv", iv_doc);
    PyModule_AddObjectRef(module, "iv", iv);

    PyObject *iv_ratio =
        xsf::numpy::ufunc({static_cast<xsf::numpy::dd_d>(xsf::iv_ratio), static_cast<xsf::numpy::ff_f>(xsf::iv_ratio)},
                          "_iv_ratio", iv_ratio_doc);
    PyModule_AddObjectRef(module, "_iv_ratio", iv_ratio);

    PyObject *iv_ratio_c = xsf::numpy::ufunc(
        {static_cast<xsf::numpy::dd_d>(xsf::iv_ratio_c), static_cast<xsf::numpy::ff_f>(xsf::iv_ratio_c)}, "_iv_ratio_c",
        iv_ratio_c_doc);
    PyModule_AddObjectRef(module, "_iv_ratio_c", iv_ratio_c);

    PyObject *ive = xsf::numpy::ufunc(
        {static_cast<xsf::numpy::ff_f>(xsf::cyl_bessel_ie), static_cast<xsf::numpy::dd_d>(xsf::cyl_bessel_ie),
         static_cast<xsf::numpy::fF_F>(xsf::cyl_bessel_ie), static_cast<xsf::numpy::dD_D>(xsf::cyl_bessel_ie)},
        "ive", ive_doc);
    PyModule_AddObjectRef(module, "ive", ive);

    PyObject *j0 = xsf::numpy::ufunc(
        {static_cast<xsf::numpy::d_d>(xsf::cyl_bessel_j0), static_cast<xsf::numpy::f_f>(xsf::cyl_bessel_j0)}, "j0",
        j0_doc);
    PyModule_AddObjectRef(module, "j0", j0);

    PyObject *j1 = xsf::numpy::ufunc(
        {static_cast<xsf::numpy::d_d>(xsf::cyl_bessel_j1), static_cast<xsf::numpy::f_f>(xsf::cyl_bessel_j1)}, "j1",
        j1_doc);
    PyModule_AddObjectRef(module, "j1", j1);

    PyObject *jv = xsf::numpy::ufunc(
        {static_cast<xsf::numpy::ff_f>(xsf::cyl_bessel_j), static_cast<xsf::numpy::dd_d>(xsf::cyl_bessel_j),
         static_cast<xsf::numpy::fF_F>(xsf::cyl_bessel_j), static_cast<xsf::numpy::dD_D>(xsf::cyl_bessel_j)},
        "jv", jv_doc);
    PyModule_AddObjectRef(module, "jv", jv);

    PyObject *jve = xsf::numpy::ufunc(
        {static_cast<xsf::numpy::ff_f>(xsf::cyl_bessel_je), static_cast<xsf::numpy::dd_d>(xsf::cyl_bessel_je),
         static_cast<xsf::numpy::fF_F>(xsf::cyl_bessel_je), static_cast<xsf::numpy::dD_D>(xsf::cyl_bessel_je)},
        "jve", jve_doc);
    PyModule_AddObjectRef(module, "jve", jve);

    PyObject *kei = xsf::numpy::ufunc({static_cast<xsf::numpy::f_f>(xsf::kei), static_cast<xsf::numpy::d_d>(xsf::kei)},
                                      "kei", kei_doc);
    PyModule_AddObjectRef(module, "kei", kei);

    PyObject *keip = xsf::numpy::ufunc(
        {static_cast<xsf::numpy::f_f>(xsf::keip), static_cast<xsf::numpy::d_d>(xsf::keip)}, "keip", keip_doc);
    PyModule_AddObjectRef(module, "keip", keip);

    PyObject *kelvin =
        xsf::numpy::ufunc({static_cast<xsf::numpy::f_FFFF>(xsf::kelvin), static_cast<xsf::numpy::d_DDDD>(xsf::kelvin)},
                          4, "kelvin", kelvin_doc);
    PyModule_AddObjectRef(module, "kelvin", kelvin);

    PyObject *ker = xsf::numpy::ufunc({static_cast<xsf::numpy::f_f>(xsf::ker), static_cast<xsf::numpy::d_d>(xsf::ker)},
                                      "ker", ker_doc);
    PyModule_AddObjectRef(module, "ker", ker);

    PyObject *kerp = xsf::numpy::ufunc(
        {static_cast<xsf::numpy::f_f>(xsf::kerp), static_cast<xsf::numpy::d_d>(xsf::kerp)}, "kerp", kerp_doc);
    PyModule_AddObjectRef(module, "kerp", kerp);

    PyObject *k0 = xsf::numpy::ufunc(
        {static_cast<xsf::numpy::d_d>(xsf::cyl_bessel_k0), static_cast<xsf::numpy::f_f>(xsf::cyl_bessel_k0)}, "k0",
        k0_doc);
    PyModule_AddObjectRef(module, "k0", k0);

    PyObject *k0e = xsf::numpy::ufunc(
        {static_cast<xsf::numpy::d_d>(xsf::cyl_bessel_k0e), static_cast<xsf::numpy::f_f>(xsf::cyl_bessel_k0e)}, "k0e",
        k0e_doc);
    PyModule_AddObjectRef(module, "k0e", k0e);

    PyObject *k1 = xsf::numpy::ufunc(
        {static_cast<xsf::numpy::d_d>(xsf::cyl_bessel_k1), static_cast<xsf::numpy::f_f>(xsf::cyl_bessel_k1)}, "k1",
        k1_doc);
    PyModule_AddObjectRef(module, "k1", k1);

    PyObject *k1e = xsf::numpy::ufunc(
        {static_cast<xsf::numpy::d_d>(xsf::cyl_bessel_k1e), static_cast<xsf::numpy::f_f>(xsf::cyl_bessel_k1e)}, "k1e",
        k1e_doc);
    PyModule_AddObjectRef(module, "k1e", k1e);

    PyObject *kv = xsf::numpy::ufunc(
        {static_cast<xsf::numpy::ff_f>(xsf::cyl_bessel_k), static_cast<xsf::numpy::dd_d>(xsf::cyl_bessel_k),
         static_cast<xsf::numpy::fF_F>(xsf::cyl_bessel_k), static_cast<xsf::numpy::dD_D>(xsf::cyl_bessel_k)},
        "kv", kv_doc);
    PyModule_AddObjectRef(module, "kv", kv);

    PyObject *kve = xsf::numpy::ufunc(
        {static_cast<xsf::numpy::ff_f>(xsf::cyl_bessel_ke), static_cast<xsf::numpy::dd_d>(xsf::cyl_bessel_ke),
         static_cast<xsf::numpy::fF_F>(xsf::cyl_bessel_ke), static_cast<xsf::numpy::dD_D>(xsf::cyl_bessel_ke)},
        "kve", kve_doc);
    PyModule_AddObjectRef(module, "kve", kve);

    PyObject *log1p =
        xsf::numpy::ufunc({static_cast<xsf::numpy::d_d>(xsf::log1p), static_cast<xsf::numpy::f_f>(xsf::log1p),
                           static_cast<xsf::numpy::D_D>(xsf::log1p), static_cast<xsf::numpy::F_F>(xsf::log1p)},
                          "log1p", log1p_doc);
    PyModule_AddObjectRef(module, "log1p", log1p);

    PyObject *_log1mexp =
        xsf::numpy::ufunc({static_cast<xsf::numpy::d_d>(xsf::log1mexp), static_cast<xsf::numpy::f_f>(xsf::log1mexp)},
                          "_log1mexp", _log1mexp_doc);
    PyModule_AddObjectRef(module, "_log1mexp", _log1mexp);

    PyObject *_log1pmx =
        xsf::numpy::ufunc({static_cast<xsf::numpy::d_d>(xsf::log1pmx), static_cast<xsf::numpy::f_f>(xsf::log1pmx)},
                          "_log1pmx", _log1pmx_doc);
    PyModule_AddObjectRef(module, "_log1pmx", _log1pmx);

    PyObject *xlogy =
        xsf::numpy::ufunc({static_cast<xsf::numpy::dd_d>(xsf::xlogy), static_cast<xsf::numpy::ff_f>(xsf::xlogy),
                           static_cast<xsf::numpy::DD_D>(xsf::xlogy), static_cast<xsf::numpy::FF_F>(xsf::xlogy)},
                          "xlogy", xlogy_doc);
    PyModule_AddObjectRef(module, "xlogy", xlogy);

    PyObject *xlog1py =
        xsf::numpy::ufunc({static_cast<xsf::numpy::dd_d>(xsf::xlog1py), static_cast<xsf::numpy::ff_f>(xsf::xlog1py),
                           static_cast<xsf::numpy::DD_D>(xsf::xlog1py), static_cast<xsf::numpy::FF_F>(xsf::xlog1py)},
                          "xlog1py", xlog1py_doc);
    PyModule_AddObjectRef(module, "xlog1py", xlog1py);

    PyObject *log_expit =
        xsf::numpy::ufunc({static_cast<xsf::numpy::d_d>(xsf::log_expit), static_cast<xsf::numpy::f_f>(xsf::log_expit),
                           static_cast<xsf::numpy::g_g>(xsf::log_expit)},
                          "log_expit", log_expit_doc);
    PyModule_AddObjectRef(module, "log_expit", log_expit);

    PyObject *log_wright_bessel = xsf::numpy::ufunc({static_cast<xsf::numpy::ddd_d>(xsf::log_wright_bessel),
                                                     static_cast<xsf::numpy::fff_f>(xsf::log_wright_bessel)},
                                                    "log_wright_bessel", log_wright_bessel_doc);
    PyModule_AddObjectRef(module, "log_wright_bessel", log_wright_bessel);

    PyObject *logit =
        xsf::numpy::ufunc({static_cast<xsf::numpy::d_d>(xsf::logit), static_cast<xsf::numpy::f_f>(xsf::logit),
                           static_cast<xsf::numpy::g_g>(xsf::logit)},
                          "logit", logit_doc);
    PyModule_AddObjectRef(module, "logit", logit);

    PyObject *loggamma =
        xsf::numpy::ufunc({static_cast<xsf::numpy::d_d>(xsf::loggamma), static_cast<xsf::numpy::D_D>(xsf::loggamma),
                           static_cast<xsf::numpy::f_f>(xsf::loggamma), static_cast<xsf::numpy::F_F>(xsf::loggamma)},
                          "loggamma", loggamma_doc);
    PyModule_AddObjectRef(module, "loggamma", loggamma);

    PyObject *legendre_p = Py_BuildValue(
        "(N, N, N)",
        xsf::numpy::gufunc({xsf::numpy::compose{xsf::numpy::autodiff(), xsf::numpy::use_long_long_int()},
                            static_cast<xsf::numpy::autodiff0_id_d>(xsf::legendre_p),
                            static_cast<xsf::numpy::autodiff0_if_f>(xsf::legendre_p)},
                           "legendre_p", nullptr, "(),()->(1)", [](const npy_intp *dims, npy_intp *new_dims) {}),
        xsf::numpy::gufunc({xsf::numpy::compose{xsf::numpy::autodiff(), xsf::numpy::use_long_long_int()},
                            static_cast<xsf::numpy::autodiff1_id_d>(xsf::legendre_p),
                            static_cast<xsf::numpy::autodiff1_if_f>(xsf::legendre_p)},
                           "legendre_p", nullptr, "(),()->(2)", [](const npy_intp *dims, npy_intp *new_dims) {}),
        xsf::numpy::gufunc({xsf::numpy::compose{xsf::numpy::autodiff(), xsf::numpy::use_long_long_int()},
                            static_cast<xsf::numpy::autodiff2_id_d>(xsf::legendre_p),
                            static_cast<xsf::numpy::autodiff2_if_f>(xsf::legendre_p)},
                           "legendre_p", nullptr, "(),()->(3)", [](const npy_intp *dims, npy_intp *new_dims) {}));
    PyModule_AddObjectRef(module, "legendre_p", legendre_p);

    PyObject *assoc_legendre_p = Py_BuildValue(
        "{(O, i): N, (O, i): N, (O, i): N, (O, i): N, (O, i): N,(O, i): N}", Py_True, 0,
        xsf::numpy::gufunc({xsf::numpy::compose{xsf::numpy::autodiff(), xsf::numpy::use_long_long_int()},
                            [](int n, int m, xsf::dual<double, 0> z, int branch_cut) {
                                return xsf::assoc_legendre_p(xsf::assoc_legendre_norm, n, m, z, branch_cut);
                            },
                            [](int n, int m, xsf::dual<float, 0> z, int branch_cut) {
                                return xsf::assoc_legendre_p(xsf::assoc_legendre_norm, n, m, z, branch_cut);
                            },
                            [](int n, int m, xsf::dual<xsf::numpy::cdouble, 0> z, int branch_cut) {
                                return xsf::assoc_legendre_p(xsf::assoc_legendre_norm, n, m, z, branch_cut);
                            },
                            [](int n, int m, xsf::dual<xsf::numpy::cfloat, 0> z, int branch_cut) {
                                return xsf::assoc_legendre_p(xsf::assoc_legendre_norm, n, m, z, branch_cut);
                            }},
                           1, "assoc_legendre_p", nullptr, "(),(),(),()->(1)",
                           [](const npy_intp *dims, npy_intp *new_dims) {}),
        Py_True, 1,
        xsf::numpy::gufunc({xsf::numpy::compose{xsf::numpy::autodiff(), xsf::numpy::use_long_long_int()},
                            [](int n, int m, xsf::dual<double, 1> z, int branch_cut) {
                                return xsf::assoc_legendre_p(xsf::assoc_legendre_norm, n, m, z, branch_cut);
                            },
                            [](int n, int m, xsf::dual<float, 1> z, int branch_cut) {
                                return xsf::assoc_legendre_p(xsf::assoc_legendre_norm, n, m, z, branch_cut);
                            },
                            [](int n, int m, xsf::dual<xsf::numpy::cdouble, 1> z, int branch_cut) {
                                return xsf::assoc_legendre_p(xsf::assoc_legendre_norm, n, m, z, branch_cut);
                            },
                            [](int n, int m, xsf::dual<xsf::numpy::cfloat, 1> z, int branch_cut) {
                                return xsf::assoc_legendre_p(xsf::assoc_legendre_norm, n, m, z, branch_cut);
                            }},
                           1, "assoc_legendre_p", nullptr, "(),(),(),()->(2)",
                           [](const npy_intp *dims, npy_intp *new_dims) {}),
        Py_True, 2,
        xsf::numpy::gufunc({xsf::numpy::compose{xsf::numpy::autodiff(), xsf::numpy::use_long_long_int()},
                            [](int n, int m, xsf::dual<double, 2> z, int branch_cut) {
                                return xsf::assoc_legendre_p(xsf::assoc_legendre_norm, n, m, z, branch_cut);
                            },
                            [](int n, int m, xsf::dual<float, 2> z, int branch_cut) {
                                return xsf::assoc_legendre_p(xsf::assoc_legendre_norm, n, m, z, branch_cut);
                            },
                            [](int n, int m, xsf::dual<xsf::numpy::cdouble, 2> z, int branch_cut) {
                                return xsf::assoc_legendre_p(xsf::assoc_legendre_norm, n, m, z, branch_cut);
                            },
                            [](int n, int m, xsf::dual<xsf::numpy::cfloat, 2> z, int branch_cut) {
                                return xsf::assoc_legendre_p(xsf::assoc_legendre_norm, n, m, z, branch_cut);
                            }},
                           1, "assoc_legendre_p", nullptr, "(),(),(),()->(3)",
                           [](const npy_intp *dims, npy_intp *new_dims) {}),
        Py_False, 0,
        xsf::numpy::gufunc({xsf::numpy::compose{xsf::numpy::autodiff(), xsf::numpy::use_long_long_int()},
                            [](int n, int m, xsf::dual<double, 0> z, int branch_cut) {
                                return xsf::assoc_legendre_p(xsf::assoc_legendre_unnorm, n, m, z, branch_cut);
                            },
                            [](int n, int m, xsf::dual<float, 0> z, int branch_cut) {
                                return xsf::assoc_legendre_p(xsf::assoc_legendre_unnorm, n, m, z, branch_cut);
                            },
                            [](int n, int m, xsf::dual<xsf::numpy::cdouble, 0> z, int branch_cut) {
                                return xsf::assoc_legendre_p(xsf::assoc_legendre_unnorm, n, m, z, branch_cut);
                            },
                            [](int n, int m, xsf::dual<xsf::numpy::cfloat, 0> z, int branch_cut) {
                                return xsf::assoc_legendre_p(xsf::assoc_legendre_unnorm, n, m, z, branch_cut);
                            }},
                           1, "assoc_legendre_p", nullptr, "(),(),(),()->(1)",
                           [](const npy_intp *dims, npy_intp *new_dims) {}),
        Py_False, 1,
        xsf::numpy::gufunc({xsf::numpy::compose{xsf::numpy::autodiff(), xsf::numpy::use_long_long_int()},
                            [](int n, int m, xsf::dual<double, 1> z, int branch_cut) {
                                return xsf::assoc_legendre_p(xsf::assoc_legendre_unnorm, n, m, z, branch_cut);
                            },
                            [](int n, int m, xsf::dual<float, 1> z, int branch_cut) {
                                return xsf::assoc_legendre_p(xsf::assoc_legendre_unnorm, n, m, z, branch_cut);
                            },
                            [](int n, int m, xsf::dual<xsf::numpy::cdouble, 1> z, int branch_cut) {
                                return xsf::assoc_legendre_p(xsf::assoc_legendre_unnorm, n, m, z, branch_cut);
                            },
                            [](int n, int m, xsf::dual<xsf::numpy::cfloat, 1> z, int branch_cut) {
                                return xsf::assoc_legendre_p(xsf::assoc_legendre_unnorm, n, m, z, branch_cut);
                            }},
                           1, "assoc_legendre_p", nullptr, "(),(),(),()->(2)",
                           [](const npy_intp *dims, npy_intp *new_dims) {}),
        Py_False, 2,
        xsf::numpy::gufunc({xsf::numpy::compose{xsf::numpy::autodiff(), xsf::numpy::use_long_long_int()},
                            [](int n, int m, xsf::dual<double, 2> z, int branch_cut) {
                                return xsf::assoc_legendre_p(xsf::assoc_legendre_unnorm, n, m, z, branch_cut);
                            },
                            [](int n, int m, xsf::dual<float, 2> z, int branch_cut) {
                                return xsf::assoc_legendre_p(xsf::assoc_legendre_unnorm, n, m, z, branch_cut);
                            },
                            [](int n, int m, xsf::dual<xsf::numpy::cdouble, 2> z, int branch_cut) {
                                return xsf::assoc_legendre_p(xsf::assoc_legendre_unnorm, n, m, z, branch_cut);
                            },
                            [](int n, int m, xsf::dual<xsf::numpy::cfloat, 2> z, int branch_cut) {
                                return xsf::assoc_legendre_p(xsf::assoc_legendre_unnorm, n, m, z, branch_cut);
                            }},
                           1, "assoc_legendre_p", nullptr, "(),(),(),()->(3)",
                           [](const npy_intp *dims, npy_intp *new_dims) {}));
    PyModule_AddObjectRef(module, "assoc_legendre_p", assoc_legendre_p);

    PyObject *mathieu_a =
        xsf::numpy::ufunc({static_cast<xsf::numpy::ff_f>(xsf::cem_cva), static_cast<xsf::numpy::dd_d>(xsf::cem_cva)},
                          "mathieu_a", mathieu_a_doc);
    PyModule_AddObjectRef(module, "mathieu_a", mathieu_a);

    PyObject *mathieu_b =
        xsf::numpy::ufunc({static_cast<xsf::numpy::ff_f>(xsf::sem_cva), static_cast<xsf::numpy::dd_d>(xsf::sem_cva)},
                          "mathieu_b", mathieu_b_doc);
    PyModule_AddObjectRef(module, "mathieu_b", mathieu_b);

    PyObject *mathieu_cem =
        xsf::numpy::ufunc({static_cast<xsf::numpy::fff_ff>(xsf::cem), static_cast<xsf::numpy::ddd_dd>(xsf::cem)}, 2,
                          "mathieu_cem", mathieu_cem_doc);
    PyModule_AddObjectRef(module, "mathieu_cem", mathieu_cem);

    PyObject *mathieu_modcem1 =
        xsf::numpy::ufunc({static_cast<xsf::numpy::fff_ff>(xsf::mcm1), static_cast<xsf::numpy::ddd_dd>(xsf::mcm1)}, 2,
                          "mathieu_modcem1", mathieu_modcem1_doc);
    PyModule_AddObjectRef(module, "mathieu_modcem1", mathieu_modcem1);

    PyObject *mathieu_modcem2 =
        xsf::numpy::ufunc({static_cast<xsf::numpy::fff_ff>(xsf::mcm2), static_cast<xsf::numpy::ddd_dd>(xsf::mcm2)}, 2,
                          "mathieu_modcem2", mathieu_modcem2_doc);
    PyModule_AddObjectRef(module, "mathieu_modcem2", mathieu_modcem2);

    PyObject *mathieu_modsem1 =
        xsf::numpy::ufunc({static_cast<xsf::numpy::fff_ff>(xsf::msm1), static_cast<xsf::numpy::ddd_dd>(xsf::msm1)}, 2,
                          "mathieu_modsem1", mathieu_modsem1_doc);
    PyModule_AddObjectRef(module, "mathieu_modsem1", mathieu_modsem1);

    PyObject *mathieu_modsem2 =
        xsf::numpy::ufunc({static_cast<xsf::numpy::fff_ff>(xsf::msm2), static_cast<xsf::numpy::ddd_dd>(xsf::msm2)}, 2,
                          "mathieu_modsem2", mathieu_modsem2_doc);
    PyModule_AddObjectRef(module, "mathieu_modsem2", mathieu_modsem2);

    PyObject *mathieu_sem =
        xsf::numpy::ufunc({static_cast<xsf::numpy::fff_ff>(xsf::sem), static_cast<xsf::numpy::ddd_dd>(xsf::sem)}, 2,
                          "mathieu_sem", mathieu_sem_doc);
    PyModule_AddObjectRef(module, "mathieu_sem", mathieu_sem);

    PyObject *modfresnelm = xsf::numpy::ufunc({static_cast<xsf::numpy::f_FF>(xsf::modified_fresnel_minus),
                                               static_cast<xsf::numpy::d_DD>(xsf::modified_fresnel_minus)},
                                              2, "modfresnelm", modfresnelm_doc);
    PyModule_AddObjectRef(module, "modfresnelm", modfresnelm);

    PyObject *modfresnelp = xsf::numpy::ufunc({static_cast<xsf::numpy::f_FF>(xsf::modified_fresnel_plus),
                                               static_cast<xsf::numpy::d_DD>(xsf::modified_fresnel_plus)},
                                              2, "modfresnelp", modfresnelp_doc);
    PyModule_AddObjectRef(module, "modfresnelp", modfresnelp);

    PyObject *modstruve =
        xsf::numpy::ufunc({static_cast<xsf::numpy::dd_d>(xsf::struve_l), static_cast<xsf::numpy::ff_f>(xsf::struve_l)},
                          "modstruve", struve_l_doc);
    PyModule_AddObjectRef(module, "modstruve", modstruve);

    PyObject *obl_ang1 = xsf::numpy::ufunc({static_cast<xsf::numpy::ffff_ff>(xsf::oblate_aswfa_nocv),
                                            static_cast<xsf::numpy::dddd_dd>(xsf::oblate_aswfa_nocv)},
                                           2, "obl_ang1", obl_ang1_doc);
    PyModule_AddObjectRef(module, "obl_ang1", obl_ang1);

    PyObject *obl_ang1_cv = xsf::numpy::ufunc(
        {static_cast<xsf::numpy::fffff_ff>(xsf::oblate_aswfa), static_cast<xsf::numpy::ddddd_dd>(xsf::oblate_aswfa)}, 2,
        "obl_ang1_cv", obl_ang1_cv_doc);
    PyModule_AddObjectRef(module, "obl_ang1_cv", obl_ang1_cv);

    PyObject *obl_cv = xsf::numpy::ufunc(
        {static_cast<xsf::numpy::fff_f>(xsf::oblate_segv), static_cast<xsf::numpy::ddd_d>(xsf::oblate_segv)}, "obl_cv",
        obl_cv_doc);
    PyModule_AddObjectRef(module, "obl_cv", obl_cv);

    PyObject *obl_rad1 = xsf::numpy::ufunc({static_cast<xsf::numpy::ffff_ff>(xsf::oblate_radial1_nocv),
                                            static_cast<xsf::numpy::dddd_dd>(xsf::oblate_radial1_nocv)},
                                           2, "obl_rad1", obl_rad1_doc);
    PyModule_AddObjectRef(module, "obl_rad1", obl_rad1);

    PyObject *obl_rad1_cv = xsf::numpy::ufunc({static_cast<xsf::numpy::fffff_ff>(xsf::oblate_radial1),
                                               static_cast<xsf::numpy::ddddd_dd>(xsf::oblate_radial1)},
                                              2, "obl_rad1_cv", obl_rad1_cv_doc);
    PyModule_AddObjectRef(module, "obl_rad1_cv", obl_rad1_cv);

    PyObject *obl_rad2 = xsf::numpy::ufunc({static_cast<xsf::numpy::ffff_ff>(xsf::oblate_radial2_nocv),
                                            static_cast<xsf::numpy::dddd_dd>(xsf::oblate_radial2_nocv)},
                                           2, "obl_rad2", obl_rad2_doc);
    PyModule_AddObjectRef(module, "obl_rad2", obl_rad2);

    PyObject *obl_rad2_cv = xsf::numpy::ufunc({static_cast<xsf::numpy::fffff_ff>(xsf::oblate_radial2),
                                               static_cast<xsf::numpy::ddddd_dd>(xsf::oblate_radial2)},
                                              2, "obl_rad2_cv", obl_rad2_cv_doc);
    PyModule_AddObjectRef(module, "obl_rad2_cv", obl_rad2_cv);

    PyObject *pbdv = xsf::numpy::ufunc(
        {static_cast<xsf::numpy::ff_ff>(xsf::pbdv), static_cast<xsf::numpy::dd_dd>(xsf::pbdv)}, 2, "pbdv", pbdv_doc);
    PyModule_AddObjectRef(module, "pbdv", pbdv);

    PyObject *pbvv = xsf::numpy::ufunc(
        {static_cast<xsf::numpy::ff_ff>(xsf::pbvv), static_cast<xsf::numpy::dd_dd>(xsf::pbvv)}, 2, "pbvv", pbvv_doc);
    PyModule_AddObjectRef(module, "pbvv", pbvv);

    PyObject *pbwa = xsf::numpy::ufunc(
        {static_cast<xsf::numpy::ff_ff>(xsf::pbwa), static_cast<xsf::numpy::dd_dd>(xsf::pbwa)}, 2, "pbwa", pbwa_doc);
    PyModule_AddObjectRef(module, "pbwa", pbwa);

    PyObject *pro_ang1 = xsf::numpy::ufunc({static_cast<xsf::numpy::ffff_ff>(xsf::prolate_aswfa_nocv),
                                            static_cast<xsf::numpy::dddd_dd>(xsf::prolate_aswfa_nocv)},
                                           2, "pro_ang1", pro_ang1_doc);
    PyModule_AddObjectRef(module, "pro_ang1", pro_ang1);

    PyObject *pro_ang1_cv = xsf::numpy::ufunc(
        {static_cast<xsf::numpy::fffff_ff>(xsf::prolate_aswfa), static_cast<xsf::numpy::ddddd_dd>(xsf::prolate_aswfa)},
        2, "pro_ang1_cv", pro_ang1_cv_doc);
    PyModule_AddObjectRef(module, "pro_ang1_cv", pro_ang1_cv);

    PyObject *pro_cv = xsf::numpy::ufunc(
        {static_cast<xsf::numpy::fff_f>(xsf::prolate_segv), static_cast<xsf::numpy::ddd_d>(xsf::prolate_segv)},
        "pro_cv", pro_cv_doc);
    PyModule_AddObjectRef(module, "pro_cv", pro_cv);

    PyObject *pro_rad1 = xsf::numpy::ufunc({static_cast<xsf::numpy::ffff_ff>(xsf::prolate_radial1_nocv),
                                            static_cast<xsf::numpy::dddd_dd>(xsf::prolate_radial1_nocv)},
                                           2, "pro_rad1", pro_rad1_doc);
    PyModule_AddObjectRef(module, "pro_rad1", pro_rad1);

    PyObject *pro_rad1_cv = xsf::numpy::ufunc({static_cast<xsf::numpy::fffff_ff>(xsf::prolate_radial1),
                                               static_cast<xsf::numpy::ddddd_dd>(xsf::prolate_radial1)},
                                              2, "pro_rad1_cv", pro_rad1_cv_doc);
    PyModule_AddObjectRef(module, "pro_rad1_cv", pro_rad1_cv);

    PyObject *pro_rad2 = xsf::numpy::ufunc({static_cast<xsf::numpy::ffff_ff>(xsf::prolate_radial2_nocv),
                                            static_cast<xsf::numpy::dddd_dd>(xsf::prolate_radial2_nocv)},
                                           2, "pro_rad2", pro_rad2_doc);
    PyModule_AddObjectRef(module, "pro_rad2", pro_rad2);

    PyObject *pro_rad2_cv = xsf::numpy::ufunc({static_cast<xsf::numpy::fffff_ff>(xsf::prolate_radial2),
                                               static_cast<xsf::numpy::ddddd_dd>(xsf::prolate_radial2)},
                                              2, "pro_rad2_cv", pro_rad2_cv_doc);
    PyModule_AddObjectRef(module, "pro_rad2_cv", pro_rad2_cv);

    PyObject *psi =
        xsf::numpy::ufunc({static_cast<xsf::numpy::d_d>(xsf::digamma), static_cast<xsf::numpy::D_D>(xsf::digamma),
                           static_cast<xsf::numpy::f_f>(xsf::digamma), static_cast<xsf::numpy::F_F>(xsf::digamma)},
                          "psi", psi_doc);
    PyModule_AddObjectRef(module, "psi", psi);

    PyObject *radian =
        xsf::numpy::ufunc({static_cast<xsf::numpy::ddd_d>(xsf::radian), static_cast<xsf::numpy::fff_f>(xsf::radian)},
                          "radian", radian_doc);
    PyModule_AddObjectRef(module, "radian", radian);

    PyObject *rgamma =
        xsf::numpy::ufunc({static_cast<xsf::numpy::d_d>(xsf::rgamma), static_cast<xsf::numpy::D_D>(xsf::rgamma),
                           static_cast<xsf::numpy::f_f>(xsf::rgamma), static_cast<xsf::numpy::F_F>(xsf::rgamma)},
                          "rgamma", rgamma_doc);
    PyModule_AddObjectRef(module, "rgamma", rgamma);

    PyObject *_riemann_zeta = xsf::numpy::ufunc(
        {static_cast<xsf::numpy::d_d>(xsf::riemann_zeta), static_cast<xsf::numpy::D_D>(xsf::riemann_zeta),
         static_cast<xsf::numpy::f_f>(xsf::riemann_zeta), static_cast<xsf::numpy::F_F>(xsf::riemann_zeta)},
        "_riemann_zeta", _riemann_zeta_doc);
    PyModule_AddObjectRef(module, "_riemann_zeta", _riemann_zeta);

    PyObject *sindg = xsf::numpy::ufunc(
        {static_cast<xsf::numpy::d_d>(xsf::sindg), static_cast<xsf::numpy::f_f>(xsf::sindg)}, "sindg", sindg_doc);
    PyModule_AddObjectRef(module, "sindg", sindg);

    PyObject *_spherical_jn = xsf::numpy::ufunc(
        {static_cast<xsf::numpy::ld_d>(xsf::sph_bessel_j), static_cast<xsf::numpy::lD_D>(xsf::sph_bessel_j),
         static_cast<xsf::numpy::lf_f>(xsf::sph_bessel_j), static_cast<xsf::numpy::lF_F>(xsf::sph_bessel_j)},
        "_spherical_jn", spherical_jn_doc);
    PyModule_AddObjectRef(module, "_spherical_jn", _spherical_jn);

    PyObject *_spherical_jn_d = xsf::numpy::ufunc(
        {static_cast<xsf::numpy::ld_d>(xsf::sph_bessel_j_jac), static_cast<xsf::numpy::lD_D>(xsf::sph_bessel_j_jac),
         static_cast<xsf::numpy::lf_f>(xsf::sph_bessel_j_jac), static_cast<xsf::numpy::lF_F>(xsf::sph_bessel_j_jac)},
        "_spherical_jn_d", spherical_jn_d_doc);
    PyModule_AddObjectRef(module, "_spherical_jn_d", _spherical_jn_d);

    PyObject *_spherical_yn = xsf::numpy::ufunc(
        {static_cast<xsf::numpy::ld_d>(xsf::sph_bessel_y), static_cast<xsf::numpy::lD_D>(xsf::sph_bessel_y),
         static_cast<xsf::numpy::lf_f>(xsf::sph_bessel_y), static_cast<xsf::numpy::lF_F>(xsf::sph_bessel_y)},
        "_spherical_yn", spherical_yn_doc);
    PyModule_AddObjectRef(module, "_spherical_yn", _spherical_yn);

    PyObject *_spherical_yn_d = xsf::numpy::ufunc(
        {static_cast<xsf::numpy::ld_d>(xsf::sph_bessel_y_jac), static_cast<xsf::numpy::lD_D>(xsf::sph_bessel_y_jac),
         static_cast<xsf::numpy::lf_f>(xsf::sph_bessel_y_jac), static_cast<xsf::numpy::lF_F>(xsf::sph_bessel_y_jac)},
        "_spherical_yn_d", spherical_yn_d_doc);
    PyModule_AddObjectRef(module, "_spherical_yn_d", _spherical_yn_d);

    PyObject *_spherical_in = xsf::numpy::ufunc(
        {static_cast<xsf::numpy::ld_d>(xsf::sph_bessel_i), static_cast<xsf::numpy::lD_D>(xsf::sph_bessel_i),
         static_cast<xsf::numpy::lf_f>(xsf::sph_bessel_i), static_cast<xsf::numpy::lF_F>(xsf::sph_bessel_i)},
        "_spherical_in", spherical_in_doc);
    PyModule_AddObjectRef(module, "_spherical_in", _spherical_in);

    PyObject *_spherical_in_d = xsf::numpy::ufunc(
        {static_cast<xsf::numpy::ld_d>(xsf::sph_bessel_i_jac), static_cast<xsf::numpy::lD_D>(xsf::sph_bessel_i_jac),
         static_cast<xsf::numpy::lf_f>(xsf::sph_bessel_i_jac), static_cast<xsf::numpy::lF_F>(xsf::sph_bessel_i_jac)},
        "_spherical_in_d", spherical_in_d_doc);
    PyModule_AddObjectRef(module, "_spherical_in_d", _spherical_in_d);

    PyObject *_spherical_kn = xsf::numpy::ufunc(
        {static_cast<xsf::numpy::ld_d>(xsf::sph_bessel_k), static_cast<xsf::numpy::lD_D>(xsf::sph_bessel_k),
         static_cast<xsf::numpy::lf_f>(xsf::sph_bessel_k), static_cast<xsf::numpy::lF_F>(xsf::sph_bessel_k)},
        "_spherical_kn", spherical_kn_doc);
    PyModule_AddObjectRef(module, "_spherical_kn", _spherical_kn);

    PyObject *_spherical_kn_d = xsf::numpy::ufunc(
        {static_cast<xsf::numpy::ld_d>(xsf::sph_bessel_k_jac), static_cast<xsf::numpy::lD_D>(xsf::sph_bessel_k_jac),
         static_cast<xsf::numpy::lf_f>(xsf::sph_bessel_k_jac), static_cast<xsf::numpy::lF_F>(xsf::sph_bessel_k_jac)},
        "_spherical_kn_d", spherical_kn_d_doc);
    PyModule_AddObjectRef(module, "_spherical_kn_d", _spherical_kn_d);

    PyObject *sph_legendre_p = Py_BuildValue(
        "(N, N, N)",
        xsf::numpy::gufunc({xsf::numpy::compose{xsf::numpy::autodiff(), xsf::numpy::use_long_long_int()},
                            static_cast<xsf::numpy::autodiff0_iid_d>(xsf::sph_legendre_p),
                            static_cast<xsf::numpy::autodiff0_iif_f>(xsf::sph_legendre_p)},
                           "sph_legendre_p", nullptr, "(),(),()->(1)", [](const npy_intp *dims, npy_intp *new_dims) {}),
        xsf::numpy::gufunc({xsf::numpy::compose{xsf::numpy::autodiff(), xsf::numpy::use_long_long_int()},
                            static_cast<xsf::numpy::autodiff1_iid_d>(xsf::sph_legendre_p),
                            static_cast<xsf::numpy::autodiff1_iif_f>(xsf::sph_legendre_p)},
                           "sph_legendre_p", nullptr, "(),(),()->(2)", [](const npy_intp *dims, npy_intp *new_dims) {}),
        xsf::numpy::gufunc({xsf::numpy::compose{xsf::numpy::autodiff(), xsf::numpy::use_long_long_int()},
                            static_cast<xsf::numpy::autodiff2_iid_d>(xsf::sph_legendre_p),
                            static_cast<xsf::numpy::autodiff2_iif_f>(xsf::sph_legendre_p)},
                           "sph_legendre_p", nullptr, "(),(),()->(3)",
                           [](const npy_intp *dims, npy_intp *new_dims) {}));
    PyModule_AddObjectRef(module, "sph_legendre_p", sph_legendre_p);

    PyObject *sph_harm_y =
        Py_BuildValue("(N, N, N)",
                      xsf::numpy::gufunc({xsf::numpy::compose{xsf::numpy::autodiff(), xsf::numpy::use_long_long_int()},
                                          static_cast<xsf::numpy::autodiff00_iidd_D>(xsf::sph_harm_y),
                                          static_cast<xsf::numpy::autodiff00_iiff_F>(xsf::sph_harm_y)},
                                         "sph_harm_y", nullptr, "(),(),(),()->(1,1)",
                                         [](const npy_intp *dims, npy_intp *new_dims) {}),
                      xsf::numpy::gufunc({xsf::numpy::compose{xsf::numpy::autodiff(), xsf::numpy::use_long_long_int()},
                                          static_cast<xsf::numpy::autodiff11_iidd_D>(xsf::sph_harm_y),
                                          static_cast<xsf::numpy::autodiff11_iiff_F>(xsf::sph_harm_y)},
                                         "sph_harm_y", nullptr, "(),(),(),()->(2,2)",
                                         [](const npy_intp *dims, npy_intp *new_dims) {}),
                      xsf::numpy::gufunc({xsf::numpy::compose{xsf::numpy::autodiff(), xsf::numpy::use_long_long_int()},
                                          static_cast<xsf::numpy::autodiff22_iidd_D>(xsf::sph_harm_y),
                                          static_cast<xsf::numpy::autodiff22_iiff_F>(xsf::sph_harm_y)},
                                         "sph_harm_y", nullptr, "(),(),(),()->(3,3)",
                                         [](const npy_intp *dims, npy_intp *new_dims) {}));
    PyModule_AddObjectRef(module, "sph_harm_y", sph_harm_y);

    PyObject *struve =
        xsf::numpy::ufunc({static_cast<xsf::numpy::dd_d>(xsf::struve_h), static_cast<xsf::numpy::ff_f>(xsf::struve_h)},
                          "struve", struve_h_doc);
    PyModule_AddObjectRef(module, "struve", struve);

    PyObject *tandg = xsf::numpy::ufunc(
        {static_cast<xsf::numpy::d_d>(xsf::tandg), static_cast<xsf::numpy::f_f>(xsf::tandg)}, "tandg", tandg_doc);
    PyModule_AddObjectRef(module, "tandg", tandg);

    PyObject *wright_bessel = xsf::numpy::ufunc(
        {static_cast<xsf::numpy::ddd_d>(xsf::wright_bessel), static_cast<xsf::numpy::fff_f>(xsf::wright_bessel)},
        "wright_bessel", wright_bessel_doc);
    PyModule_AddObjectRef(module, "wright_bessel", wright_bessel);

    PyObject *y0 = xsf::numpy::ufunc(
        {static_cast<xsf::numpy::d_d>(xsf::cyl_bessel_y0), static_cast<xsf::numpy::f_f>(xsf::cyl_bessel_y0)}, "y0",
        y0_doc);
    PyModule_AddObjectRef(module, "y0", y0);

    PyObject *y1 = xsf::numpy::ufunc(
        {static_cast<xsf::numpy::d_d>(xsf::cyl_bessel_y1), static_cast<xsf::numpy::f_f>(xsf::cyl_bessel_y1)}, "y1",
        y1_doc);
    PyModule_AddObjectRef(module, "y1", y1);

    PyObject *yv = xsf::numpy::ufunc(
        {static_cast<xsf::numpy::ff_f>(xsf::cyl_bessel_y), static_cast<xsf::numpy::dd_d>(xsf::cyl_bessel_y),
         static_cast<xsf::numpy::fF_F>(xsf::cyl_bessel_y), static_cast<xsf::numpy::dD_D>(xsf::cyl_bessel_y)},
        "yv", yv_doc);
    PyModule_AddObjectRef(module, "yv", yv);

    PyObject *yve = xsf::numpy::ufunc(
        {static_cast<xsf::numpy::ff_f>(xsf::cyl_bessel_ye), static_cast<xsf::numpy::dd_d>(xsf::cyl_bessel_ye),
         static_cast<xsf::numpy::fF_F>(xsf::cyl_bessel_ye), static_cast<xsf::numpy::dD_D>(xsf::cyl_bessel_ye)},
        "yve", yve_doc);
    PyModule_AddObjectRef(module, "yve", yve);

    PyObject *zetac = xsf::numpy::ufunc(
        {static_cast<xsf::numpy::d_d>(xsf::zetac), static_cast<xsf::numpy::f_f>(xsf::zetac)}, "zetac", zetac_doc);
    PyModule_AddObjectRef(module, "zetac", zetac);

    return 0;
}


static PyModuleDef_Slot _special_ufuncs_slots[] = {
    {Py_mod_exec, (void *)_special_ufuncs_module_exec},
    {Py_mod_multiple_interpreters, Py_MOD_PER_INTERPRETER_GIL_SUPPORTED},
#if PY_VERSION_HEX >= 0x030d00f0  /* Python 3.13+ */
    {Py_mod_gil, Py_MOD_GIL_NOT_USED},
#endif
    {0, NULL},
};


static struct PyModuleDef _special_ufuncs_def = {
    /* m_base     */ PyModuleDef_HEAD_INIT,
    /* m_name     */ "_special_ufuncs",
    /* m_doc      */ NULL,
    /* m_size     */ 0,
    /* m_methods  */ _methods,
    /* m_slots    */ _special_ufuncs_slots,
    /* m_traverse */ NULL,
    /* m_clear    */ NULL,
    /* m_free     */ NULL
};


PyMODINIT_FUNC PyInit__special_ufuncs()
{
    return PyModuleDef_Init(&_special_ufuncs_def);
}
