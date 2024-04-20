#include "ufunc.h"

#include "special.h"
#include "special/bessel.h"
#include "special/legendre.h"
#include "special/sph_harm.h"

using namespace std;

using float1_mdspan = mdspan<float, dextents<ptrdiff_t, 1>, layout_stride>;
using double1_mdspan = mdspan<double, dextents<ptrdiff_t, 1>, layout_stride>;
using complex_float1_mdspan = mdspan<complex<float>, dextents<ptrdiff_t, 1>, layout_stride>;
using complex_double1_mdspan = mdspan<complex<double>, dextents<ptrdiff_t, 1>, layout_stride>;

using func_f_f1_t = void (*)(float, mdspan<float, dextents<ptrdiff_t, 1>, layout_stride>);
using func_d_d1_t = void (*)(double, mdspan<double, dextents<ptrdiff_t, 1>, layout_stride>);
using func_F_F1_t = void (*)(complex<float>, mdspan<complex<float>, dextents<ptrdiff_t, 1>, layout_stride>);
using func_D_D1_t = void (*)(complex<double>, mdspan<complex<double>, dextents<ptrdiff_t, 1>, layout_stride>);

using func_f_f1f1_t = void (*)(float, float1_mdspan, float1_mdspan);
using func_d_d1d1_t =
    void (*)(double, mdspan<double, dextents<ptrdiff_t, 1>, layout_stride>, mdspan<double, dextents<ptrdiff_t, 1>, layout_stride>);
using func_F_F1F1_t =
    void (*)(complex<float>, mdspan<complex<float>, dextents<ptrdiff_t, 1>, layout_stride>, mdspan<complex<float>, dextents<ptrdiff_t, 1>, layout_stride>);
using func_D_D1D1_t =
    void (*)(complex<double>, mdspan<complex<double>, dextents<ptrdiff_t, 1>, layout_stride>, mdspan<complex<double>, dextents<ptrdiff_t, 1>, layout_stride>);

using func_f_f1f1f1_t =
    void (*)(float, mdspan<float, dextents<ptrdiff_t, 1>, layout_stride>, mdspan<float, dextents<ptrdiff_t, 1>, layout_stride>, mdspan<float, dextents<ptrdiff_t, 1>, layout_stride>);
using func_d_d1d1d1_t =
    void (*)(double, mdspan<double, dextents<ptrdiff_t, 1>, layout_stride>, mdspan<double, dextents<ptrdiff_t, 1>, layout_stride>, mdspan<double, dextents<ptrdiff_t, 1>, layout_stride>);
using func_F_F1F1F1_t =
    void (*)(complex<float>, mdspan<complex<float>, dextents<ptrdiff_t, 1>, layout_stride>, mdspan<complex<float>, dextents<ptrdiff_t, 1>, layout_stride>, mdspan<complex<float>, dextents<ptrdiff_t, 1>, layout_stride>);
using func_D_D1D1D1_t =
    void (*)(complex<double>, mdspan<complex<double>, dextents<ptrdiff_t, 1>, layout_stride>, mdspan<complex<double>, dextents<ptrdiff_t, 1>, layout_stride>, mdspan<complex<double>, dextents<ptrdiff_t, 1>, layout_stride>);

using func_f_f2f2_t =
    void (*)(float, mdspan<float, dextents<ptrdiff_t, 2>, layout_stride>, mdspan<float, dextents<ptrdiff_t, 2>, layout_stride>);
using func_d_d2d2_t =
    void (*)(double, mdspan<double, dextents<ptrdiff_t, 2>, layout_stride>, mdspan<double, dextents<ptrdiff_t, 2>, layout_stride>);
using func_F_F2F2_t =
    void (*)(complex<float>, mdspan<complex<float>, dextents<ptrdiff_t, 2>, layout_stride>, mdspan<complex<float>, dextents<ptrdiff_t, 2>, layout_stride>);
using func_D_D2D2_t =
    void (*)(complex<double>, mdspan<complex<double>, dextents<ptrdiff_t, 2>, layout_stride>, mdspan<complex<double>, dextents<ptrdiff_t, 2>, layout_stride>);

using func_bf_f2f2_t =
    void (*)(bool, float, mdspan<float, dextents<ptrdiff_t, 2>, layout_stride>, mdspan<float, dextents<ptrdiff_t, 2>, layout_stride>);
using func_bd_d2d2_t =
    void (*)(bool, double, mdspan<double, dextents<ptrdiff_t, 2>, layout_stride>, mdspan<double, dextents<ptrdiff_t, 2>, layout_stride>);

using func_fb_f2f2_t =
    void (*)(float, bool, mdspan<float, dextents<ptrdiff_t, 2>, layout_stride>, mdspan<float, dextents<ptrdiff_t, 2>, layout_stride>);
using func_db_d2d2_t =
    void (*)(double, bool, mdspan<double, dextents<ptrdiff_t, 2>, layout_stride>, mdspan<double, dextents<ptrdiff_t, 2>, layout_stride>);

using func_Flb_F2F2_t =
    void (*)(complex<float>, long, bool, mdspan<complex<float>, dextents<ptrdiff_t, 2>, layout_stride>, mdspan<complex<float>, dextents<ptrdiff_t, 2>, layout_stride>);
using func_Dlb_D2D2_t =
    void (*)(complex<double>, long, bool, mdspan<complex<double>, dextents<ptrdiff_t, 2>, layout_stride>, mdspan<complex<double>, dextents<ptrdiff_t, 2>, layout_stride>);

using func_ff_F2_t = void (*)(float, float, mdspan<complex<float>, dextents<ptrdiff_t, 2>, layout_stride>);
using func_dd_D2_t = void (*)(double, double, mdspan<complex<double>, dextents<ptrdiff_t, 2>, layout_stride>);

extern const char *lpn_all_doc;
extern const char *lpmn_doc;
extern const char *clpmn_doc;
extern const char *lqn_doc;
extern const char *lqmn_doc;
extern const char *rctj_doc;
extern const char *rcty_doc;
extern const char *sph_harm_all_doc;

// This is needed by sf_error, it is defined in the Cython "_ufuncs_extra_code_common.pxi" for "_generate_pyx.py".
// It exists to "call PyUFunc_getfperr in a context where PyUFunc_API array is initialized", but here we are
// already in such a context.
extern "C" int wrap_PyUFunc_getfperr() { return PyUFunc_getfperr(); }

static PyModuleDef _gufuncs_def = {
    PyModuleDef_HEAD_INIT,
    .m_name = "_gufuncs",
    .m_size = -1,
};

PyMODINIT_FUNC PyInit__gufuncs() {
    if (!SpecFun_Initialize()) {
        return nullptr;
    }

    PyObject *_gufuncs = PyModule_Create(&_gufuncs_def);
    if (_gufuncs == nullptr) {
        return nullptr;
    }

    PyObject *lpn_all = SpecFun_NewGUFunc(
        {static_cast<func_f_f1_t>(::lpn_all), static_cast<func_d_d1_t>(::lpn_all), static_cast<func_F_F1_t>(::lpn_all),
         static_cast<func_D_D1_t>(::lpn_all)},
        1, "lpn_all", lpn_all_doc, "()->(np1)"
    );
    PyModule_AddObjectRef(_gufuncs, "lpn_all", lpn_all);

    PyObject *lpn_all_until_jac = SpecFun_NewGUFunc(
        {static_cast<func_f_f1f1_t>(::lpn_all_until_jac), static_cast<func_d_d1d1_t>(::lpn_all_until_jac),
         static_cast<func_F_F1F1_t>(::lpn_all_until_jac), static_cast<func_D_D1D1_t>(::lpn_all_until_jac)},
        2, "lpn_all_until_jac", lpn_all_doc, "()->(np1),(np1)"
    );
    PyModule_AddObjectRef(_gufuncs, "lpn_all_until_jac", lpn_all_until_jac);

    PyObject *lpn_all_until_hess = SpecFun_NewGUFunc(
        {static_cast<func_f_f1f1f1_t>(::lpn_all_until_hess), static_cast<func_d_d1d1d1_t>(::lpn_all_until_hess),
         static_cast<func_F_F1F1F1_t>(::lpn_all_until_hess), static_cast<func_D_D1D1D1_t>(::lpn_all_until_hess)},
        3, "lpn_all_until_hess", lpn_all_doc, "()->(np1),(np1),(np1)"
    );
    PyModule_AddObjectRef(_gufuncs, "lpn_all_until_hess", lpn_all_until_hess);

    PyObject *_lpmn = SpecFun_NewGUFunc(
        {static_cast<func_bf_f2f2_t>(::lpmn), static_cast<func_bd_d2d2_t>(::lpmn)}, 2, "_lpmn", lpmn_doc,
        "(),()->(mp1,np1),(mp1,np1)"
    );
    PyModule_AddObjectRef(_gufuncs, "_lpmn", _lpmn);

    PyObject *_clpmn = SpecFun_NewGUFunc(
        {static_cast<func_Flb_F2F2_t>(special::clpmn), static_cast<func_Dlb_D2D2_t>(special::clpmn)}, 2, "_clpmn",
        clpmn_doc, "(),(),()->(mp1,np1),(mp1,np1)"
    );
    PyModule_AddObjectRef(_gufuncs, "_clpmn", _clpmn);

    PyObject *_lqn = SpecFun_NewGUFunc(
        {static_cast<func_f_f1f1_t>(special::lqn), static_cast<func_d_d1d1_t>(special::lqn),
         static_cast<func_F_F1F1_t>(special::lqn), static_cast<func_D_D1D1_t>(special::lqn)},
        2, "_lqn", lqn_doc, "()->(np1),(np1)"
    );
    PyModule_AddObjectRef(_gufuncs, "_lqn", _lqn);

    PyObject *_lqmn = SpecFun_NewGUFunc(
        {static_cast<func_f_f2f2_t>(special::lqmn), static_cast<func_d_d2d2_t>(special::lqmn),
         static_cast<func_F_F2F2_t>(special::lqmn), static_cast<func_D_D2D2_t>(special::lqmn)},
        2, "_lqmn", lqmn_doc, "()->(mp1,np1),(mp1,np1)"
    );
    PyModule_AddObjectRef(_gufuncs, "_lqmn", _lqmn);

    PyObject *_rctj = SpecFun_NewGUFunc(
        {static_cast<func_f_f1f1_t>(special::rctj), static_cast<func_d_d1d1_t>(special::rctj)}, 2, "_rctj", rctj_doc,
        "()->(np1),(np1)"
    );
    PyModule_AddObjectRef(_gufuncs, "_rctj", _rctj);

    PyObject *_rcty = SpecFun_NewGUFunc(
        {static_cast<func_f_f1f1_t>(special::rcty), static_cast<func_d_d1d1_t>(special::rcty)}, 2, "_rcty", rcty_doc,
        "()->(np1),(np1)"
    );
    PyModule_AddObjectRef(_gufuncs, "_rcty", _rcty);

    PyObject *_sph_harm_all = SpecFun_NewGUFunc(
        {static_cast<func_dd_D2_t>(special::sph_harm_all), static_cast<func_ff_F2_t>(special::sph_harm_all)}, 1,
        "_sph_harm_all", sph_harm_all_doc, "(),()->(mp1,np1)"
    );
    PyModule_AddObjectRef(_gufuncs, "_sph_harm_all", _sph_harm_all);

    return _gufuncs;
}
