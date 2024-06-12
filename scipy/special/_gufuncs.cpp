#include "ufunc.h"

#include "special.h"

using namespace std;

using cfloat = complex<float>;
using cdouble = complex<double>;

using float_1d = mdspan<float, dextents<ptrdiff_t, 1>, layout_stride>;
using float_2d = mdspan<float, dextents<ptrdiff_t, 2>, layout_stride>;
using double_1d = mdspan<double, dextents<ptrdiff_t, 1>, layout_stride>;
using double_2d = mdspan<double, dextents<ptrdiff_t, 2>, layout_stride>;
using cfloat_1d = mdspan<cfloat, dextents<ptrdiff_t, 1>, layout_stride>;
using cfloat_2d = mdspan<cfloat, dextents<ptrdiff_t, 2>, layout_stride>;
using cdouble_1d = mdspan<cdouble, dextents<ptrdiff_t, 1>, layout_stride>;
using cdouble_2d = mdspan<cdouble, dextents<ptrdiff_t, 2>, layout_stride>;

// 1 input, 1 output
using func_f_f1_t = void (*)(float, float_1d);
using func_f_f2_t = void (*)(float, float_2d);
using func_d_d1_t = void (*)(double, double_1d);
using func_d_d2_t = void (*)(double, double_2d);
using func_F_F1_t = void (*)(cfloat, cfloat_1d);
using func_D_D1_t = void (*)(cdouble, cdouble_1d);

// 1 input, 2 outputs
using func_f_f1f1_t = void (*)(float, float_1d, float_1d);
using func_f_f2f2_t = void (*)(float, float_2d, float_2d);
using func_d_d1d1_t = void (*)(double, double_1d, double_1d);
using func_d_d2d2_t = void (*)(double, double_2d, double_2d);
using func_F_F1F1_t = void (*)(cfloat, cfloat_1d, cfloat_1d);
using func_F_F2F2_t = void (*)(cfloat, cfloat_2d, cfloat_2d);
using func_D_D1D1_t = void (*)(cdouble, cdouble_1d, cdouble_1d);
using func_D_D2D2_t = void (*)(cdouble, cdouble_2d, cdouble_2d);

// 1 input, 3 outputs
using func_f_f1f1f1_t = void (*)(float, float_1d, float_1d, float_1d);
using func_f_f2f2f2_t = void (*)(float, float_2d, float_2d, float_2d);
using func_d_d1d1d1_t = void (*)(double, double_1d, double_1d, double_1d);
using func_d_d2d2d2_t = void (*)(double, double_2d, double_2d, double_2d);
using func_F_F1F1F1_t = void (*)(cfloat, cfloat_1d, cfloat_1d, cfloat_1d);
using func_D_D1D1D1_t = void (*)(cdouble, cdouble_1d, cdouble_1d, cdouble_1d);

// 2 inputs, 1 output
using func_ff_F2_t = void (*)(float, float, cfloat_2d);
using func_dd_D2_t = void (*)(double, double, cdouble_2d);
using func_qF_F2_t = void (*)(long long int, cfloat, cfloat_2d);
using func_qD_D2_t = void (*)(long long int, cdouble, cdouble_2d);

// 2 inputs, 2 outputs
using func_qF_F2F2_t = void (*)(long long int, cfloat, cfloat_2d, cfloat_2d);
using func_qD_D2D2_t = void (*)(long long int, cdouble, cdouble_2d, cdouble_2d);

// 2 inputs, 3 outputs
using func_qF_F2F2F2_t = void (*)(long long int, cfloat, cfloat_2d, cfloat_2d, cfloat_2d);
using func_qD_D2D2D2_t = void (*)(long long int, cdouble, cdouble_2d, cdouble_2d, cdouble_2d);

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
    import_array();
    import_umath();
    if (PyErr_Occurred()) {
        return NULL;
    }

    PyObject *_gufuncs = PyModule_Create(&_gufuncs_def);
    if (_gufuncs == nullptr) {
        return NULL;
    }

    PyObject *legendre_p_all = Py_BuildValue(
        "(N,N,N)",
        SpecFun_NewGUFunc(
            {static_cast<func_d_d1_t>(::legendre_p_all), static_cast<func_f_f1_t>(::legendre_p_all),
             static_cast<func_D_D1_t>(::legendre_p_all), static_cast<func_F_F1_t>(::legendre_p_all)},
            1, "legendre_p_all_diff_0", nullptr, "()->(np1)"
        ),
        SpecFun_NewGUFunc(
            {static_cast<func_d_d1d1_t>(::legendre_p_all), static_cast<func_f_f1f1_t>(::legendre_p_all),
             static_cast<func_D_D1D1_t>(::legendre_p_all), static_cast<func_F_F1F1_t>(::legendre_p_all)},
            2, "legendre_p_all_diff_1", nullptr, "()->(np1),(np1)"
        ),
        SpecFun_NewGUFunc(
            {static_cast<func_d_d1d1d1_t>(::legendre_p_all), static_cast<func_f_f1f1f1_t>(::legendre_p_all),
             static_cast<func_D_D1D1D1_t>(::legendre_p_all), static_cast<func_F_F1F1F1_t>(::legendre_p_all)},
            3, "legendre_p_all_diff_2", nullptr, "()->(np1),(np1),(np1)"
        )
    );
    PyModule_AddObjectRef(_gufuncs, "legendre_p_all", legendre_p_all);

    PyObject *assoc_legendre_p_all = Py_BuildValue(
        "{O:(N,N,N),O:(N,N,N)}", Py_True,
        SpecFun_NewGUFunc(
            {[](double z, double_2d res) { ::assoc_legendre_p_all(assoc_legendre_norm, z, res); },
             [](float z, float_2d res) { ::assoc_legendre_p_all(assoc_legendre_norm, z, res); }},
            1, "assoc_legendre_p_all_norm_diff_0", nullptr, "()->(np1,mpmp1)"
        ),
        SpecFun_NewGUFunc(
            {[](double z, double_2d res, double_2d res_jac) {
                 ::assoc_legendre_p_all(assoc_legendre_norm, z, res, res_jac);
             },
             [](float z, float_2d res, float_2d res_jac) {
                 ::assoc_legendre_p_all(assoc_legendre_norm, z, res, res_jac);
             }},
            2, "assoc_legendre_p_all_norm_diff_1", nullptr, "()->(np1,mpmp1),(np1,mpmp1)"
        ),
        SpecFun_NewGUFunc(
            {[](double z, double_2d res, double_2d res_jac, double_2d res_hess) {
                 ::assoc_legendre_p_all(assoc_legendre_norm, z, res, res_jac, res_hess);
             },
             [](float z, float_2d res, float_2d res_jac, float_2d res_hess) {
                 ::assoc_legendre_p_all(assoc_legendre_norm, z, res, res_jac, res_hess);
             }},
            3, "assoc_legendre_p_all_norm_diff_2", nullptr, "()->(np1,mpmp1),(np1,mpmp1),(np1,mpmp1)"
        ),
        Py_False,
        SpecFun_NewGUFunc(
            {[](double z, double_2d res) { ::assoc_legendre_p_all(assoc_legendre_unnorm, z, res); },
             [](float z, float_2d res) { ::assoc_legendre_p_all(assoc_legendre_unnorm, z, res); }},
            1, "assoc_legendre_p_all_unnorm_diff_0", nullptr, "()->(np1,mpmp1)"
        ),
        SpecFun_NewGUFunc(
            {[](double z, double_2d res, double_2d res_jac) {
                 ::assoc_legendre_p_all(assoc_legendre_unnorm, z, res, res_jac);
             },
             [](float z, float_2d res, float_2d res_jac) {
                 ::assoc_legendre_p_all(assoc_legendre_unnorm, z, res, res_jac);
             }},
            2, "assoc_legendre_p_all_unnorm_diff_1", nullptr, "()->(np1,mpmp1),(np1,mpmp1)"
        ),
        SpecFun_NewGUFunc(
            {[](double z, double_2d res, double_2d res_jac, double_2d res_hess) {
                 ::assoc_legendre_p_all(assoc_legendre_unnorm, z, res, res_jac, res_hess);
             },
             [](float z, float_2d res, float_2d res_jac, float_2d res_hess) {
                 ::assoc_legendre_p_all(assoc_legendre_unnorm, z, res, res_jac, res_hess);
             }},
            3, "assoc_legendre_p_all_unnorm_diff_1", nullptr, "()->(np1,mpmp1),(np1,mpmp1),(np1,mpmp1)"
        )
    );
    PyModule_AddObjectRef(_gufuncs, "assoc_legendre_p_all", assoc_legendre_p_all);

    PyObject *multi_assoc_legendre_p_all = Py_BuildValue(
        "{O:(N,N,N),O:(N,N,N)}", Py_True,
        SpecFun_NewGUFunc(
            {[](long long int type, cdouble z, cdouble_2d res) {
                 ::multi_assoc_legendre_p_all(assoc_legendre_norm, type, z, res);
             },
             [](long long int type, cfloat z, cfloat_2d res) {
                 ::multi_assoc_legendre_p_all(assoc_legendre_norm, type, z, res);
             }},
            1, "multi_assoc_legendre_p_all_norm_diff_0", nullptr, "(),()->(np1,mpmp1)"
        ),
        SpecFun_NewGUFunc(
            {[](long long int type, cdouble z, cdouble_2d res, cdouble_2d res_jac) {
                 ::multi_assoc_legendre_p_all(assoc_legendre_norm, type, z, res, res_jac);
             },
             [](long long int type, cfloat z, cfloat_2d res, cfloat_2d res_jac) {
                 ::multi_assoc_legendre_p_all(assoc_legendre_norm, type, z, res, res_jac);
             }},
            2, "multi_assoc_legendre_p_all_norm_diff_1", nullptr, "(),()->(np1,mpmp1),(np1,mpmp1)"
        ),
        SpecFun_NewGUFunc(
            {[](long long int type, cdouble z, cdouble_2d res, cdouble_2d res_jac, cdouble_2d res_hess) {
                 ::multi_assoc_legendre_p_all(assoc_legendre_norm, type, z, res, res_jac, res_hess);
             },
             [](long long int type, cfloat z, cfloat_2d res, cfloat_2d res_jac, cfloat_2d res_hess) {
                 ::multi_assoc_legendre_p_all(assoc_legendre_norm, type, z, res, res_jac, res_hess);
             }},
            3, "multi_assoc_legendre_p_all_norm_diff_2", nullptr, "(),()->(np1,mpmp1),(np1,mpmp1),(np1,mpmp1)"
        ),
        Py_False,
        SpecFun_NewGUFunc(
            {[](long long int type, cdouble z, cdouble_2d res) {
                 ::multi_assoc_legendre_p_all(assoc_legendre_unnorm, type, z, res);
             },
             [](long long int type, cfloat z, cfloat_2d res) {
                 ::multi_assoc_legendre_p_all(assoc_legendre_unnorm, type, z, res);
             }},
            1, "multi_assoc_legendre_p_all_unnorm_diff_0", nullptr, "(),()->(np1,mpmp1)"
        ),
        SpecFun_NewGUFunc(
            {[](long long int type, cdouble z, cdouble_2d res, cdouble_2d res_jac) {
                 ::multi_assoc_legendre_p_all(assoc_legendre_unnorm, type, z, res, res_jac);
             },
             [](long long int type, cfloat z, cfloat_2d res, cfloat_2d res_jac) {
                 ::multi_assoc_legendre_p_all(assoc_legendre_unnorm, type, z, res, res_jac);
             }},
            2, "multi_assoc_legendre_p_all_unnorm_diff_1", nullptr, "(),()->(np1,mpmp1),(np1,mpmp1)"
        ),
        SpecFun_NewGUFunc(
            {[](long long int type, cdouble z, cdouble_2d res, cdouble_2d res_jac, cdouble_2d res_hess) {
                 ::multi_assoc_legendre_p_all(assoc_legendre_unnorm, type, z, res, res_jac, res_hess);
             },
             [](long long int type, cfloat z, cfloat_2d res, cfloat_2d res_jac, cfloat_2d res_hess) {
                 ::multi_assoc_legendre_p_all(assoc_legendre_unnorm, type, z, res, res_jac, res_hess);
             }},
            3, "multi_assoc_legendre_p_all_unnorm_diff_2", nullptr, "(),()->(np1,mpmp1),(np1,mpmp1),(np1,mpmp1)"
        )
    );
    PyModule_AddObjectRef(_gufuncs, "multi_assoc_legendre_p_all", multi_assoc_legendre_p_all);

    PyObject *sph_legendre_p_all = Py_BuildValue(
        "(N,N,N)",
        SpecFun_NewGUFunc(
            {static_cast<func_d_d2_t>(::sph_legendre_p_all), static_cast<func_f_f2_t>(::sph_legendre_p_all)}, 1,
            "sph_legendre_p_all_diff_0", nullptr, "()->(np1,mpmp1)"
        ),
        SpecFun_NewGUFunc(
            {static_cast<func_d_d2d2_t>(::sph_legendre_p_all), static_cast<func_f_f2f2_t>(::sph_legendre_p_all)}, 2,
            "sph_legendre_p_all_diff_1", nullptr, "()->(np1,mpmp1),(np1,mpmp1)"
        ),
        SpecFun_NewGUFunc(
            {static_cast<func_d_d2d2d2_t>(::sph_legendre_p_all), static_cast<func_f_f2f2f2_t>(::sph_legendre_p_all)}, 3,
            "sph_legendre_p_all_diff_2", nullptr, "()->(np1,mpmp1),(np1,mpmp1),(np1,mpmp1)"
        )
    );
    PyModule_AddObjectRef(_gufuncs, "sph_legendre_p_all", sph_legendre_p_all);

    PyObject *_lqn = SpecFun_NewGUFunc(
        {static_cast<func_d_d1d1_t>(special::lqn), static_cast<func_f_f1f1_t>(special::lqn),
         static_cast<func_D_D1D1_t>(special::lqn), static_cast<func_F_F1F1_t>(special::lqn)},
        2, "_lqn", lqn_doc, "()->(np1),(np1)"
    );
    PyModule_AddObjectRef(_gufuncs, "_lqn", _lqn);

    PyObject *_lqmn = SpecFun_NewGUFunc(
        {static_cast<func_d_d2d2_t>(special::lqmn), static_cast<func_f_f2f2_t>(special::lqmn),
         static_cast<func_D_D2D2_t>(special::lqmn), static_cast<func_F_F2F2_t>(special::lqmn)},
        2, "_lqmn", lqmn_doc, "()->(mp1,np1),(mp1,np1)"
    );
    PyModule_AddObjectRef(_gufuncs, "_lqmn", _lqmn);

    PyObject *sph_harm_y_all = SpecFun_NewGUFunc(
        {static_cast<func_dd_D2_t>(::sph_harm_y_all), static_cast<func_ff_F2_t>(::sph_harm_y_all)}, 1, "sph_harm_y_all",
        sph_harm_all_doc, "(),()->(mpmp1,np1)"
    );
    PyModule_AddObjectRef(_gufuncs, "sph_harm_y_all", sph_harm_y_all);

    PyObject *_rctj = SpecFun_NewGUFunc(
        {static_cast<func_d_d1d1_t>(special::rctj), static_cast<func_f_f1f1_t>(special::rctj)}, 2, "_rctj", rctj_doc,
        "()->(np1),(np1)"
    );
    PyModule_AddObjectRef(_gufuncs, "_rctj", _rctj);

    PyObject *_rcty = SpecFun_NewGUFunc(
        {static_cast<func_d_d1d1_t>(special::rcty), static_cast<func_f_f1f1_t>(special::rcty)}, 2, "_rcty", rcty_doc,
        "()->(np1),(np1)"
    );
    PyModule_AddObjectRef(_gufuncs, "_rcty", _rcty);

    return _gufuncs;
}
