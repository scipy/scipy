#include "xsf/numpy.h"

#include "xsf_special.h"

using xsf::numpy::cdouble;
using xsf::numpy::cdouble_1d;
using xsf::numpy::cdouble_2d;
using xsf::numpy::cdouble_3d;
using xsf::numpy::cdouble_4d;
using xsf::numpy::cfloat;
using xsf::numpy::cfloat_1d;
using xsf::numpy::cfloat_2d;
using xsf::numpy::cfloat_3d;
using xsf::numpy::cfloat_4d;
using xsf::numpy::double_1d;
using xsf::numpy::double_2d;
using xsf::numpy::double_3d;
using xsf::numpy::double_4d;
using xsf::numpy::float_1d;
using xsf::numpy::float_2d;
using xsf::numpy::float_3d;
using xsf::numpy::float_4d;

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

static PyModuleDef _gufuncs_def = {PyModuleDef_HEAD_INIT, "_gufuncs", NULL, -1, NULL, NULL, NULL, NULL, NULL};

template <size_t NOut>
void legendre_map_dims(const npy_intp *dims, npy_intp *new_dims) {
    for (size_t i = 0; i < NOut; ++i) {
        new_dims[i] = dims[0];
    }
}

template <size_t NOut>
void assoc_legendre_map_dims(const npy_intp *dims, npy_intp *new_dims) {
    for (size_t i = 0; i < NOut; ++i) {
        new_dims[2 * i] = dims[0];
        new_dims[2 * i + 1] = dims[1];
    }
}

template <size_t NOut>
void sph_harm_map_dims(const npy_intp *dims, npy_intp *new_dims);

template <>
void sph_harm_map_dims<1>(const npy_intp *dims, npy_intp *new_dims) {
    new_dims[0] = dims[0];
    new_dims[1] = dims[1];
}

template <>
void sph_harm_map_dims<2>(const npy_intp *dims, npy_intp *new_dims) {
    new_dims[0] = dims[0];
    new_dims[1] = dims[1];

    new_dims[2] = 2;
    new_dims[3] = dims[0];
    new_dims[4] = dims[1];
}

template <>
void sph_harm_map_dims<3>(const npy_intp *dims, npy_intp *new_dims) {
    new_dims[0] = dims[0];
    new_dims[1] = dims[1];

    new_dims[2] = 2;
    new_dims[3] = dims[0];
    new_dims[4] = dims[1];

    new_dims[5] = 2;
    new_dims[6] = 2;
    new_dims[7] = dims[0];
    new_dims[8] = dims[1];
}

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

#if Py_GIL_DISABLED
    PyUnstable_Module_SetGIL(_gufuncs, Py_MOD_GIL_NOT_USED);
#endif

    PyObject *legendre_p_all = Py_BuildValue(
        "(N, N, N)",
        xsf::numpy::gufunc(
            {static_cast<xsf::numpy::d_d1>(::legendre_p_all), static_cast<xsf::numpy::f_f1>(::legendre_p_all),
             static_cast<xsf::numpy::D_D1>(::legendre_p_all), static_cast<xsf::numpy::F_F1>(::legendre_p_all)},
            1, "legendre_p_all", nullptr, "()->(np1)", legendre_map_dims<1>),
        xsf::numpy::gufunc(
            {static_cast<xsf::numpy::d_d1d1>(::legendre_p_all), static_cast<xsf::numpy::f_f1f1>(::legendre_p_all),
             static_cast<xsf::numpy::D_D1D1>(::legendre_p_all), static_cast<xsf::numpy::F_F1F1>(::legendre_p_all)},
            2, "legendre_p_all", nullptr, "()->(np1),(np1)", legendre_map_dims<2>),
        xsf::numpy::gufunc(
            {static_cast<xsf::numpy::d_d1d1d1>(::legendre_p_all), static_cast<xsf::numpy::f_f1f1f1>(::legendre_p_all),
             static_cast<xsf::numpy::D_D1D1D1>(::legendre_p_all), static_cast<xsf::numpy::F_F1F1F1>(::legendre_p_all)},
            3, "legendre_p_all", nullptr, "()->(np1),(np1),(np1)", legendre_map_dims<3>));
    PyModule_AddObjectRef(_gufuncs, "legendre_p_all", legendre_p_all);

    // key is norm, diff_n
    PyObject *assoc_legendre_p_all = Py_BuildValue(
        "{(O, i): N, (O, i): N, (O, i): N, (O, i): N, (O, i): N, (O, i): N}", Py_True, 0,
        xsf::numpy::gufunc({[](double z, long long int branch_cut, double_2d res) {
                                ::assoc_legendre_p_all(assoc_legendre_norm, z, branch_cut, res);
                            },
                            [](float z, long long int branch_cut, float_2d res) {
                                ::assoc_legendre_p_all(assoc_legendre_norm, z, branch_cut, res);
                            },
                            [](cdouble z, long long int branch_cut, cdouble_2d res) {
                                ::assoc_legendre_p_all(assoc_legendre_norm, z, branch_cut, res);
                            },
                            [](cfloat z, long long int branch_cut, cfloat_2d res) {
                                ::assoc_legendre_p_all(assoc_legendre_norm, z, branch_cut, res);
                            }},
                           1, "assoc_legendre_p_all", nullptr, "(),()->(np1,mpmp1)", assoc_legendre_map_dims<1>),
        Py_True, 1,
        xsf::numpy::gufunc({[](double z, long long int branch_cut, double_2d res, double_2d res_jac) {
                                ::assoc_legendre_p_all(assoc_legendre_norm, z, branch_cut, res, res_jac);
                            },
                            [](float z, long long int branch_cut, float_2d res, float_2d res_jac) {
                                ::assoc_legendre_p_all(assoc_legendre_norm, z, branch_cut, res, res_jac);
                            },
                            [](cdouble z, long long int branch_cut, cdouble_2d res, cdouble_2d res_jac) {
                                ::assoc_legendre_p_all(assoc_legendre_norm, z, branch_cut, res, res_jac);
                            },
                            [](cfloat z, long long int branch_cut, cfloat_2d res, cfloat_2d res_jac) {
                                ::assoc_legendre_p_all(assoc_legendre_norm, z, branch_cut, res, res_jac);
                            }},
                           2, "assoc_legendre_p_all", nullptr, "(),()->(np1,mpmp1),(np1,mpmp1)",
                           assoc_legendre_map_dims<2>),
        Py_True, 2,
        xsf::numpy::gufunc(
            {[](double z, long long int branch_cut, double_2d res, double_2d res_jac, double_2d res_hess) {
                 ::assoc_legendre_p_all(assoc_legendre_norm, z, branch_cut, res, res_jac, res_hess);
             },
             [](float z, long long int branch_cut, float_2d res, float_2d res_jac, float_2d res_hess) {
                 ::assoc_legendre_p_all(assoc_legendre_norm, z, branch_cut, res, res_jac, res_hess);
             },
             [](cdouble z, long long int branch_cut, cdouble_2d res, cdouble_2d res_jac, cdouble_2d res_hess) {
                 ::assoc_legendre_p_all(assoc_legendre_norm, z, branch_cut, res, res_jac, res_hess);
             },
             [](cfloat z, long long int branch_cut, cfloat_2d res, cfloat_2d res_jac, cfloat_2d res_hess) {
                 ::assoc_legendre_p_all(assoc_legendre_norm, z, branch_cut, res, res_jac, res_hess);
             }},
            3, "assoc_legendre_p_all", nullptr, "(),()->(np1,mpmp1),(np1,mpmp1),(np1,mpmp1)",
            assoc_legendre_map_dims<3>),
        Py_False, 0,
        xsf::numpy::gufunc({[](double z, long long int branch_cut, double_2d res) {
                                ::assoc_legendre_p_all(assoc_legendre_unnorm, z, branch_cut, res);
                            },
                            [](float z, long long int branch_cut, float_2d res) {
                                ::assoc_legendre_p_all(assoc_legendre_unnorm, z, branch_cut, res);
                            },
                            [](cdouble z, long long int branch_cut, cdouble_2d res) {
                                ::assoc_legendre_p_all(assoc_legendre_unnorm, z, branch_cut, res);
                            },
                            [](cfloat z, long long int branch_cut, cfloat_2d res) {
                                ::assoc_legendre_p_all(assoc_legendre_unnorm, z, branch_cut, res);
                            }},
                           1, "assoc_legendre_p_all", nullptr, "(),()->(np1,mpmp1)", assoc_legendre_map_dims<1>),
        Py_False, 1,
        xsf::numpy::gufunc({[](double z, long long int branch_cut, double_2d res, double_2d res_jac) {
                                ::assoc_legendre_p_all(assoc_legendre_unnorm, z, branch_cut, res, res_jac);
                            },
                            [](float z, long long int branch_cut, float_2d res, float_2d res_jac) {
                                ::assoc_legendre_p_all(assoc_legendre_unnorm, z, branch_cut, res, res_jac);
                            },
                            [](cdouble z, long long int branch_cut, cdouble_2d res, cdouble_2d res_jac) {
                                ::assoc_legendre_p_all(assoc_legendre_unnorm, z, branch_cut, res, res_jac);
                            },
                            [](cfloat z, long long int branch_cut, cfloat_2d res, cfloat_2d res_jac) {
                                ::assoc_legendre_p_all(assoc_legendre_unnorm, z, branch_cut, res, res_jac);
                            }},
                           2, "assoc_legendre_p_all", nullptr, "(),()->(np1,mpmp1),(np1,mpmp1)",
                           assoc_legendre_map_dims<2>),
        Py_False, 2,
        xsf::numpy::gufunc(
            {[](double z, long long int branch_cut, double_2d res, double_2d res_jac, double_2d res_hess) {
                 ::assoc_legendre_p_all(assoc_legendre_unnorm, z, branch_cut, res, res_jac, res_hess);
             },
             [](float z, long long int branch_cut, float_2d res, float_2d res_jac, float_2d res_hess) {
                 ::assoc_legendre_p_all(assoc_legendre_unnorm, z, branch_cut, res, res_jac, res_hess);
             },
             [](cdouble z, long long int branch_cut, cdouble_2d res, cdouble_2d res_jac, cdouble_2d res_hess) {
                 ::assoc_legendre_p_all(assoc_legendre_unnorm, z, branch_cut, res, res_jac, res_hess);
             },
             [](cfloat z, long long int branch_cut, cfloat_2d res, cfloat_2d res_jac, cfloat_2d res_hess) {
                 ::assoc_legendre_p_all(assoc_legendre_unnorm, z, branch_cut, res, res_jac, res_hess);
             }},
            3, "assoc_legendre_p_all", nullptr, "(),()->(np1,mpmp1),(np1,mpmp1),(np1,mpmp1)",
            assoc_legendre_map_dims<3>));
    PyModule_AddObjectRef(_gufuncs, "assoc_legendre_p_all", assoc_legendre_p_all);

    PyObject *sph_legendre_p_all = Py_BuildValue(
        "(N, N, N)",
        xsf::numpy::gufunc(
            {static_cast<xsf::numpy::d_d2>(::sph_legendre_p_all), static_cast<xsf::numpy::f_f2>(::sph_legendre_p_all)},
            1, "sph_legendre_p_all", nullptr, "()->(np1,mpmp1)", assoc_legendre_map_dims<1>),
        xsf::numpy::gufunc({static_cast<xsf::numpy::d_d2d2>(::sph_legendre_p_all),
                            static_cast<xsf::numpy::f_f2f2>(::sph_legendre_p_all)},
                           2, "sph_legendre_p_all", nullptr, "()->(np1,mpmp1),(np1,mpmp1)", assoc_legendre_map_dims<2>),
        xsf::numpy::gufunc({static_cast<xsf::numpy::d_d2d2d2>(::sph_legendre_p_all),
                            static_cast<xsf::numpy::f_f2f2f2>(::sph_legendre_p_all)},
                           3, "sph_legendre_p_all", nullptr, "()->(np1,mpmp1),(np1,mpmp1),(np1,mpmp1)",
                           assoc_legendre_map_dims<3>));
    PyModule_AddObjectRef(_gufuncs, "sph_legendre_p_all", sph_legendre_p_all);

    PyObject *_lqn =
        xsf::numpy::gufunc({static_cast<xsf::numpy::d_d1d1>(xsf::lqn), static_cast<xsf::numpy::f_f1f1>(xsf::lqn),
                            static_cast<xsf::numpy::D_D1D1>(xsf::lqn), static_cast<xsf::numpy::F_F1F1>(xsf::lqn)},
                           2, "_lqn", lqn_doc, "()->(np1),(np1)", legendre_map_dims<2>);
    PyModule_AddObjectRef(_gufuncs, "_lqn", _lqn);

    PyObject *_lqmn =
        xsf::numpy::gufunc({static_cast<xsf::numpy::d_d2d2>(xsf::lqmn), static_cast<xsf::numpy::f_f2f2>(xsf::lqmn),
                            static_cast<xsf::numpy::D_D2D2>(xsf::lqmn), static_cast<xsf::numpy::F_F2F2>(xsf::lqmn)},
                           2, "_lqmn", lqmn_doc, "()->(mp1,np1),(mp1,np1)", assoc_legendre_map_dims<2>);
    PyModule_AddObjectRef(_gufuncs, "_lqmn", _lqmn);

    PyObject *sph_harm_y_all = Py_BuildValue(
        "(N,N,N)",
        xsf::numpy::gufunc(
            {static_cast<xsf::numpy::dd_D2>(::sph_harm_y_all), static_cast<xsf::numpy::ff_F2>(::sph_harm_y_all)}, 1,
            "sph_harm_y_all", nullptr, "(),()->(np1,mpmp1)", sph_harm_map_dims<1>),
        xsf::numpy::gufunc(
            {static_cast<xsf::numpy::dd_D2D3>(::sph_harm_y_all), static_cast<xsf::numpy::ff_F2F3>(::sph_harm_y_all)}, 2,
            "sph_harm_y_all", nullptr, "(),()->(np1,mpmp1),(2,np1,mpmp1)", sph_harm_map_dims<2>),
        xsf::numpy::gufunc({static_cast<xsf::numpy::dd_D2D3D4>(::sph_harm_y_all),
                            static_cast<xsf::numpy::ff_F2F3F4>(::sph_harm_y_all)},
                           3, "sph_harm_y_all", nullptr, "(),()->(np1,mpmp1),(2,np1,mpmp1),(2,2,np1,mpmp1)",
                           sph_harm_map_dims<3>));
    PyModule_AddObjectRef(_gufuncs, "sph_harm_y_all", sph_harm_y_all);

    PyObject *_rctj =
        xsf::numpy::gufunc({static_cast<xsf::numpy::d_d1d1>(xsf::rctj), static_cast<xsf::numpy::f_f1f1>(xsf::rctj)}, 2,
                           "_rctj", rctj_doc, "()->(np1),(np1)", legendre_map_dims<2>);
    PyModule_AddObjectRef(_gufuncs, "_rctj", _rctj);

    PyObject *_rcty =
        xsf::numpy::gufunc({static_cast<xsf::numpy::d_d1d1>(xsf::rcty), static_cast<xsf::numpy::f_f1f1>(xsf::rcty)}, 2,
                           "_rcty", rcty_doc, "()->(np1),(np1)", legendre_map_dims<2>);
    PyModule_AddObjectRef(_gufuncs, "_rcty", _rcty);

    return _gufuncs;
}
