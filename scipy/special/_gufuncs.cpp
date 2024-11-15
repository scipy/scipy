#include "xsf/numpy.h"

#include "xsf_special.h"

extern const char *lpn_all_doc;
extern const char *lpmn_doc;
extern const char *clpmn_doc;
extern const char *lqn_doc;
extern const char *lqmn_doc;
extern const char *rctj_doc;
extern const char *rcty_doc;
extern const char *sph_harm_all_doc;

static PyModuleDef _gufuncs_def = {PyModuleDef_HEAD_INIT, "_gufuncs", NULL, -1, NULL, NULL, NULL, NULL, NULL};

template <size_t NOut>
void legendre_map_dims(const npy_intp *dims, npy_intp *new_dims) {
    for (size_t i = 0; i < NOut; ++i) {
        new_dims[i] = dims[0];
    }
}

void new_legendre_map_dims(const npy_intp *dims, npy_intp *new_dims) { new_dims[0] = dims[0]; }

template <size_t NOut>
void assoc_legendre_map_dims(const npy_intp *dims, npy_intp *new_dims) {
    for (size_t i = 0; i < NOut; ++i) {
        new_dims[2 * i] = dims[0];
        new_dims[2 * i + 1] = dims[1];
    }
}

void sph_harm_map_dims(const npy_intp *dims, npy_intp *new_dims) {
    new_dims[0] = dims[0];
    new_dims[1] = dims[1];
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

    PyObject *legendre_p_all =
        Py_BuildValue("(N, N, N)",
                      xsf::numpy::gufunc({xsf::numpy::compose{xsf::numpy::autodiff()},
                                          static_cast<xsf::numpy::autodiff0_d_d1>(xsf::legendre_p_all),
                                          static_cast<xsf::numpy::autodiff0_f_f1>(xsf::legendre_p_all),
                                          static_cast<xsf::numpy::autodiff0_D_D1>(xsf::legendre_p_all),
                                          static_cast<xsf::numpy::autodiff0_F_F1>(xsf::legendre_p_all)},
                                         1, "legendre_p_all", nullptr, "()->(np1,1)", new_legendre_map_dims),
                      xsf::numpy::gufunc({xsf::numpy::compose{xsf::numpy::autodiff()},
                                          static_cast<xsf::numpy::autodiff1_d_d1>(xsf::legendre_p_all),
                                          static_cast<xsf::numpy::autodiff1_f_f1>(xsf::legendre_p_all),
                                          static_cast<xsf::numpy::autodiff1_D_D1>(xsf::legendre_p_all),
                                          static_cast<xsf::numpy::autodiff1_F_F1>(xsf::legendre_p_all)},
                                         1, "legendre_p_all", nullptr, "()->(np1,2)", new_legendre_map_dims),
                      xsf::numpy::gufunc({xsf::numpy::compose{xsf::numpy::autodiff()},
                                          static_cast<xsf::numpy::autodiff2_d_d1>(xsf::legendre_p_all),
                                          static_cast<xsf::numpy::autodiff2_f_f1>(xsf::legendre_p_all),
                                          static_cast<xsf::numpy::autodiff2_D_D1>(xsf::legendre_p_all),
                                          static_cast<xsf::numpy::autodiff2_F_F1>(xsf::legendre_p_all)},
                                         1, "legendre_p_all", nullptr, "()->(np1,3)", new_legendre_map_dims));
    PyModule_AddObjectRef(_gufuncs, "legendre_p_all", legendre_p_all);

    // key is norm, diff_n
    PyObject *assoc_legendre_p_all = Py_BuildValue(
        "{(O, i): N, (O, i): N, (O, i): N, (O, i): N, (O, i): N, (O, i): N}", Py_True, 0,
        xsf::numpy::gufunc({xsf::numpy::compose{xsf::numpy::autodiff(), xsf::numpy::use_long_long_int()},
                            [](xsf::numpy::autodiff0_double z, int branch_cut, xsf::numpy::autodiff0_double_2d res) {
                                xsf::assoc_legendre_p_all(xsf::assoc_legendre_norm, z, branch_cut, res);
                            },
                            [](xsf::numpy::autodiff0_float z, int branch_cut, xsf::numpy::autodiff0_float_2d res) {
                                xsf::assoc_legendre_p_all(xsf::assoc_legendre_norm, z, branch_cut, res);
                            },
                            [](xsf::numpy::autodiff0_cdouble z, int branch_cut, xsf::numpy::autodiff0_cdouble_2d res) {
                                xsf::assoc_legendre_p_all(xsf::assoc_legendre_norm, z, branch_cut, res);
                            },
                            [](xsf::numpy::autodiff0_cfloat z, int branch_cut, xsf::numpy::autodiff0_cfloat_2d res) {
                                xsf::assoc_legendre_p_all(xsf::assoc_legendre_norm, z, branch_cut, res);
                            }},
                           1, "assoc_legendre_p_all", nullptr, "(),()->(np1,mpmp1,1)", assoc_legendre_map_dims<1>),
        Py_True, 1,
        xsf::numpy::gufunc({xsf::numpy::compose{xsf::numpy::autodiff(), xsf::numpy::use_long_long_int()},
                            [](xsf::numpy::autodiff1_double z, int branch_cut, xsf::numpy::autodiff1_double_2d res) {
                                xsf::assoc_legendre_p_all(xsf::assoc_legendre_norm, z, branch_cut, res);
                            },
                            [](xsf::numpy::autodiff1_float z, int branch_cut, xsf::numpy::autodiff1_float_2d res) {
                                xsf::assoc_legendre_p_all(xsf::assoc_legendre_norm, z, branch_cut, res);
                            },
                            [](xsf::numpy::autodiff1_cdouble z, int branch_cut, xsf::numpy::autodiff1_cdouble_2d res) {
                                xsf::assoc_legendre_p_all(xsf::assoc_legendre_norm, z, branch_cut, res);
                            },
                            [](xsf::numpy::autodiff1_cfloat z, int branch_cut, xsf::numpy::autodiff1_cfloat_2d res) {
                                xsf::assoc_legendre_p_all(xsf::assoc_legendre_norm, z, branch_cut, res);
                            }},
                           1, "assoc_legendre_p_all", nullptr, "(),()->(np1,mpmp1,2)", assoc_legendre_map_dims<1>),
        Py_True, 2,
        xsf::numpy::gufunc({xsf::numpy::compose{xsf::numpy::autodiff(), xsf::numpy::use_long_long_int()},
                            [](xsf::numpy::autodiff2_double z, int branch_cut, xsf::numpy::autodiff2_double_2d res) {
                                xsf::assoc_legendre_p_all(xsf::assoc_legendre_norm, z, branch_cut, res);
                            },
                            [](xsf::numpy::autodiff2_float z, int branch_cut, xsf::numpy::autodiff2_float_2d res) {
                                xsf::assoc_legendre_p_all(xsf::assoc_legendre_norm, z, branch_cut, res);
                            },
                            [](xsf::numpy::autodiff2_cdouble z, int branch_cut, xsf::numpy::autodiff2_cdouble_2d res) {
                                xsf::assoc_legendre_p_all(xsf::assoc_legendre_norm, z, branch_cut, res);
                            },
                            [](xsf::numpy::autodiff2_cfloat z, int branch_cut, xsf::numpy::autodiff2_cfloat_2d res) {
                                xsf::assoc_legendre_p_all(xsf::assoc_legendre_norm, z, branch_cut, res);
                            }},
                           1, "assoc_legendre_p_all", nullptr, "(),()->(np1,mpmp1,3)", assoc_legendre_map_dims<1>),
        Py_False, 0,
        xsf::numpy::gufunc({xsf::numpy::compose{xsf::numpy::autodiff(), xsf::numpy::use_long_long_int()},
                            [](xsf::numpy::autodiff0_double z, int branch_cut, xsf::numpy::autodiff0_double_2d res) {
                                xsf::assoc_legendre_p_all(xsf::assoc_legendre_unnorm, z, branch_cut, res);
                            },
                            [](xsf::numpy::autodiff0_float z, int branch_cut, xsf::numpy::autodiff0_float_2d res) {
                                xsf::assoc_legendre_p_all(xsf::assoc_legendre_unnorm, z, branch_cut, res);
                            },
                            [](xsf::numpy::autodiff0_cdouble z, int branch_cut, xsf::numpy::autodiff0_cdouble_2d res) {
                                xsf::assoc_legendre_p_all(xsf::assoc_legendre_unnorm, z, branch_cut, res);
                            },
                            [](xsf::numpy::autodiff0_cfloat z, int branch_cut, xsf::numpy::autodiff0_cfloat_2d res) {
                                xsf::assoc_legendre_p_all(xsf::assoc_legendre_unnorm, z, branch_cut, res);
                            }},
                           1, "assoc_legendre_p_all", nullptr, "(),()->(np1,mpmp1,1)", assoc_legendre_map_dims<1>),
        Py_False, 1,
        xsf::numpy::gufunc({xsf::numpy::compose{xsf::numpy::autodiff(), xsf::numpy::use_long_long_int()},
                            [](xsf::numpy::autodiff1_double z, int branch_cut, xsf::numpy::autodiff1_double_2d res) {
                                xsf::assoc_legendre_p_all(xsf::assoc_legendre_unnorm, z, branch_cut, res);
                            },
                            [](xsf::numpy::autodiff1_float z, int branch_cut, xsf::numpy::autodiff1_float_2d res) {
                                xsf::assoc_legendre_p_all(xsf::assoc_legendre_unnorm, z, branch_cut, res);
                            },
                            [](xsf::numpy::autodiff1_cdouble z, int branch_cut, xsf::numpy::autodiff1_cdouble_2d res) {
                                xsf::assoc_legendre_p_all(xsf::assoc_legendre_unnorm, z, branch_cut, res);
                            },
                            [](xsf::numpy::autodiff1_cfloat z, int branch_cut, xsf::numpy::autodiff1_cfloat_2d res) {
                                xsf::assoc_legendre_p_all(xsf::assoc_legendre_unnorm, z, branch_cut, res);
                            }},
                           1, "assoc_legendre_p_all", nullptr, "(),()->(np1,mpmp1,2)", assoc_legendre_map_dims<1>),
        Py_False, 2,
        xsf::numpy::gufunc({xsf::numpy::compose{xsf::numpy::autodiff(), xsf::numpy::use_long_long_int()},
                            [](xsf::numpy::autodiff2_double z, int branch_cut, xsf::numpy::autodiff2_double_2d res) {
                                xsf::assoc_legendre_p_all(xsf::assoc_legendre_unnorm, z, branch_cut, res);
                            },
                            [](xsf::numpy::autodiff2_float z, int branch_cut, xsf::numpy::autodiff2_float_2d res) {
                                xsf::assoc_legendre_p_all(xsf::assoc_legendre_unnorm, z, branch_cut, res);
                            },
                            [](xsf::numpy::autodiff2_cdouble z, int branch_cut, xsf::numpy::autodiff2_cdouble_2d res) {
                                xsf::assoc_legendre_p_all(xsf::assoc_legendre_unnorm, z, branch_cut, res);
                            },
                            [](xsf::numpy::autodiff2_cfloat z, int branch_cut, xsf::numpy::autodiff2_cfloat_2d res) {
                                xsf::assoc_legendre_p_all(xsf::assoc_legendre_unnorm, z, branch_cut, res);
                            }},
                           1, "assoc_legendre_p_all", nullptr, "(),()->(np1,mpmp1,3)", assoc_legendre_map_dims<1>));
    PyModule_AddObjectRef(_gufuncs, "assoc_legendre_p_all", assoc_legendre_p_all);

    PyObject *sph_legendre_p_all = Py_BuildValue(
        "(N, N, N)",
        xsf::numpy::gufunc({xsf::numpy::compose{xsf::numpy::autodiff()},
                            static_cast<xsf::numpy::autodiff0_d_d2>(xsf::sph_legendre_p_all),
                            static_cast<xsf::numpy::autodiff0_f_f2>(xsf::sph_legendre_p_all)},
                           1, "sph_legendre_p_all", nullptr, "()->(np1,mpmp1,1)", assoc_legendre_map_dims<1>),
        xsf::numpy::gufunc({xsf::numpy::compose{xsf::numpy::autodiff()},
                            static_cast<xsf::numpy::autodiff1_d_d2>(xsf::sph_legendre_p_all),
                            static_cast<xsf::numpy::autodiff1_f_f2>(xsf::sph_legendre_p_all)},
                           1, "sph_legendre_p_all", nullptr, "()->(np1,mpmp1,2)", assoc_legendre_map_dims<1>),
        xsf::numpy::gufunc({xsf::numpy::compose{xsf::numpy::autodiff()},
                            static_cast<xsf::numpy::autodiff2_d_d2>(xsf::sph_legendre_p_all),
                            static_cast<xsf::numpy::autodiff2_f_f2>(xsf::sph_legendre_p_all)},
                           1, "sph_legendre_p_all", nullptr, "()->(np1,mpmp1,3)", assoc_legendre_map_dims<1>));
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

    PyObject *sph_harm_y_all =
        Py_BuildValue("(N, N, N)",
                      xsf::numpy::gufunc({xsf::numpy::compose{xsf::numpy::autodiff()},
                                          static_cast<xsf::numpy::autodiff00_dd_D2>(xsf::sph_harm_y_all),
                                          static_cast<xsf::numpy::autodiff00_ff_F2>(xsf::sph_harm_y_all)},
                                         1, "sph_harm_y_all", nullptr, "(),()->(np1,mpmp1,1,1)", sph_harm_map_dims),
                      xsf::numpy::gufunc({xsf::numpy::compose{xsf::numpy::autodiff()},
                                          static_cast<xsf::numpy::autodiff11_dd_D2>(xsf::sph_harm_y_all),
                                          static_cast<xsf::numpy::autodiff11_ff_F2>(xsf::sph_harm_y_all)},
                                         1, "sph_harm_y_all", nullptr, "(),()->(np1,mpmp1,2,2)", sph_harm_map_dims),
                      xsf::numpy::gufunc({xsf::numpy::compose{xsf::numpy::autodiff()},
                                          static_cast<xsf::numpy::autodiff22_dd_D2>(xsf::sph_harm_y_all),
                                          static_cast<xsf::numpy::autodiff22_ff_F2>(xsf::sph_harm_y_all)},
                                         1, "sph_harm_y_all", nullptr, "(),()->(np1,mpmp1,3,3)", sph_harm_map_dims));
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
