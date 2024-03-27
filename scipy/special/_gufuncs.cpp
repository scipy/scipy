#include <cmath>

#include "special/legendre.h"
#include "ufunc.h"

using namespace std;

using lpn_fff_t = void (*)(float, mdspan<float, dextents<ptrdiff_t, 1>, layout_stride>,
                           mdspan<float, dextents<ptrdiff_t, 1>, layout_stride>);
using lpn_ddd_t = void (*)(double, mdspan<double, dextents<ptrdiff_t, 1>, layout_stride>,
                           mdspan<double, dextents<ptrdiff_t, 1>, layout_stride>);
using lpn_FFF_t = void (*)(complex<float>, mdspan<complex<float>, dextents<ptrdiff_t, 1>, layout_stride>,
                           mdspan<complex<float>, dextents<ptrdiff_t, 1>, layout_stride>);
using lpn_DDD_t = void (*)(complex<double>, mdspan<complex<double>, dextents<ptrdiff_t, 1>, layout_stride>,
                           mdspan<complex<double>, dextents<ptrdiff_t, 1>, layout_stride>);
extern const char *lpn_doc;

using lpmn_flff_t = void (*)(float, long, mdspan<float, dextents<ptrdiff_t, 2>, layout_stride>,
                             mdspan<float, dextents<ptrdiff_t, 2>, layout_stride>);
using lpmn_dldd_t = void (*)(double, long, mdspan<double, dextents<ptrdiff_t, 2>, layout_stride>,
                             mdspan<double, dextents<ptrdiff_t, 2>, layout_stride>);
extern const char *lpmn_doc;

using clpmn_FllFF_t = void (*)(complex<float>, long, long,
                               mdspan<complex<float>, dextents<ptrdiff_t, 2>, layout_stride>,
                               mdspan<complex<float>, dextents<ptrdiff_t, 2>, layout_stride>);
using clpmn_DllDD_t = void (*)(complex<double>, long, long,
                               mdspan<complex<double>, dextents<ptrdiff_t, 2>, layout_stride>,
                               mdspan<complex<double>, dextents<ptrdiff_t, 2>, layout_stride>);
extern const char *clpmn_doc;

using lqn_fff_t = void (*)(float, mdspan<float, dextents<ptrdiff_t, 1>, layout_stride>,
                           mdspan<float, dextents<ptrdiff_t, 1>, layout_stride>);
using lqn_ddd_t = void (*)(double, mdspan<double, dextents<ptrdiff_t, 1>, layout_stride>,
                           mdspan<double, dextents<ptrdiff_t, 1>, layout_stride>);
using lqn_FFF_t = void (*)(complex<float>, mdspan<complex<float>, dextents<ptrdiff_t, 1>, layout_stride>,
                           mdspan<complex<float>, dextents<ptrdiff_t, 1>, layout_stride>);
using lqn_DDD_t = void (*)(complex<double>, mdspan<complex<double>, dextents<ptrdiff_t, 1>, layout_stride>,
                           mdspan<complex<double>, dextents<ptrdiff_t, 1>, layout_stride>);
extern const char *lqn_doc;

using lqmn_fff_t = void (*)(float, mdspan<float, dextents<ptrdiff_t, 2>, layout_stride>,
                            mdspan<float, dextents<ptrdiff_t, 2>, layout_stride>);
using lqmn_ddd_t = void (*)(double, mdspan<double, dextents<ptrdiff_t, 2>, layout_stride>,
                            mdspan<double, dextents<ptrdiff_t, 2>, layout_stride>);
using lqmn_FFF_t = void (*)(complex<float>, mdspan<complex<float>, dextents<ptrdiff_t, 2>, layout_stride>,
                            mdspan<complex<float>, dextents<ptrdiff_t, 2>, layout_stride>);
using lqmn_DDD_t = void (*)(complex<double>, mdspan<complex<double>, dextents<ptrdiff_t, 2>, layout_stride>,
                            mdspan<complex<double>, dextents<ptrdiff_t, 2>, layout_stride>);
extern const char *lqmn_doc;

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

    PyObject *_lpn = SpecFun_NewGUFunc({static_cast<lpn_fff_t>(special::lpn), static_cast<lpn_ddd_t>(special::lpn),
                                        static_cast<lpn_FFF_t>(special::lpn), static_cast<lpn_DDD_t>(special::lpn)},
                                       2, "_lpn", lpn_doc, "()->(np1),(np1)");
    PyModule_AddObjectRef(_gufuncs, "_lpn", _lpn);

    PyObject *_lpmn =
        SpecFun_NewGUFunc({static_cast<lpmn_flff_t>(special::lpmn), static_cast<lpmn_dldd_t>(special::lpmn)}, 2,
                          "_lpmn", lpmn_doc, "(),()->(mp1,np1),(mp1,np1)");
    PyModule_AddObjectRef(_gufuncs, "_lpmn", _lpmn);

    PyObject *_clpmn =
        SpecFun_NewGUFunc({static_cast<clpmn_FllFF_t>(special::clpmn), static_cast<clpmn_DllDD_t>(special::clpmn)}, 2,
                          "_clpmn", clpmn_doc, "(),(),()->(mp1,np1),(mp1,np1)");
    PyModule_AddObjectRef(_gufuncs, "_clpmn", _clpmn);

    PyObject *_lqn = SpecFun_NewGUFunc({static_cast<lqn_fff_t>(special::lqn), static_cast<lqn_ddd_t>(special::lqn),
                                        static_cast<lqn_FFF_t>(special::lqn), static_cast<lqn_DDD_t>(special::lqn)},
                                       2, "_lqn", lqn_doc, "()->(np1),(np1)");
    PyModule_AddObjectRef(_gufuncs, "_lqn", _lqn);

    PyObject *_lqmn =
        SpecFun_NewGUFunc({static_cast<lqmn_fff_t>(special::lqmn), static_cast<lqmn_ddd_t>(special::lqmn),
                           static_cast<lqmn_FFF_t>(special::lqmn), static_cast<lqmn_DDD_t>(special::lqmn)},
                          2, "_lqmn", lqmn_doc, "()->(mp1,np1),(mp1,np1)");
    PyModule_AddObjectRef(_gufuncs, "_lqmn", _lqmn);

    return _gufuncs;
}
