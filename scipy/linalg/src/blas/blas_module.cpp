/**
 * @file
 * @brief The `_fblas` / `_fblas_64` extension module: assembles the per-level wrapper tables.
 *
 * The wrappers live in `blas_l1.cpp` / `blas_l2.cpp` / `blas_l3.cpp` contributing a
 * `blas::capi::l*_methods` chunk merged in the exec slot with `PyModule_AddFunctions`.
 *
 */
#define PY_ARRAY_UNIQUE_SYMBOL scipy_blas_ARRAY_API
#include <Python.h>
#include <cstddef>
#include "numpy/arrayobject.h"

namespace blas{
    namespace capi {
        extern PyMethodDef l1_methods[], l2_methods[], l3_methods[];
        PyObject *build_doc(const char *name) noexcept;
    }
}


/**
 * The LP64 build generates the module `_fblas` and the ILP64 build generates `_fblas_64`,
 * mirroring the legacy f2py naming. The ILP64 lapack dependency passes `-DHAVE_BLAS_ILP64`,
 * which already selects `CBLAS_INT = int64_t` and the ILP64 `BLAS_FUNC` symbol suffix in
 * `scipy_blas_defines.h`; the module name is the only difference so it is selected here.
 */
#ifdef HAVE_BLAS_ILP64
#define FBLAS_MODULE_NAME _fblas_64
#else
#define FBLAS_MODULE_NAME _fblas
#endif

#define FBLAS_PASTE_(a, b) a ## b
#define FBLAS_PASTE(a, b) FBLAS_PASTE_(a, b)
#define FBLAS_STR_(s) #s
#define FBLAS_STR(s) FBLAS_STR_(s)

#define FBLAS_MODULE_STR FBLAS_STR(FBLAS_MODULE_NAME)      /* "_fblas" or "_fblas_64" */
#define FBLAS_PYINIT     FBLAS_PASTE(PyInit_, FBLAS_MODULE_NAME)


/**
 * @brief  BlasFunc data structure is a PyCFunction-shaped callable plus more attributes.
 *
 * This struct and the related object type mimics the behavior of the tool f2py which,
 * historically was used to wrap Fortran symbols via a `fortran` object. It also had custom
 * attributes, `__dict__`, `__name__`, and `__doc__` such that `scipy.linalg.blas.get_blas_funcs`
 * could query  attributes and set them.
 *
 * @var meth  The wrapped BLAS routine; called via its METH_KEYWORDS 3-arg form.
 * @var name  Flavored function name, `"daxpy"` -- also exposed as `__name__`.
 * @var dict  Instance dict, lazily used by GenericGet/SetAttr.
 * @var doc   `__doc__`, built on first access and cached (see blasfunc_get_doc).
 *
 */
typedef struct {
    PyObject_HEAD
    PyCFunctionWithKeywords meth;
    PyObject *name;
    PyObject *dict;
    PyObject *doc;
} BlasFunc;


// Now fill in the methods for the BlasFunc type.

static PyObject *blasfunc_call(PyObject *self, PyObject *args, PyObject *kwds) {
    return ((BlasFunc *)self)->meth(nullptr, args, kwds);
}


static PyObject *
blasfunc_repr(PyObject *self) {
    return PyUnicode_FromFormat("<fblas function %U>", ((BlasFunc *)self)->name);
}


static int
blasfunc_traverse(PyObject *self, visitproc visit, void *arg) {
    Py_VISIT(((BlasFunc *)self)->dict);
    Py_VISIT(((BlasFunc *)self)->doc);
    return 0;
}


static int
blasfunc_clear(PyObject *self) {
    Py_CLEAR(((BlasFunc *)self)->dict);
    Py_CLEAR(((BlasFunc *)self)->doc);
    return 0;
}


static void
blasfunc_dealloc(PyObject *self) {
    PyTypeObject *tp = Py_TYPE(self);
    PyObject_GC_UnTrack(self);
    Py_CLEAR(((BlasFunc *)self)->name);
    Py_CLEAR(((BlasFunc *)self)->dict);
    Py_CLEAR(((BlasFunc *)self)->doc);
    tp->tp_free(self);
    Py_DECREF(tp);   // heap type: instances own a reference to their type
}


static PyObject *
blasfunc_get_name(PyObject *self, void *Py_UNUSED(closure)) {
    return Py_NewRef(((BlasFunc *)self)->name);
}


/**
 * @brief Getset getter for a routine's `__doc__`, built once on demand and cached.
 *
 * @param self     The `BlasFunc` whose docstring is requested.
 * @param closure  Unused getset closure.
 * @return New reference to the docstring `str`; `None` if the routine has no docstring
 *         registered; `nullptr` with an exception set on failure.
 *
 * @note Built lazily by build_doc() on first access and cached in `self->doc`. The store runs
 *       in a per-object critical section, so concurrent first accesses keep one build and drop
 *       the other. `doc` is set at most once and never cleared while `self` is live, so the
 *       cached fast path reads without locking.
 */
static PyObject *
blasfunc_get_doc(PyObject *self, void *Py_UNUSED(closure)) {
    BlasFunc *f = (BlasFunc *)self;
    if (f->doc != nullptr) { return Py_NewRef(f->doc); }

    const char *name = PyUnicode_AsUTF8(f->name);
    if (name == nullptr) { return nullptr; }

    PyObject *built = blas::capi::build_doc(name);
    if (built == nullptr) {
        if (PyErr_Occurred()) { return nullptr; }
        Py_RETURN_NONE;   // no docstring registered for this routine
    }

#if PY_VERSION_HEX >= 0x030d00f0
    Py_BEGIN_CRITICAL_SECTION(self);
#endif
    if (f->doc == nullptr) {
        f->doc = built;
        built = nullptr;
    }
#if PY_VERSION_HEX >= 0x030d00f0
    Py_END_CRITICAL_SECTION();
#endif
    if (built != nullptr) { Py_DECREF(built); }   // lost the race; keep the winner
    return Py_NewRef(f->doc);
}


static PyGetSetDef blasfunc_getset[] = {
    {"__name__", blasfunc_get_name, nullptr, nullptr, nullptr},
    {"__doc__", blasfunc_get_doc, nullptr, nullptr, nullptr},
    {nullptr, nullptr, nullptr, nullptr, nullptr},
};


static PyMemberDef blasfunc_members[] = {
    {"__dictoffset__", Py_T_PYSSIZET, offsetof(BlasFunc, dict), Py_READONLY, nullptr},
    {nullptr, 0, 0, 0, nullptr},
};


static PyType_Slot blasfunc_slots[] = {
    {Py_tp_call,      reinterpret_cast<void *>(blasfunc_call)},
    {Py_tp_repr,      reinterpret_cast<void *>(blasfunc_repr)},
    {Py_tp_traverse,  reinterpret_cast<void *>(blasfunc_traverse)},
    {Py_tp_clear,     reinterpret_cast<void *>(blasfunc_clear)},
    {Py_tp_dealloc,   reinterpret_cast<void *>(blasfunc_dealloc)},
    {Py_tp_getattro,  reinterpret_cast<void *>(PyObject_GenericGetAttr)},
    {Py_tp_setattro,  reinterpret_cast<void *>(PyObject_GenericSetAttr)},
    {Py_tp_getset,    reinterpret_cast<void *>(blasfunc_getset)},
    {Py_tp_members,   reinterpret_cast<void *>(blasfunc_members)},
    {0, nullptr},
};

static PyType_Spec blasfunc_spec = {
    "scipy.linalg." FBLAS_MODULE_STR ".blas_function", /* name      */
    sizeof(BlasFunc),                                  /* basicsize */
    0,                                                 /* itemsize  */
    Py_TPFLAGS_DEFAULT | Py_TPFLAGS_HAVE_GC,           /* flags     */
    blasfunc_slots,                                    /* slots     */
};

/** @brief Wrap every row of a PyMethodDef table in a BlasFunc and add it to the module. */
static int
add_wrapped_table(PyObject *module, PyTypeObject *tp, const PyMethodDef *defs) {

    for (const PyMethodDef *d = defs; d->ml_name != nullptr; d++) {

        BlasFunc *f = PyObject_GC_New(BlasFunc, tp);
        if (f == nullptr) { return -1; }

        Py_INCREF(tp);   // the reference the instance owns

        f->meth = reinterpret_cast<PyCFunctionWithKeywords>(reinterpret_cast<void (*)()>(d->ml_meth));
        f->dict = nullptr;
        f->doc = nullptr;
        f->name = PyUnicode_FromString(d->ml_name);
        if (f->name == nullptr) { Py_DECREF(f); return -1; }

        PyObject_GC_Track(reinterpret_cast<PyObject *>(f));
        if (PyModule_AddObject(module, d->ml_name, reinterpret_cast<PyObject *>(f)) < 0) {
            Py_DECREF(f);
            return -1;
        }
    }
    return 0;
}


static int
_blas_module_exec(PyObject *module)
{
    if (_import_array() < 0) { return -1; }

    /* Heap type, created per interpreter; the instances keep it alive, so the module
     * itself stays stateless and drops its reference right away. */
    PyObject *tp_obj = PyType_FromModuleAndSpec(module, &blasfunc_spec, nullptr);
    if (tp_obj == nullptr) { return -1; }
    PyTypeObject *tp = reinterpret_cast<PyTypeObject *>(tp_obj);

    int rc = 0;
    if (add_wrapped_table(module, tp, blas::capi::l1_methods) < 0 ||
        add_wrapped_table(module, tp, blas::capi::l2_methods) < 0 ||
        add_wrapped_table(module, tp, blas::capi::l3_methods) < 0) {
        rc = -1;
    }
    Py_DECREF(tp_obj);
    return rc;
}


static PyModuleDef_Slot _blas_slots[] = {
    {Py_mod_exec, reinterpret_cast<void *>(_blas_module_exec)},
    {Py_mod_multiple_interpreters, Py_MOD_PER_INTERPRETER_GIL_SUPPORTED},
#if PY_VERSION_HEX >= 0x030d00f0  /* Python 3.13+ */
    {Py_mod_gil, Py_MOD_GIL_NOT_USED},
#endif
    {0, nullptr},
};


// Designated initializers require C99 or C++20.
// CAUTION: Even then, C++ requires the order of the fields preserved.
static struct PyModuleDef _blas_moduledef = {
    PyModuleDef_HEAD_INIT,    /* m_base     */
    FBLAS_MODULE_STR,         /* m_name: "_fblas" or "_fblas_64" (ABI-derived above) */
    nullptr,                  /* m_doc      */
    0,                        /* m_size     */
    nullptr,                  /* m_methods, merged per level in the exec slot */
    _blas_slots,              /* m_slots    */
    nullptr,                  /* m_traverse */
    nullptr,                  /* m_clear    */
    nullptr                   /* m_free     */
};


PyMODINIT_FUNC
FBLAS_PYINIT(void)   /* PyInit__fblas (LP64) or PyInit__fblas_64 (ILP64) */
{
    return PyModuleDef_Init(&_blas_moduledef);
}
