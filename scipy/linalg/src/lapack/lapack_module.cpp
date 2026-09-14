/**
 * @file
 * @brief Module definition for the LAPACK wrappers (`_flapack` / `_flapack_64`).
 *
 * This is the only file here that calls `_import_array()`; every other .cpp defines
 * `NO_IMPORT_ARRAY` so the extension shares one numpy C-API table.
 */
#define PY_ARRAY_UNIQUE_SYMBOL scipy_lapack_ARRAY_API
#include "lapack_helpers.hpp"

/* The ILP64 build differs from the LP64 one only in the module's name. */
#ifdef HAVE_BLAS_ILP64
#define FLAPACK_MODULE_STRING "_flapack_64"
#define FLAPACK_PYINIT        PyInit__flapack_64
#else
#define FLAPACK_MODULE_STRING "_flapack"
#define FLAPACK_PYINIT        PyInit__flapack
#endif

/* One table per `.pyf.src` group, each contributed by the matching wrapper .cpp.  Adding a
 * group means one more `extern` here and one more line in the exec slot below. */
namespace lapack {
namespace capi {
    extern PyMethodDef gen_methods[];
    extern PyMethodDef gen_tri_methods[];
    extern PyMethodDef gen_banded_methods[];
    extern PyMethodDef pos_def_methods[];
    extern PyMethodDef pos_def_tri_methods[];
    extern PyMethodDef sym_herm_methods[];
    extern PyMethodDef other_methods[];
    PyObject *build_doc(const char *name) noexcept;
}
}


/**
 * @brief A PyCFunction-shaped callable that also carries attributes, like f2py's `fortran` type.
 *
 * `scipy.linalg.lapack.get_lapack_funcs` *assigns* to the routines it hands back --
 * `func.module_name, func.typecode, func.dtype = ...` in `blas.py::_get_funcs` -- which a plain
 * `PyCFunction` cannot accept, having no instance dict.  This mirrors `_fblas`'s `BlasFunc`.
 *
 * @var meth  The wrapped LAPACK routine; called via its METH_KEYWORDS 3-arg form.
 * @var name  Flavored function name, `"dgesvx"` -- also exposed as `__name__`.
 * @var dict  Instance dict, lazily used by GenericGet/SetAttr.
 * @var doc   `__doc__`, built on first access and cached (see lapackfunc_get_doc).
 */
typedef struct {
    PyObject_HEAD
    PyCFunctionWithKeywords meth;
    PyObject *name;
    PyObject *dict;
    PyObject *doc;
} LapackFunc;


static PyObject *lapackfunc_call(PyObject *self, PyObject *args, PyObject *kwds) {
    return ((LapackFunc *)self)->meth(nullptr, args, kwds);
}


static PyObject *lapackfunc_repr(PyObject *self) {
    return PyUnicode_FromFormat("<flapack function %U>", ((LapackFunc *)self)->name);
}


static int lapackfunc_traverse(PyObject *self, visitproc visit, void *arg) {
    Py_VISIT(((LapackFunc *)self)->dict);
    Py_VISIT(((LapackFunc *)self)->doc);
    return 0;
}


static int lapackfunc_clear(PyObject *self) {
    Py_CLEAR(((LapackFunc *)self)->dict);
    Py_CLEAR(((LapackFunc *)self)->doc);
    return 0;
}


static void lapackfunc_dealloc(PyObject *self) {
    PyTypeObject *tp = Py_TYPE(self);
    PyObject_GC_UnTrack(self);
    Py_CLEAR(((LapackFunc *)self)->name);
    Py_CLEAR(((LapackFunc *)self)->dict);
    Py_CLEAR(((LapackFunc *)self)->doc);
    tp->tp_free(self);
    Py_DECREF(tp);   // heap type: instances own a reference to their type
}


static PyObject *lapackfunc_get_name(PyObject *self, void *Py_UNUSED(closure)) {
    return Py_NewRef(((LapackFunc *)self)->name);
}


/**
 * @brief Getset getter for a routine's `__doc__`, built once on demand and cached.
 *
 * @note Built lazily by build_doc() on first access and cached in `self->doc`. The store runs in
 *       a per-object critical section, so concurrent first accesses keep one build and drop the
 *       other. `doc` is set at most once and never cleared while `self` is live, so the cached
 *       fast path reads without locking.  Routines with no docstring registered return None.
 */
static PyObject *lapackfunc_get_doc(PyObject *self, void *Py_UNUSED(closure)) {
    LapackFunc *f = (LapackFunc *)self;
    if (f->doc != nullptr) { return Py_NewRef(f->doc); }

    const char *name = PyUnicode_AsUTF8(f->name);
    if (name == nullptr) { return nullptr; }

    PyObject *built = lapack::capi::build_doc(name);
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


static PyGetSetDef lapackfunc_getset[] = {
    {"__name__", lapackfunc_get_name, nullptr, nullptr, nullptr},
    {"__doc__", lapackfunc_get_doc, nullptr, nullptr, nullptr},
    {nullptr, nullptr, nullptr, nullptr, nullptr},
};


static PyMemberDef lapackfunc_members[] = {
    {"__dictoffset__", Py_T_PYSSIZET, offsetof(LapackFunc, dict), Py_READONLY, nullptr},
    {nullptr, 0, 0, 0, nullptr},
};


static PyType_Slot lapackfunc_slots[] = {
    {Py_tp_call,      reinterpret_cast<void *>(lapackfunc_call)},
    {Py_tp_repr,      reinterpret_cast<void *>(lapackfunc_repr)},
    {Py_tp_traverse,  reinterpret_cast<void *>(lapackfunc_traverse)},
    {Py_tp_clear,     reinterpret_cast<void *>(lapackfunc_clear)},
    {Py_tp_dealloc,   reinterpret_cast<void *>(lapackfunc_dealloc)},
    {Py_tp_getattro,  reinterpret_cast<void *>(PyObject_GenericGetAttr)},
    {Py_tp_setattro,  reinterpret_cast<void *>(PyObject_GenericSetAttr)},
    {Py_tp_getset,    reinterpret_cast<void *>(lapackfunc_getset)},
    {Py_tp_members,   reinterpret_cast<void *>(lapackfunc_members)},
    {0, nullptr},
};

static PyType_Spec lapackfunc_spec = {
    "scipy.linalg." FLAPACK_MODULE_STRING ".lapack_function", /* name      */
    sizeof(LapackFunc),                                           /* basicsize */
    0,                                                            /* itemsize  */
    Py_TPFLAGS_DEFAULT | Py_TPFLAGS_HAVE_GC,                      /* flags     */
    lapackfunc_slots,                                             /* slots     */
};


/** @brief Wrap every row of a PyMethodDef table in a LapackFunc and add it to the module. */
static int add_wrapped_table(PyObject *module, PyTypeObject *tp, const PyMethodDef *defs)
{
    for (const PyMethodDef *d = defs; d->ml_name != nullptr; d++) {
        LapackFunc *f = PyObject_GC_New(LapackFunc, tp);
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


static int lapack_module_exec(PyObject *module)
{
    if (_import_array() < 0) { return -1; }

    /* Heap type, created per interpreter; the instances keep it alive, so the module itself
     * stays stateless and drops its reference right away. */
    PyObject *tp_obj = PyType_FromModuleAndSpec(module, &lapackfunc_spec, nullptr);
    if (tp_obj == nullptr) { return -1; }
    PyTypeObject *tp = reinterpret_cast<PyTypeObject *>(tp_obj);

    int rc = add_wrapped_table(module, tp, lapack::capi::gen_methods);
    if (rc == 0) { rc = add_wrapped_table(module, tp, lapack::capi::gen_tri_methods); }
    if (rc == 0) { rc = add_wrapped_table(module, tp, lapack::capi::gen_banded_methods); }
    if (rc == 0) { rc = add_wrapped_table(module, tp, lapack::capi::pos_def_methods); }
    if (rc == 0) { rc = add_wrapped_table(module, tp, lapack::capi::pos_def_tri_methods); }
    if (rc == 0) { rc = add_wrapped_table(module, tp, lapack::capi::sym_herm_methods); }
    if (rc == 0) { rc = add_wrapped_table(module, tp, lapack::capi::other_methods); }
    Py_DECREF(tp_obj);
    return rc;
}

/* The wrappers keep no module or process state; the one callback slot is thread-local. */
static PyModuleDef_Slot lapack_module_slots[] = {
    {Py_mod_exec, reinterpret_cast<void *>(lapack_module_exec)},
    {Py_mod_multiple_interpreters, Py_MOD_PER_INTERPRETER_GIL_SUPPORTED},
#if PY_VERSION_HEX >= 0x030d00f0  /* Python 3.13+ */
    {Py_mod_gil, Py_MOD_GIL_NOT_USED},
#endif
    {0, nullptr},
};

static PyModuleDef lapack_moduledef = {
    PyModuleDef_HEAD_INIT,
    FLAPACK_MODULE_STRING,  /* m_name    */
    nullptr,                    /* m_doc     */
    0,                          /* m_size    */
    nullptr,                    /* m_methods, added in the exec slot */
    lapack_module_slots,        /* m_slots   */
    nullptr,                    /* m_traverse */
    nullptr,                    /* m_clear   */
    nullptr,                    /* m_free    */
};

PyMODINIT_FUNC FLAPACK_PYINIT(void)
{
    return PyModuleDef_Init(&lapack_moduledef);
}
