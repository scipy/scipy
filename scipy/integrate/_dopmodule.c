#define PY_SSIZE_T_CLEAN
#include "Python.h"
#include "numpy/arrayobject.h"
#include "ccallback.h"
#include "src/dop.h"
#include <math.h>

static PyObject* dop_error;

static char doc_dopri853[] = "x,y,iwork,idid = dop853(fcn,x,y,xend,rtol,atol,solout,iout,work,iwork,nsteps,verbosity,[fcn_extra_args])";
static char doc_dopri5[] = "x,y,iwork,idid = dopri5(fcn,x,y,xend,rtol,atol,solout,iout,work,iwork,nsteps,verbosity,[fcn_extra_args])";


typedef struct {
    PyObject *func;
    PyObject *solout;
    PyObject *func_args;
} callback_info_t;


// Thread-local storage for callbacks
static SCIPY_TLS callback_info_t* current_func_callback = NULL;


static void
func_thunk(int n, double t, double *y, double *f, double *rpar, int *ipar)
{
    // SciPy does not use rpar and ipar. Silence unused parameter warnings.
    (void)rpar; (void)ipar;

    if (!current_func_callback || !current_func_callback->func) { return; }
    npy_intp dims[1] = {n};

    // Create a PyCapsule to manage the data lifetime
    PyObject *capsule = PyCapsule_New((void*)y, NULL, NULL);
    if (!capsule) { return; }

    PyObject* py_y = PyArray_SimpleNewFromData(1, dims, NPY_FLOAT64, (void*)y);
    if (!py_y) { Py_DECREF(capsule); return; }

    // Tell NumPy that Capsule owns the memory so it does not try to free it.
    if (PyArray_SetBaseObject((PyArrayObject*)py_y, capsule) < 0) {
        // capsule reference is stolen even on failure
        Py_DECREF(py_y);
        return;
    }

    // We already prepared the tuple in the callback setup
    PyTuple_SetItem(current_func_callback->func_args, 0, PyFloat_FromDouble(t));
    PyTuple_SetItem(current_func_callback->func_args, 1, py_y);

    // Call Python function
    PyObject *result = PyObject_CallObject(current_func_callback->func, current_func_callback->func_args);


    if (result) {
        // Extract result directly to f array
        PyArrayObject *result_array = (PyArrayObject*)PyArray_FROM_OTF(result, NPY_DOUBLE, NPY_ARRAY_IN_ARRAY);
        if (result_array) {
            double *result_data = (double*)PyArray_DATA(result_array);
            for (int i = 0; i < n; i++) { f[i] = result_data[i]; }
            Py_DECREF(result_array);
        }
        Py_DECREF(result);
    }
}


static void
solout_thunk(int nr, double xold, double x, double* y, int n, double* con, int* icomp, int nd, double* rpar, int* ipar, int* irtrn)
{
    // SciPy does not use rpar and ipar. Silence unused parameter warnings.
    (void)rpar; (void)ipar;

    if (!current_func_callback || !current_func_callback->solout) { return; }

    npy_intp dims[1] = {n};

    // Create a PyCapsule to manage the data lifetime
    PyObject *capsule = PyCapsule_New((void*)y, NULL, NULL);
    if (!capsule) { return; }

    PyObject* py_y = PyArray_SimpleNewFromData(1, dims, NPY_FLOAT64, (void*)y);
    if (!py_y) { Py_DECREF(capsule); return; }

    // Tell NumPy that Capsule owns the memory so it does not try to free it.
    if (PyArray_SetBaseObject((PyArrayObject*)py_y, capsule) < 0) {
        // capsule reference is stolen even on failure
        Py_DECREF((PyObject*)py_y);
        return;
    }

    // Prepare the tuple for solout
    PyTuple_SetItem(current_func_callback->func_args, 0, PyFloat_FromDouble(x));
    PyTuple_SetItem(current_func_callback->func_args, 1, py_y);

    // Call solout function
    PyObject *result = PyObject_CallObject(current_func_callback->solout, current_func_callback->func_args);
    if (result) {
        // Extract result if needed, or just check for errors
        if (PyLong_Check(result)) {
            // Until we support Python >3.13 we can't use PyLong_AsInt
            *irtrn = (int)PyLong_AsLong(result);
        } else if (Py_IsNone(result))
        {
            *irtrn = 0;
        } else {
            PyErr_SetString(PyExc_TypeError, "scipy.integrate: solout function must return -1"
                                             " for stopping, and otherwise 0, 1 or None.");
        }
        Py_DECREF(result);
    }
}


static PyObject*
dopri853_wrap(PyObject* Py_UNUSED(self), PyObject* args)
{

    double x, xend, rtol, atol;
    int ierr = 0, itol = 0, iout = 0, nsteps = 0, verbosity = 0;
    PyObject *py_func, *solout, *func_args;
    PyArrayObject *ap_work = NULL, *ap_iwork = NULL, *ap_y0 = NULL;


    if (!PyArg_ParseTuple(args, "OdO!dddOiO!O!ii|O",
        &py_func,                                // O
        &x,                                      // d
        &PyArray_Type, (PyObject **)&ap_y0,      // O!
        &xend,                                   // d
        &rtol,                                   // d
        &atol,                                   // d
        &solout,                                 // O
        &iout,                                   // i
        &PyArray_Type, (PyObject **)&ap_work,    // O!
        &PyArray_Type, (PyObject **)&ap_iwork,   // O!
        &nsteps,                                 // i
        &verbosity,                              // i
        &func_args                               // O (optional)
    ))
    {
        return NULL;
    }

    PyArrayObject *ap_y = (PyArrayObject*)PyArray_NewCopy(ap_y0, NPY_CORDER);
    if (!ap_y) { return NULL; }

    double* y = (double*)PyArray_DATA(ap_y);
    double* work = (double*)PyArray_DATA(ap_work);

    // Caller zeros out the iwork
    int* iwork = (int*)PyArray_DATA(ap_iwork);
    iwork[0] = nsteps;
    iwork[2] = verbosity;

    int n = PyArray_DIMS(ap_y)[0];

    // Prepare the callbacks with optional arguments if any.
    callback_info_t func_callback;

    func_callback.func = py_func;
    Py_INCREF(py_func);

    func_callback.solout = solout;
    if (solout) {  Py_INCREF(solout); }

    func_callback.func_args = NULL;  // Will be set later

    // Mark the callbacks as thread local callback pointers
    current_func_callback = &func_callback;

    /*************************************************************

    In the callbacks both py_func and, if provided, solout will be called with the following
    signature:

        f = py_func(t, y, *func_args)
        res = solout(t, y, *func_args)

    Instead of unpacking and building tuples inefficiently, we define the tuple once
    here and change the value of the first two elements in the callbacks.

    We initially fill the tuple with None values and hence increment the reference count then
    PyTuple_SetItem will decrement the refcount of the replaced and increase the refcount of
    the new item automatically.

    *************************************************************/
    if (func_args && PyTuple_Check(func_args)) {
        // Create (None, None, *func_args)
        Py_ssize_t extra_size = PyTuple_Size(func_args);

        PyObject *args_tuple = PyTuple_New(2 + extra_size);
        PyTuple_SetItem(args_tuple, 0, Py_None); Py_INCREF(Py_None);
        PyTuple_SetItem(args_tuple, 1, Py_None); Py_INCREF(Py_None);
        // Copy extra args once
        for (Py_ssize_t i = 0; i < extra_size; i++) {
            PyObject *item = PyTuple_GetItem(func_args, i);
            Py_INCREF(item);
            PyTuple_SetItem(args_tuple, 2 + i, item);
        }
        func_callback.func_args = args_tuple;
    } else {
        // Create (None, None)
        PyObject *args_tuple = PyTuple_New(2);
        PyTuple_SetItem(args_tuple, 0, Py_None); Py_INCREF(Py_None);
        PyTuple_SetItem(args_tuple, 1, Py_None); Py_INCREF(Py_None);
        func_callback.func_args = args_tuple;
    }

    dopri853(n, func_thunk, &x, y, &xend, &rtol, &atol, itol, solout_thunk, iout, work, iwork, NULL, NULL, &ierr);

    if (func_callback.func_args) { Py_DECREF(func_callback.func_args); }
    current_func_callback = NULL;

    return Py_BuildValue("dNi", x, PyArray_Return(ap_y), ierr);
}


static PyObject*
dopri5_wrap(PyObject* Py_UNUSED(self), PyObject* args)
{

    double x, xend, rtol, atol;
    int ierr = 0, itol = 0, iout = 0, nsteps = 0, verbosity = 0;
    PyObject *py_func, *solout, *func_args;
    PyArrayObject *ap_work = NULL, *ap_iwork = NULL, *ap_y0 = NULL;


    if (!PyArg_ParseTuple(args, "OdO!dddOiO!O!ii|O",
        &py_func,                                // O
        &x,                                      // d
        &PyArray_Type, (PyObject **)&ap_y0,      // O!
        &xend,                                   // d
        &rtol,                                   // d
        &atol,                                   // d
        &solout,                                 // O
        &iout,                                   // i
        &PyArray_Type, (PyObject **)&ap_work,    // O!
        &PyArray_Type, (PyObject **)&ap_iwork,   // O!
        &nsteps,                                 // i
        &verbosity,                              // i
        &func_args                               // O (optional)
    ))
    {
        return NULL;
    }

    PyArrayObject *ap_y = (PyArrayObject*)PyArray_NewCopy(ap_y0, NPY_CORDER);
    if (!ap_y) { return NULL; }

    double* y = (double*)PyArray_DATA(ap_y);
    double* work = (double*)PyArray_DATA(ap_work);

    // Caller zeros out the iwork
    int* iwork = (int*)PyArray_DATA(ap_iwork);
    iwork[0] = nsteps;
    iwork[2] = verbosity;

    int n = PyArray_DIMS(ap_y)[0];

    // Prepare the callbacks with optional arguments if any.
    callback_info_t func_callback;

    func_callback.func = py_func;
    Py_INCREF(py_func);

    func_callback.solout = solout;
    if (solout) {  Py_INCREF(solout); }

    func_callback.func_args = NULL;  // Will be set later

    // Mark the callbacks as thread local callback pointers
    current_func_callback = &func_callback;

    /*************************************************************

    In the callbacks both py_func and, if provided, solout will be called with the following
    signature:

        f = py_func(t, y, *func_args)
        res = solout(t, y, *func_args)

    Instead of unpacking and building tuples inefficiently, we define the tuple once
    here and change the value of the first two elements in the callbacks.

    We initially fill the tuple with None values and hence increment the reference count then
    PyTuple_SetItem will decrement the refcount of the replaced and increase the refcount of
    the new item automatically.

    *************************************************************/
    if (func_args && PyTuple_Check(func_args)) {
        // Create (None, None, *func_args)
        Py_ssize_t extra_size = PyTuple_Size(func_args);

        PyObject *args_tuple = PyTuple_New(2 + extra_size);
        PyTuple_SetItem(args_tuple, 0, Py_None); Py_INCREF(Py_None);
        PyTuple_SetItem(args_tuple, 1, Py_None); Py_INCREF(Py_None);
        // Copy extra args once
        for (Py_ssize_t i = 0; i < extra_size; i++) {
            PyObject *item = PyTuple_GetItem(func_args, i);
            Py_INCREF(item);
            PyTuple_SetItem(args_tuple, 2 + i, item);
        }
        func_callback.func_args = args_tuple;
    } else {
        // Create (None, None)
        PyObject *args_tuple = PyTuple_New(2);
        PyTuple_SetItem(args_tuple, 0, Py_None); Py_INCREF(Py_None);
        PyTuple_SetItem(args_tuple, 1, Py_None); Py_INCREF(Py_None);
        func_callback.func_args = args_tuple;
    }

    dopri5(n, func_thunk, &x, y, &xend, &rtol, &atol, itol, solout_thunk, iout, work, iwork, NULL, NULL, &ierr);

    if (func_callback.func_args) { Py_DECREF(func_callback.func_args); }
    current_func_callback = NULL;

    return Py_BuildValue("dNi", x, PyArray_Return(ap_y), ierr);
}


static struct PyMethodDef doplib_module_methods[] = {
  {"dopri853", dopri853_wrap, METH_VARARGS, doc_dopri853},
  {"dopri5", dopri5_wrap, METH_VARARGS, doc_dopri5},
  {NULL    , NULL          , 0           , NULL     }
};


static int
doplib_module_exec(PyObject *module)
{
    if (_import_array() < 0) { return -1; }

    dop_error = PyErr_NewException("_dop.error", NULL, NULL);
    if (dop_error == NULL) { return -1; }

    if (PyModule_AddObject(module, "error", dop_error) < 0) {
        Py_DECREF(dop_error);
        return -1;
    }

    return 0;
}


static struct PyModuleDef_Slot doplib_module_slots[] = {
    {Py_mod_exec, doplib_module_exec},
#if PY_VERSION_HEX >= 0x030c00f0  // Python 3.12+
    // signal that this module can be imported in isolated subinterpreters
    {Py_mod_multiple_interpreters, Py_MOD_PER_INTERPRETER_GIL_SUPPORTED},
#endif
#if PY_VERSION_HEX >= 0x030d00f0  // Python 3.13+
    // signal that this module supports running without an active GIL
    {Py_mod_gil, Py_MOD_GIL_NOT_USED},
#endif
    {0, NULL},
};


static struct
PyModuleDef moduledef = {
    .m_base = PyModuleDef_HEAD_INIT,
    .m_name = "_dop",
    .m_size = 0,
    .m_methods = doplib_module_methods,
    .m_slots = doplib_module_slots,
};


PyMODINIT_FUNC
PyInit__dop(void)
{
    return PyModuleDef_Init(&moduledef);
}


