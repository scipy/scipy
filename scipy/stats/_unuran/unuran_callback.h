#include <Python.h>
#include <setjmp.h>
#include "ccallback.h"
#include "unuran.h"

/*
 * A common UNU.RAN thunk.
 * 
 * Parameters
 * ----------
 * CAST_FUNC :
 *     CPython function to convert the C data type of the input i.e. data type of
 *     the quantile (`ARG`) into a PyObject.
 * CAST_BACK_FUNC :
 *     CPython function to convert the PyObject result to a C data type. i.e the
 *     return type of the thunk.
 * FUNCNAME :
 *     The name of the function to call in the thunk. This name must be present in
 *     the `callbacks` dictionary.
 * ARG :
 *     The argument to treat as an input to the Python function.
 * 
 * Returns
 * -------
 * result :
 *     The value returned by the Python function after converting it into a C type
 *     using `CAST_BACK_FUNC`.
 */
#define UNURAN_THUNK(CAST_FUNC, CAST_BACK_FUNC, FUNCNAME, ARG)                         \
    ccallback_t *err_callback = ccallback_obtain();                                    \
    unuran_callback_t *unur_callback = (unuran_callback_t *)(err_callback->info_p);    \
    PyObject *extra_arguments = unur_callback->params;                                 \
    PyObject *py_function = PyDict_GetItemString(unur_callback->callbacks, FUNCNAME);  \
                                                                                       \
    PyObject *arg1 = NULL, *argobj = NULL, *arglist = NULL, *res = NULL;               \
    double result = 0.;                                                                \
    int error = 0;                                                                     \
                                                                                       \
    argobj = CAST_FUNC(ARG);                                                           \
    if (argobj == NULL) {                                                              \
        error = 1;                                                                     \
        goto done;                                                                     \
    }                                                                                  \
                                                                                       \
    arg1 = PyTuple_New(1);                                                             \
    if (arg1 == NULL) {                                                                \
        error = 1;                                                                     \
        goto done;                                                                     \
    }                                                                                  \
                                                                                       \
    PyTuple_SET_ITEM(arg1, 0, argobj);                                                 \
    argobj = NULL;                                                                     \
                                                                                       \
    arglist = PySequence_Concat(arg1, extra_arguments);                                \
    if (arglist == NULL) {                                                             \
        error = 1;                                                                     \
        goto done;                                                                     \
    }                                                                                  \
                                                                                       \
    res = PyObject_CallObject(py_function, arglist);                                   \
    if (res == NULL) {                                                                 \
        error = 1;                                                                     \
        goto done;                                                                     \
    }                                                                                  \
                                                                                       \
    result = CAST_BACK_FUNC(res);                                                      \
    if (PyErr_Occurred()) {                                                            \
        error = 1;                                                                     \
        goto done;                                                                     \
    }                                                                                  \
                                                                                       \
done:                                                                                  \
    Py_XDECREF(arg1);                                                                  \
    Py_XDECREF(argobj);                                                                \
    Py_XDECREF(arglist);                                                               \
    Py_XDECREF(res);                                                                   \
                                                                                       \
    if (error) {                                                                       \
        longjmp(err_callback->error_buf, 1);                                           \
    }                                                                                  \
                                                                                       \
    return result


/* A structure to represent a UNU.RAN Python callback */
typedef struct unuran_callback {
    /* A dictionary with required callbacks */
    PyObject *callbacks;
    /* Parameters (other than the quantile parameter) that the functions take */
    PyObject *params;
} unuran_callback_t;

static ccallback_signature_t unuran_call_signatures[] = {
    {NULL}
};

/**
 * Initialize a UNU.RAN callback.
 * 
 * This function prepares a `ccallback_t` object containing the Python callbacks and
 * their parameters. This helps obtain and call Python functions inside a C thunk in
 * a thread-safe manner. Nonetheless, a callback needs to be aquired even when there
 * are no Python callbacks because the error handler needs a `ccallback_t` object to
 * jump to the Cython code where errors can be raised safely. In such a case, one can
 * pass an empty dictionary in `fcn_dict` parameter.
 * 
 * Parameters
 * ----------
 * callback : ccallback_t *
 *     Reference to an empty C callback
 * unur_callback : unuran_callback_t *
 *     Reference to an empty UNU.RAN callback
 * fcn_dict : PyObject *
 *     A Python dictionary containing the Python callbacks which are to be called in
 *     the thunk. Can be an empty dict in case of no Python callbacks.
 * extra_args : PyObject *
 *     A tuple of extra arguments that need to be passed to the functions present in
 *     the `fcn_dict` dictionary. If there are none, an empty tuple is expected.
 * 
 * Returns
 * -------
 * success : int
 *     0 if success -1 otherwise
 */
int init_unuran_callback(ccallback_t *callback, unuran_callback_t *unur_callback,
                         PyObject *fcn_dict, PyObject *extra_args)
{
    PyObject *dummy_fcn = NULL;
    int ret;
    int flags = CCALLBACK_OBTAIN;

    callback->c_function = NULL;
    callback->user_data = NULL;
    callback->signature = NULL;

    unur_callback->callbacks = fcn_dict;
    Py_INCREF(fcn_dict);
    unur_callback->params = extra_args;
    Py_INCREF(extra_args);

    /* The required callbacks are present in the `fcn_dict` dictionary. There can
       be no callbacks for UNU.RAN yet we would need a C callback for error handler
       to pass on the captured UNU.RAN errors to Python. Hence, we pack a dummy
       function in the C callback object and use `unuran_callback_t` to obtain
       required callbacks, if any. */
    dummy_fcn = PyRun_String("lambda: None", Py_eval_input, PyEval_GetGlobals(), NULL);
    if (dummy_fcn == NULL) {
        goto fail;
    }

    ret = ccallback_prepare(callback, unuran_call_signatures, dummy_fcn, flags);
    if (ret == -1) {
        goto fail;
    }

    callback->info_p = (void *)unur_callback;

    return 0;

fail:
    Py_DECREF(fcn_dict);
    Py_DECREF(extra_args);
    Py_XDECREF(dummy_fcn);
    return -1;
}

/**
 * Release a UNU.RAN callback.
 * 
 * Parameters
 * ----------
 * callback : ccallback_t *
 *     C callback to be released
 * unur_callback : unuran_callback_t *
 *     UNU.RAN callback to be released
 * 
 * Returns
 * -------
 * success : int
 *     0 if success -1 otherwise
 */
int release_unuran_callback(ccallback_t *callback, unuran_callback_t *unur_callback) {
    Py_XDECREF(unur_callback->callbacks);
    Py_XDECREF(unur_callback->params);
    Py_XDECREF(callback->py_function);
    unur_callback->callbacks = NULL;
    unur_callback->params = NULL;
    int ret = ccallback_release(callback);
    return ret;
}


void error_handler(const char *objid, const char *file, int line, const char *errortype,
                   int unur_errno, const char *reason)
{
    ccallback_t *err_callback;
    const char *errno_msg;
    char reason_[256], objid_[256], errno_msg_[256];

    if (unur_errno != UNUR_SUCCESS) {
        err_callback = ccallback_obtain();

        if (reason == NULL || strcmp(reason, "") == 0) {
            strcpy(reason_, "unknown error");
        }
        else {
            strcpy(reason_, reason);
        }
        if (objid == NULL || strcmp(objid, "") == 0) {
            strcpy(objid_, "unknown");
        }
        else {
            strcpy(objid_, objid);
        }

        errno_msg = unur_get_strerror(unur_errno);

        if (errno_msg == NULL || strcmp(errno_msg, "") == 0) {
            strcpy(errno_msg_, "unknown type");
        }
        else {
            strcpy(errno_msg_, errno_msg);
        }
        if (strcmp(errortype, "error") == 0) {
            PyErr_Format(PyExc_RuntimeError,
                         "[objid: %s] %d : %s => %s",
                         objid_, unur_errno, reason_,
                         errno_msg_);
            longjmp(err_callback->error_buf, 1);
        }
        else { /* warning */
            PyErr_WarnFormat(PyExc_UserWarning, 1,
                             "[objid: %s] %d : %s => %s",
                             objid_, unur_errno, reason_,
                             errno_msg_);
            if (PyErr_Occurred()) {
                longjmp(err_callback->error_buf, 1);
            }
        }
    }
}


/* ********************************************************************************** */
/* ********************************* UNU.RAN Thunks ********************************* */
/* ********************************************************************************** */

double pmf_thunk(int k, const struct unur_distr *distr)
{
    UNURAN_THUNK(PyLong_FromLong, PyFloat_AsDouble, "pmf", k);
}

double pdf_thunk(double x, const struct unur_distr *distr)
{
    UNURAN_THUNK(PyFloat_FromDouble, PyFloat_AsDouble, "pdf", x);
}

double dpdf_thunk(double x, const struct unur_distr *distr)
{
    UNURAN_THUNK(PyFloat_FromDouble, PyFloat_AsDouble, "dpdf", x);
}

double cont_cdf_thunk(double x, const struct unur_distr *distr)
{
    UNURAN_THUNK(PyFloat_FromDouble, PyFloat_AsDouble, "cdf", x);
}

double discr_cdf_thunk(int k, const struct unur_distr *distr)
{
    UNURAN_THUNK(PyLong_FromLong, PyFloat_AsDouble, "cdf", k);
}
