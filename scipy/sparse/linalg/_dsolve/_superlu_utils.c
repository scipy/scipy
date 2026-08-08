/* Should be imported before Python.h */

#include <Python.h>

#define NO_IMPORT_ARRAY
#define PY_ARRAY_UNIQUE_SYMBOL _scipy_sparse_superlu_ARRAY_API

#include "_superluobject.h"
#include <setjmp.h>


/* Abort to be used inside the superlu module so that memory allocation
   errors don't exit Python and memory allocated internal to SuperLU is freed.
   Calling program should deallocate (using SUPERLU_FREE) all memory that could have
   been allocated.  (It's ok to FREE unallocated memory)---will be ignored.
*/

SuperLUGlobalObject *get_tls_global(void)
{
    PyObject *thread_dict;
    SuperLUGlobalObject *obj;
    const char *key = "scipy.sparse.linalg._dsolve._superlu.__global_object";

    thread_dict = PyThreadState_GetDict();
    if (thread_dict == NULL) {
        /* Should never happen */
        PyErr_SetString(PyExc_SystemError, "no thread state obtained");
        return NULL;
    }

    obj = (SuperLUGlobalObject*)PyDict_GetItemString(thread_dict, key);
    if (obj && Py_TYPE(obj) == &SuperLUGlobalType) {
        return obj;
    }

    obj = (SuperLUGlobalObject*)PyObject_New(SuperLUGlobalObject, &SuperLUGlobalType);
    if (obj == NULL) {
        return (SuperLUGlobalObject*)PyErr_NoMemory();
    }
    obj->memory_dict = PyDict_New();
    obj->jmpbuf_valid = 0;
    if (obj->memory_dict == NULL) {
        Py_DECREF(obj);
        return NULL;
    }

    if (PyDict_SetItemString(thread_dict, key, (PyObject *)obj) < 0) {
        Py_DECREF(obj);
        return NULL;
    }

    /*
        Py_DECREF is added because get_tls_global returns
        borrowed reference. This avoids memory leak of obj.
    */
    Py_DECREF(obj);
    return obj;
}

jmp_buf *superlu_python_jmpbuf(void)
{
    SuperLUGlobalObject *g;

    g = get_tls_global();
    if (g == NULL) {
        abort();
    }
    g->jmpbuf_valid = 1;
    return &g->jmpbuf;
}

/* Install `tracker` as the thread's allocation tracker, and hand the previously
 * installed one back to the caller.  Steals a reference to `tracker` and returns
 * a new reference to the replaced tracker.
 *
 * `g` is resolved by the caller, so this is a plain pointer exchange that cannot
 * fail --- which is what lets the callers swap back on their error paths without
 * having to invent a recovery strategy for "restoring the tracker failed".
 */
PyObject *superlu_swap_memory_tracker(SuperLUGlobalObject *g, PyObject *tracker)
{
    PyObject *replaced_tracker;

    replaced_tracker = g->memory_dict;
    g->memory_dict = tracker;
    return replaced_tracker;
}

void superlu_python_module_abort(char *msg)
{
    SuperLUGlobalObject *g;
    NPY_ALLOW_C_API_DEF;

    NPY_ALLOW_C_API;
    g = get_tls_global();
    if (g == NULL) {
        /* We have to longjmp (or SEGV results), but the
           destination is not known --- no choice but abort.
           However, this should never happen.
        */
        abort();
    }
    PyErr_SetString(PyExc_RuntimeError, msg);

    if (!g->jmpbuf_valid) {
        abort();
    }

    g->jmpbuf_valid = 0;
    NPY_DISABLE_C_API;

    longjmp(g->jmpbuf, -1);
}

void *superlu_python_module_malloc(size_t size)
{
    SuperLUGlobalObject *g;
    PyObject *key = NULL;
    void *mem_ptr;
    NPY_ALLOW_C_API_DEF;

    NPY_ALLOW_C_API;
    g = get_tls_global();
    if (g == NULL) {
        return NULL;
    }
    mem_ptr = malloc(size);
    if (mem_ptr == NULL) {
        NPY_DISABLE_C_API;
        return NULL;
    }
    key = PyLong_FromVoidPtr(mem_ptr);
    if (key == NULL)
        goto fail;
    if (PyDict_SetItem(g->memory_dict, key, Py_None))
        goto fail;
    Py_DECREF(key);
    NPY_DISABLE_C_API;

    return mem_ptr;

  fail:
    Py_XDECREF(key);
    NPY_DISABLE_C_API;
    free(mem_ptr);
    superlu_python_module_abort
        ("superlu_malloc: Cannot set dictionary key value in malloc.");
    return NULL;

}

void superlu_python_module_free(void *ptr)
{
    SuperLUGlobalObject *g;
    PyObject *key;
    PyObject *exc;
    NPY_ALLOW_C_API_DEF;

    if (ptr == NULL)
        return;

    NPY_ALLOW_C_API;
    g = get_tls_global();
    if (g == NULL) {
        abort();
    }
    exc = PyErr_GetRaisedException();
    key = PyLong_FromVoidPtr(ptr);
    /* This will only free the pointer if it could find it in the dictionary
     * of already allocated pointers --- thus after abort, the module can free all
     * the memory that "might" have been allocated to avoid memory leaks on abort
     * calls. If the key cannot be created, leave the pointer in the dictionary.
     * The dictionary owner frees residual tracked allocations during deallocation.
     */
    if (key != NULL) {
        if (!PyDict_DelItem(g->memory_dict, key)) {
            free(ptr);
        }
        Py_DECREF(key);
    }
    PyErr_SetRaisedException(exc);
    NPY_DISABLE_C_API;
}


/* Free every allocation still registered in `tracker` and empty it.
 *
 * The keys are the pointers themselves (the values are unused), and every
 * pointer in a tracker came from `superlu_python_module_malloc`, so freeing all
 * of them is safe.  Callers must only invoke this on a tracker whose contents
 * they own:
 *
 *  - `Py_gstrf` gives each factorization a tracker of its own, so once it
 *    succeeds that tracker holds exactly the allocations reachable from the
 *    factor's L, U, perm_r and perm_c.  Freeing all of them in
 *    `SuperLU_dealloc` is therefore equivalent to the SUPERLU_FREE /
 *    Destroy_SuperNode_Matrix / Destroy_CompCol_Matrix sequence it replaces,
 *    and additionally reclaims anything gstrf left behind that is not
 *    reachable from L or U.
 *  - `SuperLUGlobal_dealloc` frees what is left in the thread's own tracker,
 *    which after the above can only be residue from an aborted call.
 */
void superlu_free_tracked_allocations(PyObject *tracker)
{
    PyObject *key, *value;
    Py_ssize_t pos = 0;

    while (PyDict_Next(tracker, &pos, &key, &value)) {
        void *ptr;
        ptr = PyLong_AsVoidPtr(key);
        free(ptr);
    }
    PyDict_Clear(tracker);
}


static void SuperLUGlobal_dealloc(SuperLUGlobalObject *self)
{
    if (self->memory_dict != NULL) {
        superlu_free_tracked_allocations(self->memory_dict);
    }

    Py_XDECREF(self->memory_dict);
    PyObject_Del(self);
}


PyTypeObject SuperLUGlobalType = {
    PyVarObject_HEAD_INIT(NULL, 0)
    .tp_name = "_SuperLUGlobal",
    .tp_basicsize = sizeof(SuperLUGlobalObject),
    .tp_dealloc = (destructor)SuperLUGlobal_dealloc,
    .tp_flags = Py_TPFLAGS_DEFAULT,
};


/*
 * Stub for error handling; does nothing, as we don't want to spew debug output.
 */

int input_error(char *srname, int *info)
{
    return 0;
}

/*
 * Stubs for Harwell Subroutine Library functions that SuperLU tries to call.
 */

void mc64id_(int *a)
{
    superlu_python_module_abort("chosen functionality not available");
}

void mc64ad_(int *a, int *b, int *c, int d[], int e[], double f[],
             int *g, int h[], int *i, int j[], int *k, double l[],
             int m[], int n[])
{
    superlu_python_module_abort("chosen functionality not available");
}
