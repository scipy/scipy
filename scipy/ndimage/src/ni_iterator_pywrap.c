#include "nd_image.h"
#include "ni_iterators.h"


PyObject*
PyUniformFilter1D(PyObject *self, PyObject *args, PyObject *kwds)
{
    PyArrayObject *input = NULL;
    PyArrayObject *output = NULL;
    int axis = NPY_MAXDIMS;
    npy_intp filter_size = 1;
    npy_intp filter_offset = 0;
    NI_ExtendMode extend_mode = NI_EXTEND_NEAREST;
    npy_double extend_value = 0.0;
    LineBufferIterator *lbiter = NULL;


    static char *kwlist[] = {"input", "output", "axis", "filter_size",
                             "filter_offset", "extend_mode", "extend_value",
                             NULL};

    if (!PyArg_ParseTupleAndKeywords(args, kwds, "O!O!in|nid", kwlist,
                                     &PyArray_Type, &input,
                                     &PyArray_Type, &output,
                                     &axis, &filter_size, &filter_offset,
                                     &extend_mode, &extend_value)) {
        goto fail;
    }

    lbiter = LBI_New(input, axis, output, filter_size, filter_offset,
                     extend_mode, extend_value);
    if (lbiter == NULL) {
        goto fail;
    }

    do {
        const npy_double *read_front = LBI_GetInputBuffer(lbiter);
        const npy_double *read_back = read_front;
        npy_double *write = LBI_GetOutputBuffer(lbiter);
        npy_double running_sum = 0.0;
        npy_intp filter_size_copy = filter_size;
        npy_intp line_length = LBI_GetLineLength(lbiter);

        while (filter_size_copy--) {
            running_sum += *read_front++;
        }
        *write++ = running_sum / filter_size;
        while (--line_length) {
            running_sum += *read_front++ - *read_back++;
            *write++ = running_sum / filter_size;
        }
    } while(LBI_Next(lbiter));

    LBI_Delete(lbiter);

    Py_INCREF(output);
    return (PyObject *)output;

fail:
    return (PyObject *)LBI_Delete(lbiter);
}


typedef struct {
    PyObject_HEAD
    LineBufferIterator *lbiter;
    PyObject *in_buffer;
    PyObject *out_buffer;
    int iter_started;
} PyLineBufferIterator;


static PyObject*
pylbi_new(PyTypeObject *type, PyObject *args, PyObject *kwds)
{
    PyLineBufferIterator *self;

    self = (PyLineBufferIterator *)type->tp_alloc(type, 0);
    if (self != NULL) {
        self->lbiter = NULL;
        self->in_buffer = NULL;
        self->out_buffer = NULL;
        self->iter_started = 0;
    }

    return (PyObject *)self;
}


static int
pylbi_init(PyLineBufferIterator *self, PyObject *args, PyObject *kwds)
{
    PyArrayObject *input = NULL;
    PyArrayObject *output = NULL;
    int axis = NPY_MAXDIMS;
    npy_intp filter_size = 1;
    npy_intp filter_offset = 0;
    NI_ExtendMode extend_mode = NI_EXTEND_NEAREST;
    npy_double extend_value = 0.0;
    npy_intp buffer_len;
    static char *kwlist[] = {"input", "output", "axis", "filter_size",
                             "filter_offset", "extend_mode", "extend_value",
                             NULL};

    if (!PyArg_ParseTupleAndKeywords(args, kwds, "O!O!in|nid", kwlist,
                                     &PyArray_Type, &input,
                                     &PyArray_Type, &output,
                                     &axis, &filter_size, &filter_offset,
                                     &extend_mode, &extend_value)) {
        goto fail;
    }

    LBI_Delete(self->lbiter);
    self->lbiter = LBI_New(input, axis, output, filter_size, filter_offset,
                           extend_mode, extend_value);
    if (self->lbiter == NULL) {
        goto fail;
    }

    Py_XDECREF(self->out_buffer);
    buffer_len = LBI_GetLineLength(self->lbiter);
    self->out_buffer = PyArray_SimpleNewFromData(1, &buffer_len, NPY_DOUBLE,
                                            LBI_GetOutputBuffer(self->lbiter));
    if (self->out_buffer == NULL) {
        goto fail;
    }

    Py_XDECREF(self->in_buffer);
    buffer_len += filter_size - 1;
    self->in_buffer = PyArray_SimpleNewFromData(1, &buffer_len, NPY_DOUBLE,
                                            LBI_GetInputBuffer(self->lbiter));
    if (self->in_buffer == NULL) {
        goto fail;
    }

    self->iter_started = 0;

    return 0;

fail:
    LBI_Delete(self->lbiter);
    Py_XDECREF(self->in_buffer);
    Py_XDECREF(self->out_buffer);
    return -1;
}


PyObject*
pylbi_getiter(PyLineBufferIterator *self)
{
    Py_INCREF((PyObject *)self);
    return (PyObject *)self;
}


PyObject*
pylbi_next(PyLineBufferIterator *self)
{
    if (!self->iter_started || LBI_Next(self->lbiter)) {
        self->iter_started = 1;
        return Py_BuildValue("OO", self->in_buffer, self->out_buffer);
    }
    return NULL;
}


static void
pylbi_dealloc(PyLineBufferIterator* self)
{
    LBI_Delete(self->lbiter);
    Py_XDECREF(self->in_buffer);
    Py_XDECREF(self->out_buffer);
    Py_TYPE(self)->tp_free((PyObject*)self);
}


NPY_NO_EXPORT PyTypeObject PyLineBufferIterator_Type = {
    PyVarObject_HEAD_INIT(NULL, 0)
    "ndimage.LineBufferIterator",               /* tp_name */
    sizeof(PyLineBufferIterator),               /* tp_basicsize */
    0,                                          /* tp_itemsize */
    /* methods */
    (destructor)pylbi_dealloc,                  /* tp_dealloc */
    0,                                          /* tp_print */
    0,                                          /* tp_getattr */
    0,                                          /* tp_setattr */
#if defined(NPY_PY3K)
    0,                                          /* tp_reserved */
#else
    0,                                          /* tp_compare */
#endif
    0,                                          /* tp_repr */
    0,                                          /* tp_as_number */
    0,                                          /* tp_as_sequence */
    0,                                          /* tp_as_mapping */
    0,                                          /* tp_hash */
    0,                                          /* tp_call */
    0,                                          /* tp_str */
    0,                                          /* tp_getattro */
    0,                                          /* tp_setattro */
    0,                                          /* tp_as_buffer */
    Py_TPFLAGS_DEFAULT,                         /* tp_flags */
    0,                                          /* tp_doc */
    0,                                          /* tp_traverse */
    0,                                          /* tp_clear */
    0,                                          /* tp_richcompare */
    0,                                          /* tp_weaklistoffset */
    (getiterfunc)pylbi_getiter,                 /* tp_iter */
    (iternextfunc)pylbi_next,                   /* tp_iternext */
    0,                                          /* tp_methods */
    0,                                          /* tp_members */
    0,                                          /* tp_getset */
    0,                                          /* tp_base */
    0,                                          /* tp_dict */
    0,                                          /* tp_descr_get */
    0,                                          /* tp_descr_set */
    0,                                          /* tp_dictoffset */
    (initproc)pylbi_init,                       /* tp_init */
    0,                                          /* tp_alloc */
    (newfunc)pylbi_new,                         /* tp_new */
 };


static PyMethodDef ni_iterators_methods[] = {

    {"UniformFilter1D", (PyCFunction)PyUniformFilter1D,
     METH_VARARGS | METH_KEYWORDS,
     "LineBufferIterator version of ndimage.uniform_filter1D."},
    {NULL, NULL, 0, NULL}        /* Sentinel */
};


#ifndef PyMODINIT_FUNC /* declarations for DLL import/export */
#define PyMODINIT_FUNC void
#endif
PyMODINIT_FUNC
init_ni_iterators(void)
{
    PyObject* m;

    if (PyType_Ready(&PyLineBufferIterator_Type) < 0) {
        return;
    }

    m = Py_InitModule("_ni_iterators", ni_iterators_methods);
    if (m == NULL) {
        return;
    }

    Py_INCREF(&PyLineBufferIterator_Type);
    PyModule_AddObject(m, "LineBufferIterator",
                      (PyObject *)&PyLineBufferIterator_Type);
    import_array();
}
