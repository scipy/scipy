#include "Python.h"
#include "structmember.h"
#include "numpy/arrayobject.h"

enum OpCodes {
    OP_COPY = 0,
    OP_NEG,
    OP_ADD,
    OP_SUB,
    OP_MUL,
    OP_DIV,
    OP_ADD_C,
    OP_SUB_C,
    OP_MUL_C,
    OP_DIV_C,
};

#define BLOCK_SIZE1 128
#define BLOCK_SIZE2 8

typedef struct
{
    PyObject_HEAD
    unsigned int n_inputs;
    unsigned int n_temps;
    PyObject *program;      /* a python string */
    PyObject *constants;    /* a PyArrayObject */
    PyObject *input_names;  /* tuple of strings */
    double **mem;
    double *temps;
} NumExprObject;

static void
NumExpr_dealloc(NumExprObject *self)
{
    Py_XDECREF(self->program);
    Py_XDECREF(self->constants);
    Py_XDECREF(self->input_names);
    PyMem_Del(self->mem);
    PyMem_Del(self->temps);
}

static PyObject *
NumExpr_new(PyTypeObject *type, PyObject *args, PyObject *kwds)
{
    NumExprObject *self;

    self = (NumExprObject *)type->tp_alloc(type, 0);
    if (self != NULL) {
        self->n_inputs = 0;
        self->n_temps = 0;
        self->mem = NULL;
        self->temps = NULL;
        self->program = PyString_FromString("");
        if (!self->program) {
            Py_DECREF(self);
            return NULL;
        }
        intp dims[] = {0};
        self->constants = PyArray_SimpleNew(1, dims, PyArray_DOUBLE);
        if (!self->constants) {
            Py_DECREF(self);
            return NULL;
        }
        Py_INCREF(Py_None);
        self->input_names = Py_None;
    }
    return (PyObject *)self;
}

static int
NumExpr_init(NumExprObject *self, PyObject *args, PyObject *kwds)
{
    PyObject *program = NULL, *o_constants = NULL, *input_names = NULL, *tmp;
    static char *kwlist[] = {"n_inputs", "n_temps",
                             "program", "constants",
                             "input_names", NULL};

    if (!PyArg_ParseTupleAndKeywords(args, kwds, "IIS|OO", kwlist,
                                     &self->n_inputs,
                                     &self->n_temps,
                                     &program, &o_constants, &input_names)) {
        return -1;
    }
    if (self->n_inputs + self->n_temps + 1 > 256) {
        PyErr_SetString(PyExc_ValueError,
                        "number inputs + outputs + temporaries must be <= 256");
        return -1;
    }
    if (o_constants) {
        PyObject *constants;
        constants = PyArray_FromAny(o_constants,
                                    PyArray_DescrFromType(PyArray_DOUBLE),
                                    1, 1,
                                    ENSURECOPY | CARRAY_FLAGS_RO,
                                    NULL);
        if (!constants) {
            return -1;
        }
        if (PyArray_DIM(constants, 0) > 256) {
            PyErr_SetString(PyExc_ValueError,
                            "number of constants must be <= 256");
            Py_DECREF(constants);
            return -1;
        }
        tmp = self->constants;
        self->constants = constants;
        Py_XDECREF(tmp);
    }

    tmp = self->program;
    Py_INCREF(program);
    self->program = program;
    Py_XDECREF(tmp);

    tmp = self->input_names;
    if (!input_names) {
        input_names = Py_None;
    }
    Py_INCREF(input_names);
    self->input_names = input_names;
    Py_XDECREF(tmp);

    PyMem_Del(self->mem);
    self->mem = PyMem_New(double *, 1 + self->n_inputs + self->n_temps);
    PyMem_Del(self->temps);
    self->temps = PyMem_New(double, self->n_temps * BLOCK_SIZE1);

    return 0;
}

static PyMemberDef NumExpr_members[] = {
    {"n_inputs", T_UINT, offsetof(NumExprObject, n_inputs), READONLY, NULL},
    {"n_temps", T_UINT, offsetof(NumExprObject, n_temps), READONLY, NULL},
    {"program", T_OBJECT_EX, offsetof(NumExprObject, program), READONLY, NULL},
    {"constants", T_OBJECT_EX, offsetof(NumExprObject, constants),
     READONLY, NULL},
    {"input_names", T_OBJECT, offsetof(NumExprObject, input_names), 0, NULL},
    {NULL},
};

static int
run_interpreter(NumExprObject *self, int len, double *output, double **inputs)
{
    double **mem, *constants;
    char *program;
    unsigned int n_inputs, prog_len, t, blen1, index;

    n_inputs = self->n_inputs;
    mem = self->mem;
    for (t = 0; t < self->n_temps; t++) {
        mem[1+n_inputs+t] = self->temps + BLOCK_SIZE1 * t;
    }

    if (PyString_AsStringAndSize(self->program, &program, &prog_len) < 0) {
        return -1;
    }
    constants = PyArray_DATA(self->constants);

    blen1 = len - len % BLOCK_SIZE1;
    for (index = 0; index < blen1; index += BLOCK_SIZE1) {
#define VECTOR_SIZE BLOCK_SIZE1
#include "interp_body.c"
#undef VECTOR_SIZE
    }
    if (len != blen1) {
        int blen2 = len - len % BLOCK_SIZE2;
        for (index = blen1; index < blen2; index += BLOCK_SIZE2) {
#define VECTOR_SIZE BLOCK_SIZE2
#include "interp_body.c"
#undef VECTOR_SIZE
        }
        if (len != blen2) {
            int rest = len - blen2;
            index = blen2;
#define VECTOR_SIZE rest
#include "interp_body.c"
#undef VECTOR_SIZE
        }
    }
    return 0;
}

/* keyword arguments are ignored! */
static PyObject *
NumExpr_run(NumExprObject *self, PyObject *args, PyObject *kwds)
{
    PyObject *output = NULL, *a_inputs;
    unsigned int n_inputs;
    int i, len = -1;
    double **inputs;

    n_inputs = PyTuple_Size(args);
    if (self->n_inputs != n_inputs) {
        return PyErr_Format(PyExc_ValueError,
                            "number of inputs doesn't match program");
    }
    a_inputs = PyTuple_New(n_inputs);
    if (!a_inputs) {
        return NULL;
    }
    inputs = PyMem_New(double *, n_inputs);
    for (i = 0; i < n_inputs; i++) {
        PyObject *o = PyTuple_GetItem(args, i); /* borrowed ref */
        PyObject *a = PyArray_ContiguousFromAny(o, PyArray_DOUBLE, 0, 0);
        PyTuple_SET_ITEM(a_inputs, i, a);  /* steals reference */
        if (len == -1) {
            len = PyArray_SIZE(a);
            output = PyArray_SimpleNew(PyArray_NDIM(a),
                                       PyArray_DIMS(a),
                                       PyArray_DOUBLE);
        } else {
            if (len != PyArray_SIZE(a)) {
                Py_XDECREF(a_inputs);
                Py_XDECREF(output);
                PyMem_Del(inputs);
                return PyErr_Format(PyExc_ValueError,
                                    "all inputs must be the same size");
            }
        }
        inputs[i] = PyArray_DATA(a);
    }
    if (run_interpreter(self, len, PyArray_DATA(output), inputs) < 0) {
        Py_XDECREF(output);
        output = NULL;
        PyErr_SetString(PyExc_RuntimeError,
                        "an error occurred while running the program");
    }
    Py_XDECREF(a_inputs);
    PyMem_Del(inputs);
    return output;
}

static PyMethodDef NumExpr_methods[] = {
    {"run", (PyCFunction) NumExpr_run, METH_VARARGS, NULL},
    {NULL, NULL}
};

static PyTypeObject NumExprType = {
    PyObject_HEAD_INIT(NULL)
    0,                         /*ob_size*/
    "numexpr.NumExpr",         /*tp_name*/
    sizeof(NumExprObject),     /*tp_basicsize*/
    0,                         /*tp_itemsize*/
    (destructor)NumExpr_dealloc, /*tp_dealloc*/
    0,                         /*tp_print*/
    0,                         /*tp_getattr*/
    0,                         /*tp_setattr*/
    0,                         /*tp_compare*/
    0,                         /*tp_repr*/
    0,                         /*tp_as_number*/
    0,                         /*tp_as_sequence*/
    0,                         /*tp_as_mapping*/
    0,                         /*tp_hash */
    (ternaryfunc)NumExpr_run,  /*tp_call*/
    0,                         /*tp_str*/
    0,                         /*tp_getattro*/
    0,                         /*tp_setattro*/
    0,                         /*tp_as_buffer*/
    Py_TPFLAGS_DEFAULT | Py_TPFLAGS_BASETYPE, /*tp_flags*/
    "NumExpr objects",         /* tp_doc */
    0,		               /* tp_traverse */
    0,		               /* tp_clear */
    0,		               /* tp_richcompare */
    0,		               /* tp_weaklistoffset */
    0,		               /* tp_iter */
    0,		               /* tp_iternext */
    NumExpr_methods,           /* tp_methods */
    NumExpr_members,           /* tp_members */
    0,                         /* tp_getset */
    0,                         /* tp_base */
    0,                         /* tp_dict */
    0,                         /* tp_descr_get */
    0,                         /* tp_descr_set */
    0,                         /* tp_dictoffset */
    (initproc)NumExpr_init,    /* tp_init */
    0,                         /* tp_alloc */
    NumExpr_new,               /* tp_new */
};

static PyMethodDef module_methods[] = {
    {NULL}
};

void
initinterpreter(void)
{
    PyObject *m, *d, *o;
    int r;

    if (PyType_Ready(&NumExprType) < 0)
        return;

    m = Py_InitModule3("interpreter", module_methods, NULL);
    if (m == NULL)
        return;

    Py_INCREF(&NumExprType);
    PyModule_AddObject(m, "NumExpr", (PyObject *)&NumExprType);

    import_array();

    d = PyDict_New();
    if (!d) return;

#define add_op(sname, name) o = PyInt_FromLong(name);   \
    r = PyDict_SetItemString(d, sname, o);              \
    Py_XDECREF(o);                                      \
    if (r < 0) {PyErr_SetString(PyExc_RuntimeError, "add_op"); return;}

    add_op("copy", OP_COPY);
    add_op("neg", OP_NEG);
    add_op("add", OP_ADD);
    add_op("sub", OP_SUB);
    add_op("mul", OP_MUL);
    add_op("div", OP_DIV);
    add_op("add_c", OP_ADD_C);
    add_op("sub_c", OP_SUB_C);
    add_op("mul_c", OP_MUL_C);
    add_op("div_c", OP_DIV_C);
#undef add_op

    if (PyModule_AddObject(m, "opcodes", d) < 0) return;

}
