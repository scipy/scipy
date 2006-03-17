#include "Python.h"
#include "structmember.h"
#include "numpy/arrayobject.h"
#include "math.h"
#ifdef _WIN32
#define inline __inline
#endif


enum OpCodes {
    OP_NOOP = 0,
    OP_COPY,
    OP_ONES_LIKE,
    OP_NEG,
    OP_ADD,
    OP_SUB,
    OP_MUL,
    OP_DIV,
    OP_POW,
    OP_MOD,
    OP_GT,
    OP_GE,
    OP_EQ,
    OP_NE,
    OP_SIN,
    OP_COS,
    OP_TAN,
    OP_SQRT,
    OP_ARCTAN2,
    OP_WHERE,
    OP_FUNC_1,
    OP_FUNC_2,

    OP_LAST_OP
};

/*
   Lots of functions still to be added: exp, ln, log10, etc, etc. Still not
   sure which get there own opcodes and which get relegated to loopup table.
   Some functions at least (sin and arctan2 for instance) seem to have a large
   slowdown when run through lookup table. Not entirely sure why.

   To add a function to the lookup table, add to FUNC_CODES (first
   group is 1-arg functions, second is 2-arg functions), also to
   functions_1 or functions_2 as appropriate. Finally, use add_func
   down below to add to funccodes. Functions with more arguments
   aren't implemented at present, but should be easy; just copy the 1-
   or 2-arg case.

   To add a function opcode, just copy OP_SIN or OP_ARCTAN2.

*/

enum Func1Codes {
    FUNC_SINH = 0,
    FUNC_COSH,
    FUNC_TANH,

    FUNC1_LAST
};

enum Func2Codes {
    FUNC_FMOD = 0,

    FUNC2_LAST
};

typedef double (*Func1Ptr)(double);

Func1Ptr functions_1[] = {
    sinh,
    cosh,
    tanh,
};

typedef double (*Func2Ptr)(double, double);

Func2Ptr functions_2[] = {
    fmod,
};


#define BLOCK_SIZE1 128
#define BLOCK_SIZE2 8

typedef struct
{
    PyObject_HEAD
    unsigned int n_inputs;
    unsigned int n_temps;
    unsigned int r_constants;
    unsigned int r_temps;
    unsigned int r_end;
    PyObject *program;      /* a python string */
    PyObject *constants;    /* a PyArrayObject */
    PyObject *input_names;  /* tuple of strings */
    double **mem;
    PyObject *temps;        /* a PyArrayObject */
} NumExprObject;

static void
NumExpr_dealloc(NumExprObject *self)
{
    Py_XDECREF(self->program);
    Py_XDECREF(self->constants);
    Py_XDECREF(self->input_names);
    Py_XDECREF(self->temps);
    PyMem_Del(self->mem);
}

static PyObject *
NumExpr_new(PyTypeObject *type, PyObject *args, PyObject *kwds)
{
    NumExprObject *self;
    intp dims[] = {0}, dims2[] = {0,0};
    self = (NumExprObject *)type->tp_alloc(type, 0);
    if (self != NULL) {
        self->n_inputs = 0;
        self->n_temps = 0;
        self->r_constants = 0;
        self->r_temps = 0;
        self->r_end = 0;
        self->mem = NULL;
        self->program = PyString_FromString("");
        if (!self->program) {
            Py_DECREF(self);
            return NULL;
        }
        self->constants = PyArray_SimpleNew(1, dims, PyArray_DOUBLE);
        if (!self->constants) {
            Py_DECREF(self);
            return NULL;
        }
        self->temps = PyArray_SimpleNew(2, dims2, PyArray_DOUBLE);
        if (!self->temps) {
            Py_DECREF(self);
            return NULL;
        }
        Py_INCREF(Py_None);
        self->input_names = Py_None;
    }
    return (PyObject *)self;
}

static int
check_program(NumExprObject *self)
{
    unsigned char *program;
    int prog_len, pc;

    if (PyString_AsStringAndSize(self->program, (char **)&program,
                                 &prog_len) < 0) {
        return -1;
    }
    if (prog_len % 4 != 0) {
        return -1;
    }
    for (pc = 0; pc < prog_len; pc += 4) {
        unsigned int op = program[pc];
        unsigned int store = program[pc+1];
        unsigned int arg1 = program[pc+2];
        unsigned int arg2 = program[pc+3];
        if (op >= OP_LAST_OP) return -1;
        if (op == OP_NOOP) {
            continue;
        }
        if (store >= self->r_end) return -1;
        if (arg1 >= self->r_end && arg1 != 255) return -1;
        if (op == OP_WHERE) {
            unsigned int arg3;
            if (pc + 5 > prog_len) return -1;
            if (program[pc+4] != OP_NOOP) return -1;
            arg3 = program[pc+5];
            if (arg2 >= self->r_end) return -1;
        }
        if (op == OP_FUNC_1) {
            if (arg2 >= FUNC1_LAST) return -1;
        } else if (arg2 >= self->r_end && arg2 != 255) return -1;
        if (op == OP_FUNC_2) {
            unsigned int arg3;
            if (pc + 5 > prog_len) return -1;
            if (program[pc+4] != OP_NOOP) return -1;
            arg3 = program[pc+5];
            if (arg3 >= FUNC2_LAST) return -1;
        }
    }
    return 0;
}

static int
NumExpr_init(NumExprObject *self, PyObject *args, PyObject *kwds)
{
    int i, j, n_constants;
    PyObject *program = NULL, *o_constants = NULL, *input_names = NULL;
    PyObject *constants, *tmp;
    static char *kwlist[] = {"n_inputs", "n_temps",
                             "program", "constants",
                             "input_names", NULL};

    if (!PyArg_ParseTupleAndKeywords(args, kwds, "IIS|OO", kwlist,
                                     &self->n_inputs,
                                     &self->n_temps,
                                     &program, &o_constants, &input_names)) {
        return -1;
    }
    if (o_constants) {
        constants = PyArray_FromAny(o_constants,
                                    PyArray_DescrFromType(PyArray_DOUBLE),
                                    1, 1,
                                    ENSURECOPY | CARRAY_FLAGS_RO,
                                    NULL);
        if (!constants) {
            return -1;
        }
    } else {
        intp dims[] = {0};
        constants = PyArray_SimpleNew(1, dims, PyArray_DOUBLE);
        if (!constants) {
            return -1;
        }
    }
    n_constants = PyArray_DIM(constants, 0);
    self->r_constants = 1 + self->n_inputs;
    self->r_temps = self->r_constants + n_constants;
    self->r_end = self->r_temps + self->n_temps;

    if (self->r_end > 255) {
        PyErr_SetString(PyExc_ValueError,
                        "expression too complicated");
        Py_DECREF(constants);
        return -1;
    }

    tmp = self->constants;
    self->constants = constants;
    Py_XDECREF(tmp);

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
    self->mem = PyMem_New(double *, self->r_end);

    {
        intp dims[2];
        double *tdata;

        tmp = self->temps;
        dims[0] = n_constants + self->n_temps;
        dims[1] = BLOCK_SIZE1;
        self->temps = PyArray_SimpleNew(2, dims, PyArray_DOUBLE);
        Py_XDECREF(tmp);

        PyArray_FILLWBYTE(self->temps, 0);
        tdata = PyArray_DATA(self->temps);
        for (i = 0; i < n_constants; i++ ) {
            double c = ((double *)PyArray_DATA(self->constants))[i];
            double *p = tdata + i * BLOCK_SIZE1;
            for (j = 0; j < BLOCK_SIZE1; j++) {
                p[j] = c;
            }
        }
    }

    if (check_program(self) < 0) {
        PyErr_SetString(PyExc_ValueError, "invalid program");
        return -1;
    }

    return 0;
}

static PyMemberDef NumExpr_members[] = {
    {"n_inputs", T_UINT, offsetof(NumExprObject, n_inputs), READONLY, NULL},
    {"n_temps", T_UINT, offsetof(NumExprObject, n_temps), READONLY, NULL},
    {"r_constants", T_UINT, offsetof(NumExprObject, r_constants),
     READONLY, NULL},
    {"r_temps", T_UINT, offsetof(NumExprObject, r_temps), READONLY, NULL},
    {"r_end", T_UINT, offsetof(NumExprObject, r_end), READONLY, NULL},
    {"program", T_OBJECT_EX, offsetof(NumExprObject, program), READONLY, NULL},
    {"constants", T_OBJECT_EX, offsetof(NumExprObject, constants),
     READONLY, NULL},
    {"temps", T_OBJECT_EX, offsetof(NumExprObject, temps), READONLY, NULL},
    {"input_names", T_OBJECT, offsetof(NumExprObject, input_names), 0, NULL},
    {NULL},
};

struct vm_params {
    int prog_len;
    unsigned char *program;
    unsigned int n_inputs;
    unsigned int r_end;
    double *output;
    double **inputs;
    double **mem;
};

#if DO_BOUNDS_CHECK
#define BOUNDS_CHECK(arg) if (((arg) >= params.r_end) && ((arg) != 255)) { \
        *pc_error = pc;                                                 \
        return -2;                                                      \
    }
#else
#define BOUNDS_CHECK(arg)
#endif

static inline int
vm_engine_1(int start, int blen, struct vm_params params, int *pc_error)
{
    unsigned int index;
    for (index = start; index < blen; index += BLOCK_SIZE1) {
#define UNROLL 0
#define VECTOR_SIZE BLOCK_SIZE1
#include "interp_body.c"
#undef VECTOR_SIZE
#undef UNROLL
    }
    return 0;
}

static inline int
vm_engine_2(int start, int blen, struct vm_params params, int *pc_error)
{
    unsigned int index;
    for (index = start; index < blen; index += BLOCK_SIZE2) {
#define UNROLL 0
#define VECTOR_SIZE BLOCK_SIZE2
#include "interp_body.c"
#undef VECTOR_SIZE
#undef UNROLL
    }
    return 0;
}

static inline int
vm_engine_rest(int start, int blen, struct vm_params params, int *pc_error)
{
    unsigned int index = start;
    unsigned int rest = blen - start;
#define UNROLL 0
#define VECTOR_SIZE rest
#include "interp_body.c"
#undef VECTOR_SIZE
#undef UNROLL
    return 0;
}

static int
run_interpreter(NumExprObject *self, int len, double *output, double **inputs,
                int *pc_error)
{
    int r;
    unsigned int t, blen1, blen2, n_constants;
    struct vm_params params;
    double *tdata;

    *pc_error = -1;
    if (PyString_AsStringAndSize(self->program, (char **)&(params.program),
                                 &(params.prog_len)) < 0) {
        return -1;
    }
    params.n_inputs = self->n_inputs;
    params.output = output;
    params.inputs = inputs;
    params.mem = self->mem;
    n_constants = PyArray_DIM(self->constants, 0);
    tdata = PyArray_DATA(self->temps);
    for (t = 0; t < n_constants + self->n_temps; t++) {
        params.mem[self->r_constants + t] = tdata + BLOCK_SIZE1 * t;
    }
    params.r_end = self->r_end;

    blen1 = len - len % BLOCK_SIZE1;
    r = vm_engine_1(0, blen1, params, pc_error);
    if (r < 0) return r;
    if (len != blen1) {
        blen2 = len - len % BLOCK_SIZE2;
        r = vm_engine_2(blen1, blen2, params, pc_error);
        if (r < 0) return r;
        if (len != blen2) {
            r = vm_engine_rest(blen2, len, params, pc_error);
            if (r < 0) return r;
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
    int i, len = -1, r, pc_error;
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
    r = run_interpreter(self, len, PyArray_DATA(output), inputs, &pc_error);
    if (r < 0) {
        Py_XDECREF(output);
        output = NULL;
        if (r == -1) {
            PyErr_SetString(PyExc_RuntimeError,
                            "an error occurred while running the program");
        } else if (r == -2) {
            PyErr_Format(PyExc_RuntimeError,
                         "bad argument at pc=%d", pc_error);
        } else if (r == -3) {
            PyErr_Format(PyExc_RuntimeError,
                         "bad opcode at pc=%d", pc_error);
        } else {
            PyErr_SetString(PyExc_RuntimeError,
                            "unknown error occurred while running the program");
        }
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

    add_op("noop", OP_NOOP);
    add_op("copy", OP_COPY);
    add_op("ones_like", OP_ONES_LIKE);
    add_op("neg", OP_NEG);
    add_op("add", OP_ADD);
    add_op("sub", OP_SUB);
    add_op("mul", OP_MUL);
    add_op("div", OP_DIV);
    add_op("pow", OP_POW);
    add_op("mod", OP_MOD);
    add_op("gt", OP_GT);
    add_op("ge", OP_GE);
    add_op("eq", OP_EQ);
    add_op("ne", OP_NE);
    add_op("sin", OP_SIN);
    add_op("cos", OP_COS);
    add_op("tan", OP_TAN);
    add_op("sqrt", OP_SQRT);
    add_op("arctan2", OP_ARCTAN2);
    add_op("where", OP_WHERE);
    add_op("func_1", OP_FUNC_1);
    add_op("func_2", OP_FUNC_2);
#undef add_op

    if (PyModule_AddObject(m, "opcodes", d) < 0) return;

    d = PyDict_New();
    if (!d) return;

#define add_func(sname, name) o = PyInt_FromLong(name); \
    r = PyDict_SetItemString(d, sname, o);              \
    Py_XDECREF(o);                                      \
    if (r < 0) {PyErr_SetString(PyExc_RuntimeError, "add_func"); return;}

    add_func("sinh", FUNC_SINH);
    add_func("cosh", FUNC_COSH);
    add_func("tanh", FUNC_TANH);

    add_func("fmod", FUNC_FMOD);

#undef add_func

   if (PyModule_AddObject(m, "funccodes", d) < 0) return;

}
