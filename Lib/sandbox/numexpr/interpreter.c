#include "Python.h"
#include "structmember.h"
#include "numpy/arrayobject.h"
#include "math.h"

#include "complex_functions.inc"

#ifdef _WIN32
#define inline __inline
#endif

#define BLOCK_SIZE1 128
#define BLOCK_SIZE2 8

/* This file and interp_body should really be generated from a description of
   the opcodes -- there's too much repetition here for manually editing */
   

enum OpCodes {
    OP_NOOP = 0,

    OP_COPY_II,
    OP_ONES_LIKE_II,
    OP_NEG_II,
    OP_ADD_III,
    OP_SUB_III,
    OP_MUL_III,
    OP_DIV_III,
    OP_POW_III,
    OP_MOD_III,
    OP_WHERE_IFII,

    OP_GT_IFF,
    OP_GE_IFF,
    OP_EQ_IFF,
    OP_NE_IFF,
    
    OP_CAST_FI,
    OP_COPY_FF,
    OP_ONES_LIKE_FF,
    OP_NEG_FF,
    OP_ADD_FFF,
    OP_SUB_FFF,
    OP_MUL_FFF,
    OP_DIV_FFF,
    OP_POW_FFF,
    OP_MOD_FFF,
    OP_SIN_FF,
    OP_COS_FF,
    OP_TAN_FF,
    OP_SQRT_FF,
    OP_ARCTAN2_FFF,
    OP_WHERE_FFFF, 
    OP_FUNC_FF,
    OP_FUNC_FFF,
    
    OP_EQ_ICC,
    OP_NE_ICC,
    
    OP_CAST_CI,
    OP_ONES_LIKE_CC,
    OP_CAST_CF,  
    OP_COPY_CC,
    OP_NEG_CC,
    OP_ADD_CCC,
    OP_SUB_CCC,
    OP_MUL_CCC,
    OP_DIV_CCC,
    OP_WHERE_CFCC, 
    OP_FUNC_CC,
    OP_FUNC_CCC,
    
    OP_REAL_FC,
    OP_IMAG_FC,
    OP_COMPLEX_CFF,
    
};

/* returns the sig of the nth op, '\0' if no more ops -1 on failure */
static char op_signature(int op, int n) {
    switch (op) {
        case OP_NOOP:
            break;
        case OP_COPY_II:
        case OP_ONES_LIKE_II:
        case OP_NEG_II:
            if (n == 0 || n == 1) return 'i';
            break;
        case OP_ADD_III:
        case OP_SUB_III:
        case OP_MUL_III:
        case OP_DIV_III:
        case OP_MOD_III:
        case OP_POW_III:
            if (n == 0 || n == 1 || n == 2) return 'i';
            break;
        case OP_WHERE_IFII:
            if (n == 0 || n == 2 || n == 3) return 'i';
            if (n == 1) return 'f';
            break;
        case OP_GT_IFF:
        case OP_GE_IFF:
        case OP_EQ_IFF:
        case OP_NE_IFF:
            if (n == 0) return 'i';
            if (n == 1 || n == 2) return 'f';
            break;
        case OP_CAST_FI:
            if (n == 0) return 'f';
            if (n == 1) return 'i';
            break;
        case OP_COPY_FF:
        case OP_ONES_LIKE_FF:
        case OP_NEG_FF:
        case OP_SIN_FF:
        case OP_COS_FF:
        case OP_TAN_FF:
        case OP_SQRT_FF:
            if (n == 0 || n == 1) return 'f';
            break;
        case OP_ADD_FFF:
        case OP_SUB_FFF:
        case OP_MUL_FFF:
        case OP_DIV_FFF:
        case OP_POW_FFF:
        case OP_MOD_FFF:
        case OP_ARCTAN2_FFF:
            if (n == 0 || n == 1 || n == 2) return 'f';
            break;
        case OP_WHERE_FFFF:
            if (n == 0 || n == 1 || n == 2 || n == 3) return 'f';
            break;
        case OP_FUNC_FF:
            if (n == 0 || n == 1) return 'f';
            if (n == 2) return 'n';
            break;
        case OP_FUNC_FFF:
            if (n == 0 || n == 1 || n == 2) return 'f';
            if (n == 3) return 'n';
            break;
        case OP_EQ_ICC:
        case OP_NE_ICC:
            if (n == 0) return 'i';
            if (n == 1 || n == 2) return 'c';
            break;
        case OP_CAST_CI:
            if (n == 0) return 'c';
            if (n == 1) return 'i';
            break;
        case OP_CAST_CF:
            if (n == 0) return 'c';
            if (n == 1) return 'f';
            break;
        case OP_COPY_CC:
        case OP_ONES_LIKE_CC:
        case OP_NEG_CC:
            if (n == 0 || n == 1) return 'c';
            break;
        case OP_ADD_CCC:
        case OP_SUB_CCC:
        case OP_MUL_CCC:
        case OP_DIV_CCC:
            if (n == 0 || n == 1 || n == 2) return 'c';
            break;
        case OP_WHERE_CFCC:
            if (n == 0 || n == 2 || n == 3) return 'c';
            if (n == 1) return 'f';
            break;
        case OP_FUNC_CC:
            if (n == 0 || n == 1) return 'c';
            if (n == 2) return 'n';
            break;
        case OP_FUNC_CCC:
            if (n == 0 || n == 1 || n == 2) return 'c';
            if (n == 3) return 'n';
            break;
        case OP_REAL_FC:
        case OP_IMAG_FC:
            if (n == 0) return 'f';
            if (n == 1) return 'c';
            break;
        case OP_COMPLEX_CFF:
            if (n == 0) return 'c';
            if (n == 1 || n == 2) return 'f';
            break;
        default:
            return -1;
            break;
    }
    return 0;
}



/*
   Lots of functions still to be added: exp, ln, log10, etc, etc. Still not
   sure which get there own opcodes and which get relegated to loopup table.
   Some functions at least (sin and arctan2 for instance) seem to have a large
   slowdown when run through lookup table. Not entirely sure why.

   To add a function to the lookup table, add to FUNC_CODES (first
   group is 1-arg functions, second is 2-arg functions), also to
   functions_f or functions_ff as appropriate. Finally, use add_func
   down below to add to funccodes. Functions with more arguments
   aren't implemented at present, but should be easy; just copy the 1-
   or 2-arg case.

   To add a function opcode, just copy OP_SIN or OP_ARCTAN2.

*/

enum FuncFFCodes { 
    FUNC_SQRT_FF = 0,
    FUNC_SIN_FF,
    FUNC_COS_FF,
    FUNC_TAN_FF,
    FUNC_ARCSIN_FF,
    FUNC_ARCCOS_FF,
    FUNC_ARCTAN_FF,
    FUNC_SINH_FF,
    FUNC_COSH_FF,
    FUNC_TANH_FF,

    FUNC_FF_LAST
};

typedef double (*FuncFFPtr)(double);

FuncFFPtr functions_f[] = {
    sqrt,
    sin,
    cos, 
    tan,
    asin,
    acos,
    atan,
    sinh,
    cosh,
    tanh,
};

enum FuncFFFCodes {
    FUNC_FMOD_FFF = 0,

    FUNC_FFF_LAST
};

typedef double (*FuncFFFPtr)(double, double);

FuncFFFPtr functions_ff[] = {
    fmod,
};

enum FuncCCCodes {
    FUNC_SQRT_CC = 0,
    FUNC_SIN_CC,
    FUNC_COS_CC,
    FUNC_TAN_CC,
    FUNC_ARCSIN_CC,
    FUNC_ARCCOS_CC,
    FUNC_ARCTAN_CC,
    FUNC_SINH_CC,
    FUNC_COSH_CC,
    FUNC_TANH_CC,
    
    FUNC_CC_LAST 
};


typedef void (*FuncCCPtr)(cdouble*, cdouble*);

FuncCCPtr functions_cc[] = {
    nc_sqrt,
    nc_sin,
    nc_cos, 
    nc_tan,
    nc_asin,
    nc_acos,
    nc_atan,
    nc_sinh,
    nc_cosh,
    nc_tanh,
};

enum FuncCCCCodes {
    FUNC_POW_CCC = 0,
    
    FUNC_CCC_LAST 
};

typedef void (*FuncCCCPtr)(cdouble*, cdouble*, cdouble*);

FuncCCCPtr functions_ccc[] = {
    nc_pow,
};

typedef struct
{
    PyObject_HEAD
    PyObject *signature;    /* a python string */
    PyObject *tempsig;
    PyObject *constsig;
    PyObject *fullsig;
    PyObject *program;      /* a python string */
    PyObject *constants;    /* a tuple of int/float/complex */
    PyObject *input_names;  /* tuple of strings */
    char **mem;             /* pointers to registers */
    char *rawmem;           /* a chunks of raw memory for storing registers */
    int  *memsteps;
    int  rawmemsize;
} NumExprObject;

static void
NumExpr_dealloc(NumExprObject *self)
{
    Py_XDECREF(self->signature);
    Py_XDECREF(self->tempsig);
    Py_XDECREF(self->constsig);
    Py_XDECREF(self->fullsig);
    Py_XDECREF(self->program);
    Py_XDECREF(self->constants);
    Py_XDECREF(self->input_names);
    PyMem_Del(self->mem);
    PyMem_Del(self->rawmem);
    PyMem_Del(self->memsteps);
}

static PyObject *
NumExpr_new(PyTypeObject *type, PyObject *args, PyObject *kwds)
{
    NumExprObject *self;
    intp dims[] = {0}, dims2[] = {0,0};
    self = (NumExprObject *)type->tp_alloc(type, 0);
    if (self != NULL) {
#define INIT_WITH(name, object) \
        self->name = object; \
        if (!self->name) { \
            Py_DECREF(self); \
            return NULL; \
        }

        INIT_WITH(signature, PyString_FromString(""));
        INIT_WITH(tempsig, PyString_FromString(""));
        INIT_WITH(constsig, PyString_FromString(""));
        INIT_WITH(fullsig, PyString_FromString(""));
        INIT_WITH(program, PyString_FromString(""));
        INIT_WITH(constants, PyTuple_New(0));
        Py_INCREF(Py_None);
        self->input_names = Py_None;
        self->mem = NULL;
        self->rawmem = NULL;
        self->memsteps = NULL;
        self->rawmemsize = 0;
#undef INIT_WITH
    }
    return (PyObject *)self;
}

static char
get_return_sig(PyObject* program) {
    char last_opcode, sig;
    int end = PyString_Size(program);
    do {
        end -= 4;
        if (end < 0) return 'X';
    }
    while ((last_opcode = PyString_AS_STRING(program)[end]) == OP_NOOP);
    sig = op_signature(last_opcode, 0);
    if (sig <= 0) return 'X';
    return sig;
}

static intp
size_from_char(char c) 
{
    switch (c) {
        case 'i': return sizeof(long);
        case 'f': return sizeof(double);
        case 'c': return 2*sizeof(double);
        default:
            PyErr_SetString(PyExc_TypeError, "signature value not in 'ifc'");
            return -1;
    }
}

static intp
size_from_sig(PyObject *o)
{
    intp size = 0;
    char *s = PyString_AsString(o);
    if (!s) return -1;
    for (; *s != NULL; s++) {
        intp x = size_from_char(*s);
        if (x == -1) return -1;
        size += x;
    }
    return size;
}

static intp
typecode_from_char(char c)
{
    switch (c) {
        case 'i': return PyArray_LONG;
        case 'f': return PyArray_DOUBLE;
        case 'c': return PyArray_CDOUBLE;
        default:
            PyErr_SetString(PyExc_TypeError, "signature value not in 'ifc'");
            return -1;
    }
}

static int
check_program(NumExprObject *self)
{
    unsigned char *program;
    int prog_len, rno, pc, arg, argloc, argno, n_buffers, n_inputs;
    char sig, *fullsig, *signature;

    if (PyString_AsStringAndSize(self->program, (char **)&program,
                                 &prog_len) < 0) {
        PyErr_Format(PyExc_RuntimeError, "invalid program: can't read program");
        return -1;
    }
    if (prog_len % 4 != 0) {
        PyErr_Format(PyExc_RuntimeError, "invalid program: prog_len % 4 != 0");
        return -1;
    }
    if (PyString_AsStringAndSize(self->fullsig, (char **)&fullsig,
                                 &n_buffers) < 0) {
        PyErr_Format(PyExc_RuntimeError, "invalid program: can't read fullsig");
        return -1;
    }
    if (PyString_AsStringAndSize(self->signature, (char **)&signature,
                                 &n_inputs) < 0) {
        PyErr_Format(PyExc_RuntimeError, "invalid program: can't read signature");
        return -1;
    }
    if (n_buffers > 255) {
        PyErr_Format(PyExc_RuntimeError, "invalid program: too many buffers");
        return -1;
    }
    for (rno = n_inputs+1; rno < n_buffers; rno++) {
        char *bufend = self->mem[rno] + BLOCK_SIZE1 * size_from_char(fullsig[rno]);
        if ( (bufend - self->rawmem) > self->rawmemsize) {
            PyErr_Format(PyExc_RuntimeError, "invalid program: too many buffers");
            return -1;
        }
    }
    for (pc = 0; pc < prog_len; pc += 4) {
        unsigned int op = program[pc];
        if (op == OP_NOOP) {
            continue;
        }
        for (argno = 0; ; argno++) {
            sig = op_signature(op, argno);
            if (sig == -1) {
                PyErr_Format(PyExc_RuntimeError, "invalid program: illegal opcode at %i (%c)", pc, op);
                return -1;
            }
            if (sig == 0) break;
            if (argno < 3) {
                argloc = pc+argno+1;
            }
            if (argno >= 3) {
                if (pc + 1 >= prog_len) {
                    PyErr_Format(PyExc_RuntimeError, "invalid program: double opcode (%c) at end (%i)", pc, sig);
                    return -1;
                }
                argloc = pc+argno+2;
            }
            arg = program[argloc];
            
            if (sig != 'n' && arg >= n_buffers) {
                PyErr_Format(PyExc_RuntimeError, "invalid program: buffer out of range (%i) at %i", arg, argloc);
                return -1;
            }
            if (sig == 'n') {
                if (op == OP_FUNC_FF) {
                    if (arg < 0 || arg >= FUNC_FF_LAST) {
                        PyErr_Format(PyExc_RuntimeError, "invalid program: funccode out of range (%i) at %i", arg, argloc);
                        return -1;
                    }
                } else if (op == OP_FUNC_FFF) {
                    if (arg < 0 || arg >= FUNC_FFF_LAST) {
                        PyErr_Format(PyExc_RuntimeError, "invalid program: funccode out of range (%i) at %i", arg, argloc);
                        return -1;
                    }
                } else if (op == OP_FUNC_CC) {
                    if (arg < 0 || arg >= FUNC_CC_LAST) {
                        PyErr_Format(PyExc_RuntimeError, "invalid program: funccode out of range (%i) at %i", arg, argloc);
                        return -1;
                    }
                } else if (op == OP_FUNC_CCC) {
                    if (arg < 0 || arg >= FUNC_CCC_LAST) {
                        PyErr_Format(PyExc_RuntimeError, "invalid program: funccode out of range (%i) at %i", arg, argloc);
                        return -1;
                    }
                } else {
                    PyErr_Format(PyExc_RuntimeError, "invalid program: internal checker errror processing %i", argloc);
                    return -1;
                }
            } else if (sig != fullsig[arg]) {
                PyErr_Format(PyExc_RuntimeError, 
                "invalid : opcode signature doesn't match buffer (%c vs %c) at %i", sig, fullsig[arg], argloc);
                return -1;
            }
        }
    }
    return 0;
}



static int
NumExpr_init(NumExprObject *self, PyObject *args, PyObject *kwds)
{
    int i, j, mem_offset;
    int n_constants, n_inputs, n_temps;
    PyObject *signature = NULL, *tempsig = NULL, *constsig = NULL;
    PyObject *fullsig = NULL, *program = NULL, *constants = NULL;
    PyObject *input_names = NULL, *o_constants = NULL;
    char **mem = NULL, *rawmem = NULL;
    int *memsteps, rawmemsize;
    static char *kwlist[] = {"signature", "tempsig",
                             "program",  "constants",
                             "input_names", NULL};

    if (!PyArg_ParseTupleAndKeywords(args, kwds, "SSS|OO", kwlist,
                                     &signature,
                                     &tempsig,
                                     &program, &o_constants, 
                                     &input_names)) {
        return -1;
    }
    
    n_inputs = PyString_Size(signature);
    n_temps = PyString_Size(tempsig);
    
    if (o_constants) {
        if (!PySequence_Check(o_constants) ) {
                PyErr_SetString(PyExc_TypeError, "constants must be a sequence");
                return -1;
        }
        n_constants = PySequence_Length(o_constants);
        if (!(constants = PyTuple_New(n_constants)))
            return -1;
        if (!(constsig = PyString_FromStringAndSize(NULL, n_constants))) {
            Py_DECREF(constants);
            return -1;
        }
        for (i = 0; i < n_constants; i++) {
            PyObject *o;
            if (!(o = PySequence_GetItem(o_constants, i))) {
                Py_DECREF(constants);
                Py_DECREF(constsig);
                return -1;
            }
            Py_INCREF(o);
            PyTuple_SET_ITEM(constants, i, o);
            if (PyInt_Check(o)) {
                PyString_AS_STRING(constsig)[i] = 'i';
                continue;
            }
            if (PyFloat_Check(o)) {
                PyString_AS_STRING(constsig)[i] = 'f';
                continue;
            }
            if (PyComplex_Check(o)) {
                PyString_AS_STRING(constsig)[i] = 'c';
                continue;
            }
            PyErr_SetString(PyExc_TypeError, "constants must be of type int/float/complex");
            Py_DECREF(constsig);
            Py_DECREF(constants);
            return -1;
        }
    } else {
        n_constants = 0;
        if (!(constants = PyTuple_New(0)))
            return -1;
        if (!(constsig = PyString_FromString(""))) {
            Py_DECREF(constants);
            return -1;
        }
    }

    fullsig = PyString_FromFormat("%c%s%s%s", get_return_sig(program),
        PyString_AS_STRING(signature), PyString_AS_STRING(constsig), 
        PyString_AS_STRING(tempsig));
    if (!fullsig) {
        Py_DECREF(constants);
        Py_DECREF(constsig);
    }

    if (!input_names) {
        input_names = Py_None;
    }
    
    rawmemsize = BLOCK_SIZE1 * (size_from_sig(constsig) + size_from_sig(tempsig));
    mem = PyMem_New(void *, 1 + n_inputs + n_constants + n_temps);
    rawmem = PyMem_New(char, rawmemsize);
    memsteps = PyMem_New(int, 1 + n_inputs);
    if (!mem || !rawmem || !memsteps) {
        Py_DECREF(constants);
        Py_DECREF(constsig);
        Py_DECREF(fullsig);
        PyMem_Del(mem);    
        PyMem_Del(rawmem);
        PyMem_Del(memsteps);
        return -1;
    }
    /* 
       0                                                  -> output 
       [1, n_inputs+1)                                    -> inputs 
       [n_inputs+1, n_inputs+n_consts+1)                  -> constants
       [n_inputs+n_consts+1, n_inputs+n_consts+n_temps+1) -> temps
    */
    /* Fill in 'memsteps' for  output and inputs. 'mem' filled on call. */
    memsteps[0] = size_from_char(get_return_sig(program));
    for (i = 0; i < n_inputs; i++) {
        memsteps[i+1] = size_from_char(PyString_AS_STRING(signature)[i]);
    }
    /* Fill in 'mem' and 'rawmem' for constants */
    mem_offset = 0;
    for (i = 0; i < n_constants; i++) {
        char c = PyString_AS_STRING(constsig)[i];
        mem[i+n_inputs+1] = rawmem + mem_offset;
        mem_offset += BLOCK_SIZE1 * size_from_char(c);
        /* fill in the constants */
        if (c == 'i') {
            long *imem = (long*)mem[i+n_inputs+1];
            long value = PyInt_AS_LONG(PyTuple_GET_ITEM(constants, i));
            for (j = 0; j < BLOCK_SIZE1; j++) 
                imem[j] = value;
        } else if (c == 'f') {
            double *dmem = (double*)mem[i+n_inputs+1];
            double value = PyFloat_AS_DOUBLE(PyTuple_GET_ITEM(constants, i));
            for (j = 0; j < BLOCK_SIZE1; j++) 
                dmem[j] = value;
        } else if (c == 'c') {
            double *cmem = (double*)mem[i+n_inputs+1];
            Py_complex value = PyComplex_AsCComplex(PyTuple_GET_ITEM(constants, i));
            for (j = 0; j < 2*BLOCK_SIZE1; j+=2) { 
                cmem[j] = value.real;
                cmem[j+1] = value.imag;
            }
        }          
    }    
    /* Fill in 'mem' for temps */
    for (i = 0; i < n_temps; i++) {
        mem[i+n_inputs+n_constants+1] = rawmem + mem_offset;
        mem_offset += BLOCK_SIZE1 * size_from_char(PyString_AS_STRING(tempsig)[i]);
    }     
    /* See if any errors occured (e.g., in size_from_char) or if mem_offset is wrong */
    if (PyErr_Occurred() || mem_offset != rawmemsize) {
        if (mem_offset != rawmemsize) {
            PyErr_Format(PyExc_RuntimeError, "mem_offset does not match rawmemsize");
        }
        Py_DECREF(constants);
        Py_DECREF(constsig);
        Py_DECREF(fullsig);
        PyMem_Del(mem);    
        PyMem_Del(rawmem);
        PyMem_Del(memsteps);
        return -1;
    }

    
    #define REPLACE_OBJ(arg) \
    {PyObject *tmp = self->arg; \
     self->arg = arg; \
     Py_XDECREF(tmp);} 
    #define INCREF_REPLACE_OBJ(arg) {Py_INCREF(arg); REPLACE_OBJ(arg);}
    #define REPLACE_MEM(arg) {PyMem_Del(self->arg); self->arg=arg;}
    
    INCREF_REPLACE_OBJ(signature);
    INCREF_REPLACE_OBJ(tempsig);
    REPLACE_OBJ(constsig);
    REPLACE_OBJ(fullsig);
    INCREF_REPLACE_OBJ(program);
    REPLACE_OBJ(constants);
    INCREF_REPLACE_OBJ(input_names);
    REPLACE_MEM(mem);
    REPLACE_MEM(rawmem);
    REPLACE_MEM(memsteps);
    self->rawmemsize = rawmemsize;
    
    #undef REPLACE_OBJ
    #undef INCREF_REPLACE_OBJ
    #undef REPLACE_MEM

    #define REPLACE_MEM
    
    return check_program(self);
}

static PyMemberDef NumExpr_members[] = {  
    {"signature", T_OBJECT_EX, offsetof(NumExprObject, signature), READONLY, NULL},
    {"constsig", T_OBJECT_EX, offsetof(NumExprObject, constsig), READONLY, NULL},
    {"tempsig", T_OBJECT_EX, offsetof(NumExprObject, tempsig), READONLY, NULL},
    {"fullsig", T_OBJECT_EX, offsetof(NumExprObject, fullsig), READONLY, NULL},

    {"program", T_OBJECT_EX, offsetof(NumExprObject, program), READONLY, NULL},
    {"constants", T_OBJECT_EX, offsetof(NumExprObject, constants),
     READONLY, NULL},
    {"input_names", T_OBJECT, offsetof(NumExprObject, input_names), 0, NULL},
    {NULL},
};

struct vm_params {
    int prog_len;
    unsigned char *program;
    unsigned int n_inputs;
    unsigned int r_end;
    char *output;
    char **inputs;
    char **mem;
    int  *memsteps;
};

#define DO_BOUNDS_CHECK 1

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
#define VECTOR_SIZE BLOCK_SIZE1
#include "interp_body.c"
#undef VECTOR_SIZE
    }
    return 0;
}

static inline int
vm_engine_2(int start, int blen, struct vm_params params, int *pc_error)
{
    unsigned int index;
    for (index = start; index < blen; index += BLOCK_SIZE2) {
#define VECTOR_SIZE BLOCK_SIZE2
#include "interp_body.c"
#undef VECTOR_SIZE
    }
    return 0;
}

static inline int
vm_engine_rest(int start, int blen, struct vm_params params, int *pc_error)
{
    unsigned int index = start;
    unsigned int rest = blen - start;
#define VECTOR_SIZE rest
#include "interp_body.c"
#undef VECTOR_SIZE
    return 0;
}

static int
run_interpreter(NumExprObject *self, int len, char *output, char **inputs,
                int *pc_error)
{
    int r;
    unsigned int t, blen1, blen2, n_constants;
    struct vm_params params;

    *pc_error = -1;
    if (PyString_AsStringAndSize(self->program, (char **)&(params.program),
                                 &(params.prog_len)) < 0) {
        return -1;
    }
    if ((params.n_inputs = PyObject_Length(self->signature)) == -1)
        return -1;
    params.output = output;
    params.inputs = inputs;
    params.mem = self->mem;
    params.memsteps = self->memsteps;
    params.r_end = PyString_Size(self->fullsig);
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
    PyObject *output = NULL, *a_inputs = NULL;
    unsigned int n_inputs;
    int i, len = -1, r, pc_error;
    char **inputs;

    n_inputs = PyTuple_Size(args);
    if (PyString_Size(self->signature) != n_inputs) {
        return PyErr_Format(PyExc_ValueError,
                            "number of inputs doesn't match program");
    }
    if (kwds && PyObject_Length(kwds) > 0) {
        return PyErr_Format(PyExc_ValueError,
                            "keyword arguments are not accepted");
    }
    a_inputs = PyTuple_New(n_inputs);
    if (!a_inputs) {
        return NULL;
    }
    inputs = PyMem_New(char *, n_inputs);
    if (!inputs) {
        Py_DECREF(a_inputs);
        return NULL;
    }

    if (n_inputs == 0) {
        /* allocate one space for scalar result */
        char retsig = get_return_sig(self->program);
        intp dims[1];
        len = 1;
        output = PyArray_SimpleNew(0,
                                   &dims,
                                   typecode_from_char(retsig));    
        if (!output) {
            Py_XDECREF(a_inputs);
            PyMem_Del(inputs);
        }
    }
    for (i = 0; i < n_inputs; i++) {
        PyObject *o = PyTuple_GetItem(args, i); /* borrowed ref */
        PyObject *a;
        intp typecode = typecode_from_char(PyString_AS_STRING(self->signature)[i]);
        if (typecode == -1) return NULL;
        a = PyArray_ContiguousFromAny(o, typecode, 0, 0);
        if (!a) return NULL;
        PyTuple_SET_ITEM(a_inputs, i, a);  /* steals reference */
        if (len == -1) { 
            char retsig = get_return_sig(self->program);
            if (retsig != 'i' && retsig != 'f' && retsig != 'c') {
                Py_XDECREF(a_inputs);
                PyMem_Del(inputs);
                return PyErr_Format(PyExc_RuntimeError,
                                    "bad return format %c", retsig);
            } 
            len = PyArray_SIZE(a);
            output = PyArray_SimpleNew(PyArray_NDIM(a),
                                       PyArray_DIMS(a),
                                       typecode_from_char(retsig));    
            if (!output) {
                Py_XDECREF(a_inputs);
                PyMem_Del(inputs);
            }
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
    
    add_op("ones_like_ii", OP_ONES_LIKE_II);
    add_op("copy_ii", OP_COPY_II);
    add_op("neg_ii", OP_NEG_II);
    add_op("add_iii", OP_ADD_III);
    add_op("sub_iii", OP_SUB_III);
    add_op("mul_iii", OP_MUL_III);
    add_op("div_iii", OP_DIV_III);
    add_op("pow_iii", OP_POW_III);
    add_op("mod_iii", OP_MOD_III);
    add_op("where_ifii", OP_WHERE_IFII);
    
    add_op("gt_iff", OP_GT_IFF);
    add_op("ge_iff", OP_GE_IFF);
    add_op("eq_iff", OP_EQ_IFF);
    add_op("ne_iff", OP_NE_IFF);
    
    add_op("cast_fi", OP_CAST_FI);
    add_op("copy_ff", OP_COPY_FF);
    add_op("ones_like_ff", OP_ONES_LIKE_FF);
    add_op("neg_cc", OP_NEG_CC);
    add_op("neg_ff", OP_NEG_FF);
    add_op("add_fff", OP_ADD_FFF);
    add_op("sub_fff", OP_SUB_FFF);
    add_op("mul_fff", OP_MUL_FFF);
    add_op("div_fff", OP_DIV_FFF);
    add_op("pow_fff", OP_POW_FFF);
    add_op("mod_fff", OP_MOD_FFF);
    add_op("sin_ff", OP_SIN_FF);
    add_op("cos_ff", OP_COS_FF);
    add_op("tan_ff", OP_TAN_FF);
    add_op("sqrt_ff", OP_SQRT_FF);
    add_op("arctan2_fff", OP_ARCTAN2_FFF);
    add_op("where_ffff", OP_WHERE_FFFF);
    add_op("func_ff", OP_FUNC_FF);
    add_op("func_fff", OP_FUNC_FFF);
    
    add_op("eq_icc", OP_EQ_ICC);
    add_op("ne_icc", OP_NE_ICC);
    
    add_op("cast_ci", OP_CAST_CI);
    add_op("cast_cf", OP_CAST_CF);
    add_op("copy_cc", OP_COPY_CC);
    add_op("ones_like_cc", OP_ONES_LIKE_CC);
    add_op("neg_cc", OP_NEG_CC);
    add_op("add_ccc", OP_ADD_CCC);
    add_op("sub_ccc", OP_SUB_CCC);
    add_op("mul_ccc", OP_MUL_CCC);
    add_op("div_ccc", OP_DIV_CCC);
    add_op("where_cfcc", OP_WHERE_CFCC);
    add_op("func_cc", OP_FUNC_CC);
    add_op("func_ccc", OP_FUNC_CCC);
    
    add_op("real_fc", OP_REAL_FC);
    add_op("imag_fc", OP_IMAG_FC);
    add_op("complex_cff", OP_COMPLEX_CFF);
    
#undef add_op

    if (PyModule_AddObject(m, "opcodes", d) < 0) return;

    d = PyDict_New();
    if (!d) return;
    
#define add_func(sname, name) o = PyInt_FromLong(name); \
    r = PyDict_SetItemString(d, sname, o);              \
    Py_XDECREF(o);                                      \
    if (r < 0) {PyErr_SetString(PyExc_RuntimeError, "add_func"); return;}

    add_func("sqrt_ff", FUNC_SQRT_FF);
    add_func("sin_ff", FUNC_SIN_FF);
    add_func("cos_ff", FUNC_COS_FF);
    add_func("tan_ff", FUNC_TAN_FF);
    add_func("arcsin_ff", FUNC_ARCSIN_FF);
    add_func("arccos_ff", FUNC_ARCCOS_FF);
    add_func("arctan_ff", FUNC_ARCTAN_FF);
    add_func("sinh_ff", FUNC_SINH_FF);
    add_func("cosh_ff", FUNC_COSH_FF);
    add_func("tanh_ff", FUNC_TANH_FF);
    
    add_func("fmod_fff", FUNC_FMOD_FFF);
    
    add_func("sqrt_cc", FUNC_SQRT_CC);
    add_func("sin_cc", FUNC_SIN_CC);
    add_func("cos_cc", FUNC_COS_CC);
    add_func("tan_cc", FUNC_TAN_CC);
    add_func("arcsin_cc", FUNC_ARCSIN_CC);
    add_func("arccos_cc", FUNC_ARCCOS_CC);
    add_func("arctan_cc", FUNC_ARCTAN_CC);
    add_func("sinh_cc", FUNC_SINH_CC);
    add_func("cosh_cc", FUNC_COSH_CC);
    add_func("tanh_cc", FUNC_TANH_CC);
    
    add_func("pow_ccc", FUNC_POW_CCC);

#undef add_func

   if (PyModule_AddObject(m, "funccodes", d) < 0) return;

}
