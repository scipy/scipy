
/* Written by Charles Harris charles.harris@sdl.usu.edu */

/* Modifications by Travis Oliphant to separate Python code from C
   routines */

#include <setjmp.h>
#include "Python.h"

typedef struct {
    int funcalls;
    int iterations;
    int error_num; 
    PyObject *function;
    PyObject *args;
    jmp_buf env;
} scipy_zeros_parameters;

/*
 * Storage for the relative precision of doubles. This is computed when the module
 * is initialized.
 */

#include "Zeros/zeros.h"

/*
extern double brentq();
extern double brenth();
extern double ridder();
extern double bisect();
*/
#define SIGNERR -1
#define CONVERR -2

static double scipy_zeros_rtol=0;

double 
scipy_zeros_functions_func(double x, void *params)
{
    scipy_zeros_parameters *myparams = params;
    PyObject *args, *f, *retval=NULL;
    double val;

    args = myparams->args;
    f = myparams->function;
    PyTuple_SetItem(args,0,Py_BuildValue("d",x));
    val = PyFloat_AsDouble(retval=PyObject_CallObject(f,args));
    Py_XDECREF(retval);
    if (PyErr_Occurred()) {
        fprintf(stderr, "Internal Error occured.\n");
        longjmp(myparams->env, 1);
    }
    return val;    
}

/*
 * Helper function that calls a Python function with extended arguments
 */

static PyObject *
call_solver(double (*solver)(), PyObject *self, PyObject *args)
{    
    double a,b,xtol,zero;
    int iter,i, len, fulloutput, disp=1, flag=0;
    scipy_zeros_parameters params;
    jmp_buf env;
    PyObject *f,*xargs,*item,*fargs=NULL;

    if (!PyArg_ParseTuple(args,"OdddiOi|i",&f,&a,&b,&xtol,&iter,&xargs,&fulloutput,&disp)) 
        {
            PyErr_SetString(PyExc_RuntimeError,"Unable to parse arguments");
            return NULL;
        }
    if (xtol < 0) {
        PyErr_SetString(PyExc_ValueError,"xtol must be >= 0");
        return NULL;
    }
    if (iter < 0) {
        PyErr_SetString(PyExc_ValueError,"maxiter should be > 0");
        return NULL;
    }
    
    len = PyTuple_Size(xargs);
    fargs = PyTuple_New(len + 1);  /* room for the double
                                      as the first argument */
    if (fargs == NULL) {
        PyErr_SetString(PyExc_RuntimeError,"Failed to allocate argument tuple");
        return NULL;
    }

    for (i = 0; i < len; i++) {
        item = PyTuple_GetItem(xargs, i);
        if (item == NULL) { Py_DECREF(fargs); return NULL;}
        Py_INCREF(item);
        PyTuple_SET_ITEM(fargs,i+1,item);
    }

    params.function = f;
    params.args = fargs;

    if (!setjmp(env)) {  /* direct return */
        memcpy(params.env,env,sizeof(jmp_buf));
        params.error_num = 0;
        zero = solver(scipy_zeros_functions_func,a,b,xtol,scipy_zeros_rtol,iter, &params);    
        Py_DECREF(fargs);
        if (params.error_num != 0) {
            if (params.error_num == SIGNERR) {
                PyErr_SetString(PyExc_ValueError,"f(a) and f(b) must have different signs");
                return NULL;
            }
            if (params.error_num == CONVERR) {
                if (disp) {
                    fprintf(stderr, "Warning: failed to converge after %d iterations.\n", params.iterations);
                    flag = 1;
                }
            }
        }
        if (fulloutput) return Py_BuildValue("diii",zero,params.funcalls,params.iterations,flag);
        else return Py_BuildValue("d",zero);
    }
    else {  /* error return from Python function */
        fprintf(stderr, "Error when calling Python function.  See traceback.\n");
        Py_DECREF(fargs);
        return NULL;        
    }
}

/*
 * These routines interface with the solvers through call_solver
 */

static PyObject *
_bisect(PyObject *self, PyObject *args)
{
        return call_solver(bisect,self,args);
}

static PyObject *
_ridder(PyObject *self, PyObject *args)
{
        return call_solver(ridder,self,args);
}

static PyObject *
_brenth(PyObject *self, PyObject *args)
{
        return call_solver(brenth,self,args);
}

static PyObject *
_brentq(PyObject *self, PyObject *args)
{
        return call_solver(brentq,self,args);
}

/*
 * Standard Python module inteface
 */

static PyMethodDef 
Zerosmethods[] = {
	{"_bisect", _bisect, METH_VARARGS, "a"},
	{"_ridder", _ridder, METH_VARARGS, "a"},
	{"_brenth", _brenth, METH_VARARGS, "a"},
	{"_brentq", _brentq, METH_VARARGS, "a"},
	{NULL, NULL}
};

void
init_zeros(void)
{
        double tol;

        /* Determine relative precision of doubles, assumes binary */
        for(tol = 1; tol + 1 != 1; tol /= 2);
        scipy_zeros_rtol = 2*tol;

        Py_InitModule("_zeros", Zerosmethods);
}


