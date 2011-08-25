#include "Python.h"
#include "math.h"
#include "numpy/ndarraytypes.h"
#include "numpy/ufuncobject.h"


/* 
    ufuncs to compute logit(p) = log(p/(1-p)) and
    expit(x) = 1/(1+exp(-x)).
*/


static PyMethodDef LogitMethods[] = {
        {NULL, NULL, 0, NULL}
};


/* 
    The three different inner loops for logit, corresponding to 
    long_double, double, and float dtypes.
*/

static void long_double_logit(char **args, npy_intp *dimensions,
                            npy_intp* steps, void* data){

    npy_intp i;
    npy_intp n=dimensions[0];
    char *in=args[0], *out=args[1];
    npy_intp in_step=steps[0], out_step=steps[1];
    
    long double tmp;
    
    for(i=0; i<n; i++){
        tmp = *(long double *)in;
        tmp /= 1-tmp;
        *((long double *)out) = logl(tmp);
        
        in += in_step;
        out += out_step;
    }
}

static void double_logit(char **args, npy_intp *dimensions,
                            npy_intp* steps, void* data){

    npy_intp i;
    npy_intp n=dimensions[0];
    char *in=args[0], *out=args[1];
    npy_intp in_step=steps[0], out_step=steps[1];
    
    double tmp;
    
    for(i=0; i<n; i++){
        tmp = *(double *)in;
        tmp /= 1-tmp;
        *((double *)out) = log(tmp);
        
        in += in_step;
        out += out_step;
    }
}

static void float_logit(char **args, npy_intp *dimensions,
                            npy_intp* steps, void* data){

    npy_intp i;
    npy_intp n=dimensions[0];
    char *in=args[0], *out=args[1];
    npy_intp in_step=steps[0], out_step=steps[1];
    
    float tmp;
    
    for(i=0; i<n; i++){
        tmp = *(float *)in;
        tmp /= 1-tmp;
        *((float *)out) = logf(tmp);
        
        in += in_step;
        out += out_step;
    }
}

/*
    The three different inner loops for expit, corresponding 
    to long_double, double, and float dtypes.
*/


static void long_double_expit(char **args, npy_intp *dimensions,
                            npy_intp* steps, void* data){

    npy_intp i;
    npy_intp n=dimensions[0];
    char *in=args[0], *out=args[1];
    npy_intp in_step=steps[0], out_step=steps[1];
    
    long double tmp;
    
    for(i=0; i<n; i++){
        tmp = *(long double *)in;
         if( tmp>0 ){
            tmp = expl(tmp);
            *((long double *)out) = tmp/(1+tmp);
        }
        else{
            *((long double *)out) = 1/(1+expl(-tmp));
        }
        in += in_step;
        out += out_step;
    }
}

static void double_expit(char **args, npy_intp *dimensions,
                            npy_intp* steps, void* data){

    npy_intp i;
    npy_intp n=dimensions[0];
    char *in=args[0], *out=args[1];
    npy_intp in_step=steps[0], out_step=steps[1];
    
    double tmp;
    
    for(i=0; i<n; i++){
        tmp = *(double *)in;
        if( tmp>0 ){
            tmp = expl(tmp);
            *((double *)out) = tmp/(1+tmp);
        }
        else{
            *((double *)out) = 1/(1+exp(-tmp));
        }
        
        in += in_step;
        out += out_step;
    }
}

static void float_expit(char **args, npy_intp *dimensions,
                            npy_intp* steps, void* data){

    npy_intp i;
    npy_intp n=dimensions[0];
    char *in=args[0], *out=args[1];
    npy_intp in_step=steps[0], out_step=steps[1];
    
    float tmp;
    
    for(i=0; i<n; i++){
        tmp = *(float *)in;
           if( tmp>0 ){
            tmp = expl(tmp);
            *((float *)out) = tmp/(1+tmp);
        }
        else{
            *((float *)out) = 1/(1+expf(-tmp));
        }
        
        in += in_step;
        out += out_step;
    }
}

/*Definitions for generating the ufuncs. */

PyUFuncGenericFunction expit_funcs[3] = {&float_expit, 
                                    &double_expit, &long_double_expit};

PyUFuncGenericFunction logit_funcs[3] = {&float_logit, 
                                    &double_logit, &long_double_logit};
char types[6] = {NPY_FLOAT, NPY_FLOAT, 
                NPY_DOUBLE,NPY_DOUBLE, 
                NPY_LONGDOUBLE, NPY_LONGDOUBLE};
void *data[3] = {NULL, NULL, NULL};

/* Module init */

#if PY_VERSION_HEX >= 0x03000000

static struct PyModuleDef moduledef = {
    PyModuleDef_HEAD_INIT,
    "logit",
    NULL,
    -1,
    LogitMethods,
    NULL,
    NULL,
    NULL,
    NULL
};

PyObject *PyInit_logit(void)
{
    PyObject *m, *f, *d;
    m = PyModule_Create(&moduledef);
    if(!m){
        return NULL;
    }

    import_array();
    import_umath();

    f = PyUFunc_FromFuncAndData(logit_funcs,data, types, 3, 1, 1, 
                                PyUFunc_None, "logit",NULL , 0);
    PyDict_SetItemString(d , "logit", f);
    Py_DECREF(f);


    f = PyUFunc_FromFuncAndData(expit_funcs,data, types, 3, 1, 1, 
                                PyUFunc_None, "expit",NULL , 0);
    PyDict_SetItemString(d , "expit", f);
    Py_DECREF(f);

    return m;
}

#else

PyMODINIT_FUNC initlogit(void){
    PyObject *m, *f,  *d;
    
    
    m  =Py_InitModule("logit", LogitMethods);
    if( m==NULL ){
        return;
    }

    d = PyModule_GetDict(m);
    
    import_array();
    import_umath();
   
    f = PyUFunc_FromFuncAndData(logit_funcs,data, types, 3, 1, 1, 
                                PyUFunc_None, "logit",NULL , 0);
    PyDict_SetItemString(d , "logit", f);
    Py_DECREF(f);


    f = PyUFunc_FromFuncAndData(expit_funcs,data, types, 3, 1, 1, 
                                PyUFunc_None, "expit",NULL , 0);
    PyDict_SetItemString(d , "expit", f);
    Py_DECREF(f);
}
#endif
