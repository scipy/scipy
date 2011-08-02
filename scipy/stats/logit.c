#include "Python.h"
#include "math.h"
#include "numpy/ndarraytypes.h"
#include "numpy/ufuncobject.h"
#include "numpy/halffloat.h"


#define logit_doc "A logit ufunc.\nReturns logit(p) = log(p/(1-p))."

#define expit_doc "A expit ufunc.\nReturns returns expit(x) = 1/(1+exp(-x)). expit is the inverse of the logit function."


/* 
    ufuncs to compute logit(p) = log(p/(1-p)) and
    expit(x) = 1/(1+exp(-x)).
*/


static PyMethodDef LogitMethods[] = {
        {NULL, NULL, 0, NULL}
};


/* 
    The four different inner loops for logit, corresponding to 
    long_double, double, float, and half_float dtypes.
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


static void half_float_logit(char **args, npy_intp *dimensions,
                            npy_intp* steps, void* data){

    npy_intp i;
    npy_intp n=dimensions[0];
    char *in=args[0], *out=args[1];
    npy_intp in_step=steps[0], out_step=steps[1];
    
    float tmp;
    
    for(i=0; i<n; i++){
        tmp = *(npy_half *)in;

        tmp = npy_half_to_float(tmp);
        tmp /= 1-tmp;
        tmp = logf(tmp);

        *((npy_half *)out) = npy_float_to_half(tmp);
        
        in += in_step;
        out += out_step;
    }
}



/*
    The four different inner loops for expit, corresponding 
    to long_double, double, float, and half_float dtypes.
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


static void half_float_expit(char **args, npy_intp *dimensions,
                            npy_intp* steps, void* data){

    npy_intp i;
    npy_intp n=dimensions[0];
    char *in=args[0], *out=args[1];
    npy_intp in_step=steps[0], out_step=steps[1];
    
    float tmp;
    
    for(i=0; i<n; i++){
        tmp = *(npy_half *)in;
        tmp = npy_half_to_float(tmp);

        if( tmp>0 ){
            tmp = expl(tmp);
            tmp = tmp/(1+tmp);
        }
        else{
            tmp = 1/(1+expl(-tmp));
        }

        *((npy_half *)out) = npy_float_to_half(tmp);
        
        in += in_step;
        out += out_step;
    }
}

/*Definitions for generating the ufuncs. */

PyUFuncGenericFunction expit_funcs[4] = {&half_float_expit, &float_expit, 
                                    &double_expit, &long_double_expit};

PyUFuncGenericFunction logit_funcs[4] = {&half_float_logit, &float_logit, 
                                    &double_logit, &long_double_logit};
char types[8] = {NPY_HALF, NPY_HALF, NPY_FLOAT, NPY_FLOAT, 
                NPY_DOUBLE,NPY_DOUBLE, NPY_LONGDOUBLE, NPY_LONGDOUBLE};
void *data[4] = {NULL, NULL, NULL, NULL};

/* Module init */

PyMODINIT_FUNC initlogit(void){
    PyObject *m, *f,  *d;
    
    
    m  =Py_InitModule("logit", LogitMethods);
    if( m==NULL ){
        return;
    }

    d = PyModule_GetDict(m);
    
    import_array();
    import_umath();
   
    f = PyUFunc_FromFuncAndData(logit_funcs,data, types, 4, 1, 1, 
                                PyUFunc_None, "logit", logit_doc, 0);
    PyDict_SetItemString(d , "logit", f);
    Py_DECREF(f);


    f = PyUFunc_FromFuncAndData(expit_funcs,data, types, 4, 1, 1, 
                                PyUFunc_None, "expit", expit_doc, 0);
    PyDict_SetItemString(d , "expit", f);
    Py_DECREF(f);


}
