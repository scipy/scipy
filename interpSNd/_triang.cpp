/*This will become an extension module used to speed up the slow steps
    in InterpolateSNd.
    
    The best step to speed up is calculation of the circumcircle.
*/

#include "Python.h"
#include <stdlib.h>
#include <map>
#include <iostream>

#include "numpy/arrayobject.h" // required to use arrays
#include "numpy/ufuncobject.h"

using namespace std;

extern "C" {
    
/* test function */
PyObject * sayHello(PyObject *self, PyObject *args)
{
    int x;
    PyObject *other_arg;
    if (!PyArg_ParseTuple(args, "iO", x, &other_arg)){
        return NULL;
    }
    printf("\nhello world\n");
    return other_arg;
}



/* merely return first element of a matrix */
PyObject * upperleft(PyObject *self, PyObject *args)
{
    PyArrayObject *data_coordinates;
    
    // parse tuple
    if (!PyArg_ParseTuple(args, "O!", &PyArray_Type, &data_coordinates)){
        return NULL;
    }
    
    // make sure floats
    if (data_coordinates->descr->type_num != PyArray_DOUBLE){
        PyErr_SetString(PyExc_ValueError, 
            "array must be of type float");
        return NULL;
    }
    
    //first element of array
    double ans;
    ans = *(double *)(data_coordinates->data);
    
    return PyFloat_FromDouble(ans);
    
}


/* tell Python which methods are available */
static PyMethodDef triang_methods[]={
    {"sayHello", (PyCFunction)sayHello, METH_VARARGS, "says hi"},
    {"upperleft", (PyCFunction)upperleft, METH_VARARGS, ""},
    //{"lookup", (PyCFunction)lookup, METH_VARARGS, ""},*/
    {NULL, NULL, 0, NULL}
};


/* initialize the module */
PyMODINIT_FUNC init_triang(void)
{
    PyObject* m;
    m = Py_InitModule3("_triang", triang_methods,
        "stupid, useless module that will get better\n"
    );
    if (m == NULL)
        return;
    import_array(); //required at initialization in order to use arrays
}


    
    
    
    
} //extern "C"