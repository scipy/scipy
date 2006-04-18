/* Python module for fast variate generation from a discrete distribution.  Uses
 * Marsaglia's compact 5-table method, as described in his paper 'Fast
 * generation of discrete random variables' in the Journal of Statistical
 * Software, July 2004, vol 11, issue 3.  
 *
 * The underlying algorithms are in the file sampler5tbl.c.  This code is based
 * upon the C implementation that accompanies that paper, but is simpler, and
 * uses a different random number generator, the Mersenne Twister in
 * Jean-Sebastien Roy's RandomKit.
 *
 * Copyright: Ed Schofield, 2005-6
 * License: BSD-style (see LICENSE.txt at root of scipy tree)
 */


#include "Python.h"
#include "numpy/arrayobject.h"
#include "compact5table.h"

// Function prototypes:
static PyArrayObject* PyArray_Intsample(Sampler* mysampler, unsigned long size);
static PyObject* sample(PyObject *self, PyObject *args, PyObject *keywords);


// IntSampler type
typedef struct {
    PyObject_HEAD
    Sampler* pSampler;
} IntSampler;

//    "destroy" should be called automatically upon deletion
static void IntSampler_destroy(IntSampler* self)
{
    // printf("[Started destroying sampler]\n");
    if (self->pSampler != NULL)
        destroy_sampler5tbl(self->pSampler);
    self->ob_type->tp_free((PyObject*)self);
    // printf("[Finished destroying sampler]\n");
}

static PyObject *
IntSampler_new(PyTypeObject *type, PyObject *args, PyObject *kwds)
{
    IntSampler *self;

    // printf("[Started new() method]\n");
    self = (IntSampler *)type->tp_alloc(type, 0);
    if (self != NULL)
        self->pSampler = NULL;
    // printf("[Finished new() method]\n");
    return (PyObject *) self;
}

static int
IntSampler_init(IntSampler *self, PyObject *args, PyObject *kwds)
{
    PyArrayObject *pmf_table;
    int k;                              /* size of the table */
    unsigned long seed = 0;             /* Initialize with random seed. This can
                                         * be set explicitly by calling seed().
                                         */

    // printf("[Started initializing sampler]\n");

    /* parse the arguments  */
    if (!PyArg_ParseTuple(args, "O!", &PyArray_Type, &pmf_table)) {
      return -1;  /* Error indicator */
    }

    // printf("[Parsed arguments]\n");

    /* check that we have a one-dimensional array */
    if (pmf_table->nd != 1) {
      /* throw a Python exception (ValueError) */
      PyErr_SetString(PyExc_ValueError,
                      "the pmf table must be a one-dimensional array");
      return -1;
    }

    // printf("[Array is 1D]\n");

    /* check that the data type is float64, (C double) */
    if (PyArray_DESCR(pmf_table)->type_num != PyArray_DOUBLE) {
      PyErr_SetString(PyExc_ValueError,
                      "the pmf table must be of type float64");
      return -1;
    }

    // printf("[Data type is float64]\n");

    k = pmf_table->dimensions[0];    /* length of the array */

    self->pSampler = init_sampler5tbl((double*) pmf_table->data, k, seed);
    if (self->pSampler == NULL)
    {
        Py_DECREF(self);
        return -1;
    }
    
    // printf("[Finished initializing sampler]\n");

    // Eventually, we need to do this:
    // destroy(self->pSampler);

    return 0;
}



static char sample__doc__[] = \
  "sample(size): return an array with a random discrete sample\n"\
  "of the given size from the probability mass function specified when\n"\
  "initializing the sampler.\n";
                                                      
static PyObject*
sample(PyObject *self, PyObject *args, PyObject *keywords)
{
    PyArrayObject *samplearray;

    static char *kwlist[] = {"size",NULL};
    int size;
    
    /* parse the arguments  */
    if (!PyArg_ParseTupleAndKeywords(args, keywords, "i", kwlist, &size))
        return NULL;

    // printf("[Parsed arguments]\n");

    // Check that size > 0
    if (size <= 0)
    {
        PyErr_SetString(PyExc_ValueError, "sample size must be positive");
        return NULL;
    }

    // printf("[size is > 0]\n");

    samplearray = PyArray_Intsample(((IntSampler*)self)->pSampler, size);

    return (PyObject*) samplearray;
}



static PyArrayObject*
PyArray_Intsample(Sampler* mysampler, unsigned long size)
{
    PyArrayObject* samplearray;
    unsigned long* ptr;
    
    int ndim = 1;
    int dims[1] = {size};
    int typenum = PyArray_LONG;
    samplearray = (PyArrayObject*) PyArray_SimpleNew(ndim, dims,
                                                     typenum);
    if (samplearray == NULL) return NULL;

    ptr = (unsigned long*) PyArray_DATA(samplearray);
    Dran_array(mysampler, ptr, size);

    return samplearray;
}

 

/* Doc strings: */
static char intsampler__doc__[] = \
  "A module allowing fast sampling from a discrete distribution given its\n"\
  "probability mass function.\n"\
  "\n"\
  "Use the syntax:\n"\
  ">>> s = _intsampler(table)\n"\
  "to create an object for sampling from a distribution with probability\n"\
  "mass function given by 'table', where:\n"\
  "\n"\
  "x       0           1           ...     k-1\n"\
  "p(x)    table[0]    table[1]            table[k-1]\n"\
  "\n"\
  "The values of table[i] need not be normalized to sum to 1, but must be\n"\
  "non-negative.\n";



/* Module functions */
static PyMethodDef module_methods[] = {
    {NULL}  /* Sentinel */
};

/* IntSampler class methods */
static PyMethodDef IntSampler_methods[] = {
//        {"setup",    setup, METH_VARARGS, setup__doc__},
        {"sample", 
         (PyCFunction)sample,
         METH_VARARGS | METH_KEYWORDS,
         sample__doc__},
        {NULL,          NULL}           /* sentinel */
};


static PyTypeObject IntSamplerType = {
    PyObject_HEAD_INIT(NULL)
    0,                         /*ob_size*/
    "_intsampler._intsampler",             /*tp_name*/
    sizeof(IntSampler), /*tp_basicsize*/
    0,                         /*tp_itemsize*/
    (destructor)IntSampler_destroy,  /*tp_dealloc*/
    0,                         /*tp_print*/
    0,                         /*tp_getattr*/
    0,                         /*tp_setattr*/
    0,                         /*tp_compare*/
    0,                         /*tp_repr*/
    0,                         /*tp_as_number*/
    0,                         /*tp_as_sequence*/
    0,                         /*tp_as_mapping*/
    0,                         /*tp_hash */
    0,                         /*tp_call*/
    0,                         /*tp_str*/
    0,                         /*tp_getattro*/
    0,                         /*tp_setattro*/
    0,                         /*tp_as_buffer*/
    Py_TPFLAGS_DEFAULT | Py_TPFLAGS_BASETYPE,        /*tp_flags*/
    "IntSampler objects",      /*tp_doc*/
    0,                         /* tp_traverse */
    0,                         /* tp_clear */
    0,                         /* tp_richcompare */
    0,                         /* tp_weaklistoffset */
    0,                         /* tp_iter */
    0,                         /* tp_iternext */
    IntSampler_methods,        /* tp_methods */
    0,                         /* tp_members */
    0,                         /* tp_getset */
    0,                         /* tp_base */
    0,                         /* tp_dict */
    0,                         /* tp_descr_get */
    0,                         /* tp_descr_set */
    0,                         /* tp_dictoffset */
    (initproc)IntSampler_init,      /* tp_init */
    0,                         /* tp_alloc */
    IntSampler_new,            /* tp_new */
};



#ifndef PyMODINIT_FUNC  /* declarations for shared library import/export */
#define PyMODINIT_FUNC void
#endif
PyMODINIT_FUNC
init_intsampler(void)
{
    PyObject *m, *d, *s;
    
    /* Initialize scipy */
    import_array();   
    /* Import the ufunc objects */
    // import_ufunc();

    // IntSamplerType.tp_new = PyType_GenericNew;
    if (PyType_Ready(&IntSamplerType) < 0)
        return;

    /* Create the module and add the functions */
    m = Py_InitModule3("_intsampler", module_methods, intsampler__doc__); 
    
    if (m == NULL)
        return;
    
    // /* Add our class */
    // PyStructSequence_InitType(&SamplerType, &sampler_type_desc);
    Py_INCREF(&IntSamplerType);
    PyModule_AddObject(m, "_intsampler", (PyObject*) &IntSamplerType);


    /* Add some symbolic constants to the module */
    d = PyModule_GetDict(m);

    s = PyString_FromString("2.0-alpha5");
    PyDict_SetItemString(d, "__version__", s);
    Py_DECREF(s);

  /* Check for errors */
  if (PyErr_Occurred())
    Py_FatalError("can't initialize module _intsampler");
}






