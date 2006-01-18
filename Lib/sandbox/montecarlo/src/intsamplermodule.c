#include "Python.h"
#include "numpy/arrayobject.h"
#include "sampler5tbl.h"

// Function prototypes:
static PyArrayObject* PyArray_Intsample(Sampler* mysampler, int size);
static PyObject* sample(PyObject *self, PyObject *args, PyObject *keywords);


// Global variable:
// Sampler* sampler_global;

// IntSampler type
typedef struct {
    PyObject_HEAD
    Sampler* pSampler;
} IntSampler;

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
    int k;  // size of the table

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

    /* check that the datatype is float64, (C double) */
    if (PyArray_DESCR(pmf_table)->type_num != PyArray_DOUBLE) {
      PyErr_SetString(PyExc_ValueError,
                      "the pmf table must be of type float64");
      return -1;
    }

    // printf("[Data type is float64]\n");

    k = pmf_table->dimensions[0];    /* length of the array */

    /* Init sampler */
    unsigned t = (unsigned) time( NULL );
    srand48(t);
    // printf("Seeded C RNG with time %u\n",t);

    // Check that a sample hasn't already been created but not destroyed
    // if (sampler_global != NULL)
    // {
    //     PyErr_SetString(PyExc_ValueError, "a sample exists that must be destroyed first.");
    //     return NULL;
    // }
    
    // if (sampler_global != NULL)
    //     destroy_sampler5tbl(sampler_global);


    self->pSampler = init_sampler5tbl((double*) pmf_table->data, k);
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
  // "\n"\
  // "If the optional argument 'probs' is zero (the default), this function\n"\
  // "returns only the generated sample.\n"
  // "\n"\
  // "If 'probs' is 1, it returns a tuple (sample, prob), where prob is an\n"\
  // "array of length len(sample) for which probs[i] is the (normalized) pmf\n"\
  // "value of sample[i] in the table used to initialize the sampler.\n"\
  // "\n"\
  // "If 'probs' is 2, it returns a tuple (sample, logprob), where logprob\n"\
  // "contains the natural logarithms of the probabilities of the sample points.";
                                                      
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
PyArray_Intsample(Sampler* mysampler, int size)
{
    PyArrayObject* samplearray;
    long* ptr;
    
    int ndim = 1;
    int dims[1] = {size};
    int typenum = PyArray_LONG;
    samplearray = (PyArrayObject*) PyArray_SimpleNew(ndim, dims,
                                                     typenum);
    if (samplearray == NULL) return NULL;

    ptr = (long*) PyArray_DATA(samplearray);
    Dran_array(mysampler, ptr, size);

    return samplearray;
}

 

/* Doc strings: */
// static char sumarray__doc__[] = "sumarray(a)";
static char intsampler__doc__[] = \
  "A module allowing fast sampling from a given discrete distribution.\n"\
  "\n"\
  "Use the syntax:\n"\
  ">>> s = intsampler(table)\n"\
  "to create an object for sampling from a distribution with probability\n"\
  "mass function given by 'table', where:\n"\
  "\n"\
  "x       0           1           ...     k-1\n"\
  "p(x)    table[0]    table[1]            table[k-1]\n"\
  "\n"\
  "The values of table[i] are normalized by dividing by sum(table).\n";



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
        //    "destroy" should be called automatically upon deletion
//     "Free the memory and destroy an existing sample of sentences."},
//        {"destroy",  destroy, METH_VARARGS, destroy__doc__},
        {NULL,          NULL}           /* sentinel */
};


static PyTypeObject IntSamplerType = {
    PyObject_HEAD_INIT(NULL)
    0,                         /*ob_size*/
    "intsampler.intsampler",             /*tp_name*/
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



#ifndef PyMODINIT_FUNC  /* declarations for DLL import/export */
#define PyMODINIT_FUNC void
#endif
PyMODINIT_FUNC
initintsampler(void)
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
    m = Py_InitModule3("intsampler", module_methods, intsampler__doc__); 
    
    if (m == NULL)
        return;
    
    // /* Add our class */
    // PyStructSequence_InitType(&SamplerType, &sampler_type_desc);
    Py_INCREF(&IntSamplerType);
    PyModule_AddObject(m, "intsampler", (PyObject*) &IntSamplerType);


    /* Add some symbolic constants to the module */
    d = PyModule_GetDict(m);

    s = PyString_FromString("2.0-alpha4");
    PyDict_SetItemString(d, "__version__", s);
    Py_DECREF(s);

  /* Check for errors */
  if (PyErr_Occurred())
    Py_FatalError("can't initialize module intsampler");
}






