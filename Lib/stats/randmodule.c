#include "Python.h"
#include "Numeric/arrayobject.h"
#include "ranlib.h"
#include "stdio.h"

static PyObject *ErrorObject;

/* ----------------------------------------------------- */

static PyObject*
get_continuous_random(int num_dist_params, PyObject* self, PyObject* args, void* fun) {
  PyArrayObject *op;
  double *out_ptr;
  int i, n=-1;
  float a, b, c;

  switch(num_dist_params) {
  case 0:
    if( !PyArg_ParseTuple(args, "|i", &n) )
      return NULL;
    break;
  case 1:
    if( !PyArg_ParseTuple(args, "f|i", &a, &n) )
      return NULL;
    break;
  case 2:
    if( !PyArg_ParseTuple(args, "ff|i", &a, &b, &n) )
      return NULL;
    break;
  case 3:
    if( !PyArg_ParseTuple(args, "fff|i", &a, &b, &c, &n) )
      return NULL;
    break;
  }
  if( n == -1 )
    n = 1;

  /* Create a 1 dimensional array of length n of type double */
  op = (PyArrayObject*) PyArray_FromDims(1, &n, PyArray_DOUBLE);
  if( op == NULL )
    return NULL;

  out_ptr = (double *) op->data;
  for(i=0; i<n; i++) {
    switch(num_dist_params) {
    case 0:
      *out_ptr = (double) ((float (*)(void)) fun)();
      break;
    case 1:
      *out_ptr = (double) ((float (*)(float)) fun)(a);
      break;
    case 2:
      *out_ptr = (double) ((float (*)(float, float)) fun)(a,b);
      break;
    case 3:
      *out_ptr = (double) ((float (*)(float, float, float)) fun)(a,b,c);
      break;
    }
    out_ptr++;
  }

  return PyArray_Return(op);
}


static PyObject*
get_discrete_scalar_random(int num_integer_args, PyObject* self, PyObject* args, void* fun) {
  long int_arg;
  int n=-1, i;
  long* out_ptr;
  PyArrayObject* op;
  float float_arg;

  switch( num_integer_args ) {
  case 0:
    if( !PyArg_ParseTuple(args, "f|i", &float_arg, &n) ) {
      return NULL;
    }
    break;
  case 1:
    if( !PyArg_ParseTuple(args, "lf|i", &int_arg, &float_arg, &n) ) {
      return NULL;
    }
    break;
  }
  if( n==-1 ) {
    n = 1;
  }
  
  /* Create a 1 dimensional array of length n of type long */
  op = (PyArrayObject*) PyArray_FromDims(1, &n, PyArray_LONG);
  if( op == NULL ) {
    return NULL;
  }
  
  out_ptr = (long*) op->data;
  for(i=0; i<n; i++) {
    switch( num_integer_args ) {
    case 0:
      *out_ptr = ((long (*)(float)) fun)(float_arg);
      break;
    case 1:
      *out_ptr = ((long (*)(long, float)) fun)(int_arg, float_arg);
      break;
    }
    out_ptr++;
  }

  return PyArray_Return(op);
}


static char random_sample__doc__[] ="";

static PyObject *
random_sample(PyObject *self, PyObject *args) {
  return get_continuous_random(0, self, args, ranf);
}


static char standard_normal__doc__[] ="";

static PyObject *
standard_normal(PyObject *self, PyObject *args) {
  return get_continuous_random(0, self, args, snorm);
}


static char beta__doc__[] ="";

static PyObject *
beta(PyObject *self, PyObject *args) {
  return get_continuous_random(2, self, args, genbet);
}


static char exponential__doc__[] = "";

static PyObject *
exponential(PyObject *self, PyObject *args) { 
  return get_continuous_random(1, self, args, genexp);
}

static char normal__doc__[] = "";

static PyObject *
normal(PyObject *self, PyObject *args) { 
  return get_continuous_random(2, self, args, gennor);
}

static char uniform__doc__[] = "";

static PyObject *
uniform(PyObject *self, PyObject *args) { 
  return get_continuous_random(2, self, args, genunf);
}


static char standard_exp__doc__[] = "";

static PyObject *
standard_exp(PyObject *self, PyObject *args) { 
  return get_continuous_random(0, self, args, sexpo);
}


static char standard_gamma__doc__[] = "";

static PyObject *
standard_gamma(PyObject *self, PyObject *args) { 
  return get_continuous_random(1, self, args, sgamma);
}

static char gamma__doc__[] ="";

static PyObject *
/* there is a function named `gamma' in some libm's */
_gamma(PyObject *self, PyObject *args) { 
  return get_continuous_random(2, self, args, gengam);
}


static char f__doc__[] ="";

static PyObject *
f(PyObject *self, PyObject *args) {
  return get_continuous_random(2, self, args, genf);
}



static char noncentral_f__doc__[] ="";

static PyObject *
noncentral_f(PyObject *self, PyObject *args) {
  return get_continuous_random(3, self, args, gennf);
}


static char noncentral_chisquare__doc__[] ="";

static PyObject *
noncentral_chisquare(PyObject *self, PyObject *args) {
  return get_continuous_random(2, self, args, gennch);
}


static char chisquare__doc__[] ="";

static PyObject *
chisquare(PyObject *self, PyObject *args) {
  return get_continuous_random(1, self, args, genchi);
}


static char binomial__doc__[] ="";

static PyObject *
binomial(PyObject *self, PyObject *args) {
  return get_discrete_scalar_random(1, self, args, ignbin);
}


static char negative_binomial__doc__[] ="";

static PyObject *
negative_binomial(PyObject *self, PyObject *args) {
  return get_discrete_scalar_random(1, self, args, ignnbn);
}

static char poisson__doc__[] ="";

static PyObject *
poisson(PyObject *self, PyObject *args) {
  return get_discrete_scalar_random(0, self, args, ignpoi);
}


static char multinomial__doc__[] ="";

static PyObject*
multinomial(PyObject* self, PyObject* args) {
  int n=-1, i;
  long num_trials, num_categories;
  char* out_ptr;
  PyArrayObject* priors_array;
  PyObject* priors_object;
  PyArrayObject* op;
  int out_dimensions[2];

  if( !PyArg_ParseTuple(args, "lO|i", &num_trials, &priors_object, &n) ) {
    return NULL;
  }
  priors_array = (PyArrayObject*) PyArray_ContiguousFromObject(priors_object, PyArray_FLOAT, 1, 1);
  if( priors_array == NULL ) {
    return NULL;
  }
  num_categories = priors_array->dimensions[0]+1;
  if( n==-1 ) {
    n = 1;
  }
  
  /* Create an n by num_categories array of long */
  out_dimensions[0] = n;
  out_dimensions[1] = num_categories;
  op = (PyArrayObject*) PyArray_FromDims(2, out_dimensions, PyArray_LONG);
  if( op == NULL ) {
    return NULL;
  }
  
  out_ptr = op->data;
  for(i=0; i<n; i++) {
    genmul(num_trials, (float*)(priors_array->data), num_categories, (long*) out_ptr);
    out_ptr += op->strides[0];
  }

  return PyArray_Return(op);
}


static char permutation__doc__[] ="";

static PyObject*
permutation(PyObject* self, PyObject* args) {
  PyArrayObject* vector_array;
  PyObject* vector;
  int len;

  if( !PyArg_ParseTuple(args, "O", &vector) ) {
    return NULL;
  }
  vector_array = (PyArrayObject*) PyArray_CopyFromObject(vector, PyArray_LONG, 1, 0);
  if( vector_array == NULL ) {
    return NULL;
  }

  len = _PyArray_multiply_list(vector_array->dimensions, vector_array->nd);
  genprm((long *) vector_array->data, (long) len);
  
  return PyArray_Return(vector_array);
}

static char multivariate_normal__doc__[] ="";

static PyObject*
multivariate_normal(PyObject* self, PyObject* args) {
  PyArrayObject *amean=NULL, *acov=NULL, *ameanf=NULL, *acovf=NULL;
  PyObject *omean=NULL, *ocov=NULL;
  PyArrayObject *outarr=NULL, *outarr_f=NULL;
  float *parm=NULL, *work=NULL, *outptr;
  long p;
  int n=1, dims[2], k;

  if( !PyArg_ParseTuple(args, "OO|i", &omean, &ocov, &n) ) {
    return NULL;
  }
  amean = (PyArrayObject*) PyArray_ContiguousFromObject(omean, PyArray_DOUBLE, 1, 1);
  if( amean == NULL ) goto fail;

  ameanf = (PyArrayObject *)PyArray_Cast(amean, PyArray_FLOAT);
  if (ameanf == NULL) goto fail;

  acov = (PyArrayObject*) PyArray_ContiguousFromObject(ocov, PyArray_DOUBLE, 2, 2);
  if( acov == NULL ) goto fail;

  acovf = (PyArrayObject *)PyArray_Cast(acov, PyArray_FLOAT);
  if (acovf == NULL) goto fail;

  p = (long) ameanf->dimensions[0]; /* length of mean vector */
  if ((((long) acovf->dimensions[0]) != p) ||
      (((long) acovf->dimensions[1]) != p)) {
    fprintf(stderr, "%d %d %d %d", ameanf->dimensions[0], acovf->dimensions[0], acovf->dimensions[1], acovf->nd);
    PyErr_SetString(PyExc_ValueError, "Covariance matrix must be square with dimension size equal to the length\n of the mean vector.");
    goto fail;
  }

  work = PyMem_New(float, ((int) (p*(p+3)/2 + 1 + p)));
  if (work == NULL) {
    PyErr_SetString(PyExc_MemoryError, "Could not allocate needed memory.");
    goto fail;
  }
  parm = work + p;
  
  setgmn((float *)ameanf->data, (float *)acovf->data, p, parm);

  if (PyErr_Occurred()) goto fail;

  if (n == 1) {
    int tmp = p; /* avoid "invalid pointer type" warnings from compiler */
    outarr_f = (PyArrayObject *)PyArray_FromDims(1, &tmp, PyArray_FLOAT);
  }
  else {
    dims[0] = n;
    dims[1] = p;
    outarr_f = (PyArrayObject *)PyArray_FromDims(2, dims, PyArray_FLOAT);
  }
  if (outarr_f == NULL) goto fail;

  outptr = (float *) outarr_f->data;
  
  for (k=0; k<n; k++) {
    genmn(parm, outptr, work);
    outptr += p;
  }

  outarr = (PyArrayObject *)PyArray_Cast(outarr_f, PyArray_DOUBLE);
  if (outarr == NULL) goto fail;
    
  PyMem_Del(work);
  Py_DECREF(amean);
  Py_DECREF(acov);
  Py_DECREF(ameanf);
  Py_DECREF(acovf);
  return PyArray_Return(outarr);

 fail:
  if (work != NULL) PyMem_Del(work);
  Py_XDECREF(amean);
  Py_XDECREF(ameanf);
  Py_XDECREF(acov);
  Py_XDECREF(acovf);
  Py_XDECREF(outarr_f);
  Py_XDECREF(outarr);
  return NULL;
}


static PyObject *
random_set_seeds(PyObject *self, PyObject *args)
{
  long seed1, seed2;

  if (!PyArg_ParseTuple(args, "ll", &seed1, &seed2)) return NULL;


  setall(seed1, seed2);
  if (PyErr_Occurred ()) return NULL;
  Py_INCREF(Py_None);
  return (PyObject *)Py_None;
}


static PyObject *
random_get_seeds(PyObject *self, PyObject *args)
{
  long seed1, seed2;

  if (!PyArg_ParseTuple(args, "")) return NULL;

  getsd(&seed1, &seed2);

  return Py_BuildValue("ll", seed1, seed2);
}


/* Missing interfaces to */
/* multivariate normal (genmn), 
   permutation (genprm),
*/

/* List of methods defined in the module */

static struct PyMethodDef random_methods[] = {
 {"sample",     random_sample,          1,      random_sample__doc__},
 {"standard_normal", standard_normal,   1,      standard_normal__doc__},
 {"beta",	beta,                   1,      beta__doc__},
 {"exponential", exponential,           1,      exponential__doc__},
 {"gamma",	_gamma,                  1,      gamma__doc__},
 {"f",	        f,                      1,      f__doc__},
 {"noncentral_f", noncentral_f,         1,      noncentral_f__doc__},
 {"chi2",	chisquare,              1,      chisquare__doc__},
 {"noncentral_chi2", noncentral_chisquare,
                                        1,      noncentral_chisquare__doc__},
 {"binomial",	binomial,               1,      binomial__doc__},
 {"negative_binomial",  negative_binomial,
                                        1,      negative_binomial__doc__},
 {"multinomial", multinomial,           1,      multinomial__doc__},
 {"multivariate_normal", multivariate_normal,         1,      multivariate_normal__doc__},
 {"normal", normal,                     1,      normal__doc__},
 {"permutation", permutation,           1,      permutation__doc__},
 {"poisson",    poisson,                1,      poisson__doc__},
 {"uniform",    uniform,                1,      uniform__doc__},
 {"standard_exp",      standard_exp,         1,      standard_exp__doc__},
 {"standard_gamma",    standard_gamma,     1,      standard_gamma__doc__},
 {"set_seeds",  random_set_seeds,       1, },
 {"get_seeds",  random_get_seeds,       1, },
 {NULL,		NULL}		/* sentinel */
};


/* Initialization function for the module (*must* be called initranlib) */

static char random_module_documentation[] = 
""
;

DL_EXPORT(void)
initrand(void)
{
	PyObject *m, *d;

	/* Create the module and add the functions */
	m = Py_InitModule4("rand", random_methods,
		random_module_documentation,
		(PyObject*)NULL,PYTHON_API_VERSION);

	/* Import the array object */
	import_array();

	/* Add some symbolic constants to the module */
	d = PyModule_GetDict(m);
	ErrorObject = PyString_FromString("rand.error");
	PyDict_SetItemString(d, "error", ErrorObject);

	/* XXXX Add constants here */
	
	/* Check for errors */
	if (PyErr_Occurred())
		Py_FatalError("can't initialize module rand");
}
