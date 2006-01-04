#include "Python.h"
#include "numpy/arrayobject.h"

#if defined(NO_APPEND_FORTRAN)
#if defined(UPPERCASE_FORTRAN)
#define F_FUNC(f,F) F
#else
#define F_FUNC(f,F) f
#endif
#else
#if defined(UPPERCASE_FORTRAN)
#define F_FUNC(f,F) F##_
#else
#define F_FUNC(f,F) f##_
#endif
#endif

#define PYERR(errobj,message) {PyErr_SetString(errobj,message); goto fail;}
#define PYERR2(errobj,message) {PyErr_Print(); PyErr_SetString(errobj, message); goto fail;}
#define ISCONTIGUOUS(m) ((m)->flags & CONTIGUOUS)

#define MAX(n1,n2) ((n1) > (n2))?(n1):(n2);
#define MIN(n1,n2) ((n1) > (n2))?(n2):(n1);

struct ODR_info_ {
  PyObject* fcn;
  PyObject* fjacb;
  PyObject* fjacd;
  PyObject* pyBeta;
  PyObject* extra_args;
};

typedef struct ODR_info_ ODR_info;

static ODR_info odr_global;

static PyObject *odr_error=NULL;
static PyObject *odr_stop=NULL;

void fcn_callback(int *n, int *m, int *np, int *nq, int *ldn, int *ldm,
		  int *ldnp, double *beta, double *xplusd, int *ifixb,
		  int *ifixx, int *ldfix, int *ideval, double *f,
		  double *fjacb, double *fjacd, int *istop);

PyObject *gen_output(int n, int m, int np, int nq, int ldwe, int ld2we,
		     PyArrayObject *beta, PyArrayObject *work, PyArrayObject *iwork,
		     int isodr, int info, int full_output);

PyObject *odr(PyObject *self, PyObject *args, PyObject *kwds);

#define PyArray_CONTIGUOUS(m) (ISCONTIGUOUS(m) ? Py_INCREF(m), m : \
(PyArrayObject *)(PyArray_ContiguousFromObject((PyObject *)(m), \
(m)->descr->type_num, 0,0)))
#define D(dbg) printf("we're here: %i\n", dbg)
#define EXIST(name,obj) if (obj==NULL){printf("%s\n",name);}
static void check_args(int n, int m, int np, int nq, 
		       PyArrayObject *beta,
		       PyArrayObject *y, int ldy, 
		       PyArrayObject *x, int ldx, 
		       PyArrayObject *we, int ldwe, int ld2we, 
		       PyArrayObject *wd, int ldwd, int ld2wd, 
		       PyArrayObject *ifixb, PyArrayObject *ifixx, int ldifx,
		       int job, int ndigit, double taufac, double sstol, 
		       double partol, int maxit,
		       PyArrayObject *stpb, PyArrayObject *stpd, int ldstpd,
		       PyArrayObject *sclb, PyArrayObject *scld, int ldscld,
		       PyArrayObject *work, int lwork, 
		       PyArrayObject *iwork, int liwork,
		       int info);
