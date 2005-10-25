#define JDSYM_OP_UNSYM 1
#define JDSYM_OP_SYM   2

typedef struct CorrEqObject {
  PyObject_VAR_HEAD
  int op_type;
  int n;
  int k;
  PyObject *amat;
  PyObject *mmat;
  PyObject *prec;
  double theta;
  double *Q;
  double *Qm;
  double *Y;
  /* ldQ and ldQm are equal to n */
  int *Hpiv;
  double *Hlu;
  int ldh;
  double *work;
  void (*right)(struct CorrEqObject *, double *);
  void (*update)(struct CorrEqObject *, int, double);
} CorrEqObject;


static void 
CorrEq_A_theta_M(CorrEqObject *self, int n, double *x, double *y, double *w) {

  double alpha = -self->theta;

  SpMatrix_Matvec(self->amat, n, x, n, y);	   /* y := A*x; */
  if (alpha != 0.0) {
    SpMatrix_Matvec(self->mmat, n, x, n, w);	   /* work := M*x */
    F77(daxpy)(&n, &alpha, w, &ONE, y, &ONE);      /* y := -theta*work + y */
  }
}

static void 
CorrEq_A_theta_I(CorrEqObject *self, int n, double *x, double *y) {

  double alpha = -self->theta;

  SpMatrix_Matvec(self->amat, n, x, n, y);	   /* y := A*x; */
  F77(daxpy)(&n, &alpha, x, &ONE, y, &ONE);        /* y := -theta*x + y*/
}

/* PROJECT1 - compute y := (I - A*B')*x
 *
 * sizes of matrices and vectors:
 *   A and B     --  r-by-c matrices, ldim = r
 *   x, y and w  --  r vectors
 *
 * (x, y) and also (A, B) may point to the same memory locations
 */
static void 
CorrEq_project1(int r, int c, double *A, double *B, double *x, double *y, double *w){
  /* if x and y do not refer to the same vector, we copy x to y */
  if (x != y)
    F77(dcopy)(&r, x, &ONE, y, &ONE);
  /* w = B'*y */
  F77(dgemv)("t", &r, &c, &DONE, B, &r, y, &ONE, &DZER, w, &ONE, 1);
  /* y = y - A*w */
  F77(dgemv)("n", &r, &c, &DMONE, A, &r, w, &ONE, &DONE, y, &ONE, 1);
}

/* PROJECT2 - compute y := (I - Y*inv(H)*Qm') * x
 *                     or
 *                    y := (I - Y*inv(H)*Q' ) * x
 *
 *            depending on parameter Qx                  
 */
static void 
CorrEq_project2(CorrEqObject *self, double *Qx, double *x, double *y, double *w){
  int info;
  int n = self->n;
  int k = self->k;

  /* if x and y do not refer to the same vector, we copy x to y */
  if (x != y)
    F77(dcopy)(&n, x, &ONE, y, &ONE);
  /* w = Qx'*x */
  F77(dgemv)("t", &n, &k, &DONE, Qx, &n,
	     x, &ONE, &DZER, w, &ONE, 1);
  /* w = inv(H)*w */
  F77(dgetrs)("n", &k, &ONE, self->Hlu, &(self->ldh),
	      self->Hpiv, w, &n, &info, 1);
  assert(info == 0);
  /* y = y - Y*w */
  F77(dgemv)("n", &n, &k, &DMONE, self->Y, &n,
	     w, &ONE, &DONE, y, &ONE, 1);
}

static void
CorrEq_update(CorrEqObject *self, int k, double theta) {
  self->k = k;
  self->theta = theta;
}

static void
CorrEq_right(CorrEqObject *self, double *r) {
  int n = self->n;
  int k = self->k;

  if (self->op_type == JDSYM_OP_SYM)
    if (self->mmat) {
      /* r = r - Qb*(Q'*r); */
      CorrEq_project1(n, k, self->Qm, self->Q, r, r, self->work); 
    } else {
      /* r :=  (I - Q*Q')*r */
      mgs(r, n, k, self->Q);
    }
  else /* op_type = JDSYM_OP_UNSYM */
    if (self->mmat)
      if (self->prec) {
	/* mmat && prec */
	/* r := (I - Y'*inv(H)*Qb')*inv(K)*r */
	SpMatrix_Precon(self->prec, n, r, self->work);
	CorrEq_project2(self, self->Qm, self->work, r, self->work + n);
      } else {
	/* mmat && !prec */
	/* r :=  (I - Q*Qb')*r */
	CorrEq_project1(n, k, self->Q, self->Qm, r, r, self->work);
      }
    else
      if (self->prec) {
	/* !mmat && prec */
	/* r := (I - Y'*inv(H)*Q')*inv(K)*r */
	SpMatrix_Precon(self->prec, n, r, self->work);
	CorrEq_project2(self, self->Q, self->work, r, self->work + n);
      } else {
	/* !mmat && !prec */
	/* r :=  (I - Q*Q')*r */
	mgs(r, n, k, self->Q);
      }
}

static PyObject *
CorrEq_matvec(CorrEqObject *self, PyObject *args) {
  PyArrayObject *xp, *yp;
  double *x, *y;
  int n = self->n;
  int k = self->k;

  /* parse input arguments */
  SPMATRIX_PARSE_ARGS_ARR_ARR(args, xp, yp, self->n, self->n);
  x = (double *)(xp->data);
  y = (double *)(yp->data);

  if (self->op_type == JDSYM_OP_SYM) {
    if (self->mmat) {
      /* y := (I - Qm*Q')*(A - theta*M)*x */
      CorrEq_A_theta_M(self, n, x, y, self->work);
      CorrEq_project1(n, k, self->Qm, self->Q, y, y, self->work);
    } else {
      /* y := (I - Q*Q')*(A - theta*I)*x */
      CorrEq_A_theta_I(self, n, x, y);
      CorrEq_project1(n, k, self->Q, self->Q, y, y, self->work);
    }
  } else { /* self->op_type = JDSYM_OP_UNSYM */
    if (self->mmat)
      if (self->prec) {
	/* mmat && prec */
	/* y := (I - Y'*inv(H)*QM')*inv(K)*(A-theta*M)*x */
	CorrEq_A_theta_M(self, n, x, self->work + n, self->work);
	SpMatrix_Precon(self->prec, n, self->work + n, y);
	CorrEq_project2(self, self->Qm, y, y, self->work);
      } else { 
	/* mmat && !prec */
	/* y := (I - Q*Qm')*(A-theta*M)*x */
	CorrEq_A_theta_M(self, n, x, y, self->work);
	CorrEq_project1(n, k, self->Q, self->Qm, y, y, self->work);
      }
    else /* ! self->mmat */
      if (self->prec) {
	/* !mmat && prec */
	/* y := (I - Y'*inv(H)*Q')*inv(K)*(A-theta*I)*x */
	CorrEq_A_theta_I(self, n, x, self->work + n);
	SpMatrix_Precon(self->prec, n, self->work + n, y);
	CorrEq_project2(self, self->Q, y, y, self->work);
      } else { 
	/* !mmat && !prec */
	/* y := (I - Q*Q')*(A-theta*I)*x */
	CorrEq_A_theta_I(self, n, x, y);
	CorrEq_project1(n, k, self->Q, self->Q, y, y, self->work);
      }
  }

  /* return Py_None */
  Py_INCREF(Py_None); 
  return Py_None;
}

static PyObject *
CorrEq_precon(CorrEqObject *self, PyObject *args) {
  PyArrayObject *xp, *yp;
  double *x, *y;
  int n = self->n;
  int k = self->k;

  /* parse input arguments */
  SPMATRIX_PARSE_ARGS_ARR_ARR(args, xp, yp, self->n, self->n);
  x = (double *)(xp->data);
  y = (double *)(yp->data);

  if (self->op_type == JDSYM_OP_SYM)
    if (self->mmat)
      if (self->prec) {
	/* mmat && prec */
	/* y = (I-Y*H\Qm')*inv(K)*x */
	SpMatrix_Precon(self->prec, n, x, y);
	CorrEq_project2(self, self->Qm, y, y, self->work);
      } else {
	/* mmat && !prec */
	/* y = (I - Q*Qm')*x */
	CorrEq_project1(n, k, self->Q, self->Qm, x, y, self->work);
      }
    else
      if (self->prec) {
	/* !mmat && prec */
	/* y = (I-Y*H\Q')*inv(K)*x */
	SpMatrix_Precon(self->prec, n, x, y);
	CorrEq_project2(self, self->Q, y, y, self->work);
      } else {
	/* !mmat && !prec */
	/* y = x */
	F77(dcopy)(&n, x, &ONE, y, &ONE);
      }
  else { /* self->optype == JDSYM_OP_UNSYM */
    /* the Preconditioner is contained in the operator */
    F77(dcopy)(&n, x, &ONE, y, &ONE);
  }

  /* return Py_None */
  Py_INCREF(Py_None); 
  return Py_None;
}

/** table of object methods
 */
PyMethodDef CorrEq_methods[] = {
  {"matvec", (PyCFunction)CorrEq_matvec, METH_VARARGS},
  {"precon", (PyCFunction)CorrEq_precon, METH_VARARGS},
  {NULL, NULL}			/* sentinel */
};

/*********************************************************************** 
 * CorrEqType methods
 */

static void
CorrEq_dealloc(CorrEqObject *self)
{
  PyMem_DEL(self->work);
  PyObject_Del(self);
}

static PyObject *
CorrEq_getattr(CorrEqObject *self, char *name)
{
  if (strcmp(name, "shape") == 0)
    return Py_BuildValue("(i,i)", self->n, self->n);
  if (strcmp(name, "__members__") == 0) {
    char *members[] = {"shape"};
    int i;

    PyObject *list = PyList_New(sizeof(members)/sizeof(char *));
    if (list != NULL) {
      for (i = 0; i < sizeof(members)/sizeof(char *); i ++)
	PyList_SetItem(list, i, PyString_FromString(members[i]));
      if (PyErr_Occurred()) {
	Py_DECREF(list);
	list = NULL;
      }
    }
    return list;
  }
  return Py_FindMethod(CorrEq_methods, (PyObject *)self, name);
}

/***********************************************************************
 * Type structures
 */

PyTypeObject CorrEqType = {
  PyObject_HEAD_INIT(NULL)
  0,
  "CorrEqSystem",
  sizeof(CorrEqObject),
  0,
  (destructor)CorrEq_dealloc,   /* tp_dealloc */
  0,				/* tp_print */
  (getattrfunc)CorrEq_getattr,  /* tp_getattr */
  0,				/* tp_setattr */
  0,				/* tp_compare */
  0,				/* tp_repr */
  0,				/* tp_as_number*/
  0,				/* tp_as_sequence*/
  0,				/* tp_as_mapping*/
  0,				/* tp_hash */
};

/*********************************************************************** 
 * Object construction functions
 */

PyObject *
newCorrEqObject(int op_type, int n, 
		PyObject *amat, PyObject *mmat, PyObject *prec,
		double *Q, double *Qm, double *Y, int *Hpiv, double *Hlu, int ldh) {
  CorrEqObject *self;
  
  /* create new CorrEqObject */
  self = PyObject_New(CorrEqObject, &CorrEqType);
  if (self == NULL)
    return PyErr_NoMemory();

  /* setup work arrays */
  self->work = PyMem_New(double, 2*n);
  if (self->work == NULL) {
    PyObject_Del(self);
    return PyErr_NoMemory();
  }

  /* set function pointers */
  self->right = CorrEq_right;
  self->update = CorrEq_update;

  self->op_type = op_type;
  self->n = n;
  self->amat = amat;
  self->mmat = mmat;
  self->prec = prec;
  self->Q = Q;
  self->Qm = Qm;
  self->Y = Y;
  self->Hpiv = Hpiv;
  self->Hlu = Hlu;
  self->ldh = ldh;
  
  return (PyObject *)self;  
}
