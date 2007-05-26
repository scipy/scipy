/*
  Python-C wrapper of FITPACK (by P. Dierckx) (in netlib known as dierckx)
  Author: Pearu Peterson <pearu@ioc.ee>
  June 1.-4., 1999
  June 7. 1999
  $Revision$
  $Date$
 */

/*  module_methods:
  {"_curfit", fitpack_curfit, METH_VARARGS, doc_curfit},
  {"_spl_", fitpack_spl_, METH_VARARGS, doc_spl_},
  {"_splint", fitpack_splint, METH_VARARGS, doc_splint},
  {"_sproot", fitpack_sproot, METH_VARARGS, doc_sproot},
  {"_spalde", fitpack_spalde, METH_VARARGS, doc_spalde},
  {"_parcur", fitpack_parcur, METH_VARARGS, doc_parcur},
  {"_surfit", fitpack_surfit, METH_VARARGS, doc_surfit},
  {"_bispev", fitpack_bispev, METH_VARARGS, doc_bispev},
  {"_insert", fitpack_insert, METH_VARARGS, doc_insert},
 */
/* link libraries: (one item per line)
   ddierckx
 */
/* python files: (to be imported to Multipack.py)
   fitpack.py
 */
#if defined(NO_APPEND_FORTRAN)
#define CURFIT curfit
#define PERCUR percur
#define SPALDE spalde
#define SPLDER splder
#define SPLEV splev
#define SPLINT splint
#define SPROOT sproot
#define PARCUR parcur
#define CLOCUR clocur
#define SURFIT surfit
#define BISPEV bispev
#define PARDER parder
#define INSERT insert
#else
#define CURFIT curfit_
#define PERCUR percur_
#define SPALDE spalde_
#define SPLDER splder_
#define SPLEV splev_
#define SPLINT splint_
#define SPROOT sproot_
#define PARCUR parcur_
#define CLOCUR clocur_
#define SURFIT surfit_
#define BISPEV bispev_
#define PARDER parder_
#define INSERT insert_
#endif
void CURFIT(int*,int*,double*,double*,double*,double*,double*,int*,double*,int*,int*,double*,double*,double*,double*,int*,int*,int*);
void PERCUR(int*,int*,double*,double*,double*,int*,double*,int*,int*,double*,double*,double*,double*,int*,int*,int*);
void SPALDE(double*,int*,double*,int*,double*,double*,int*);
void SPLDER(double*,int*,double*,int*,int*,double*,double*,int*,double*,int*);
void SPLEV(double*,int*,double*,int*,double*,double*,int*,int*);
double SPLINT(double*,int*,double*,int*,double*,double*,double*);
void SPROOT(double*,int*,double*,double*,int*,int*,int*);
void PARCUR(int*,int*,int*,int*,double*,int*,double*,double*,double*,double*,int*,double*,int*,int*,double*,int*,double*,double*,double*,int*,int*,int*);
void CLOCUR(int*,int*,int*,int*,double*,int*,double*,double*,int*,double*,int*,int*,double*,int*,double*,double*,double*,int*,int*,int*);
void SURFIT(int*,int*,double*,double*,double*,double*,double*,double*,double*,double*,int*,int*,double*,int*,int*,int*,double*,int*,double*,int*,double*,double*,double*,double*,int*,double*,int*,int*,int*,int*);
void BISPEV(double*,int*,double*,int*,double*,int*,int*,double*,int*,double*,int*,double*,double*,int*,int*,int*,int*);
void PARDER(double*,int*,double*,int*,double*,int*,int*,int*,int*,double*,int*,double*,int*,double*,double*,int*,int*,int*,int*);
void INSERT(int*,double*,int*,double*,int*,double*,double*,int*,double*,int*,int*);

/* Note that curev, cualde need no interface. */

static char doc_bispev[] = " [z,ier] = _bispev(tx,ty,c,kx,ky,x,y,nux,nuy)";
static PyObject *fitpack_bispev(PyObject *dummy, PyObject *args) {
  int nx,ny,kx,ky,mx,my,lwrk,*iwrk,kwrk,ier,lwa,mxy,nux,nuy;
  double *tx,*ty,*c,*x,*y,*z,*wrk,*wa = NULL;
  PyArrayObject *ap_x = NULL,*ap_y = NULL,*ap_z = NULL,*ap_tx = NULL,\
    *ap_ty = NULL,*ap_c = NULL;
  PyObject *x_py = NULL,*y_py = NULL,*c_py = NULL,*tx_py = NULL,*ty_py = NULL;
  if (!PyArg_ParseTuple(args, "OOOiiOOii",&tx_py,&ty_py,&c_py,&kx,&ky,
			&x_py,&y_py,&nux,&nuy))
    return NULL;
  ap_x = (PyArrayObject *)PyArray_ContiguousFromObject(x_py, PyArray_DOUBLE, 0, 1);
  ap_y = (PyArrayObject *)PyArray_ContiguousFromObject(y_py, PyArray_DOUBLE, 0, 1);
  ap_c = (PyArrayObject *)PyArray_ContiguousFromObject(c_py, PyArray_DOUBLE, 0, 1);
  ap_tx = (PyArrayObject *)PyArray_ContiguousFromObject(tx_py, PyArray_DOUBLE, 0, 1);
  ap_ty = (PyArrayObject *)PyArray_ContiguousFromObject(ty_py, PyArray_DOUBLE, 0, 1);
  if (ap_x == NULL || ap_y == NULL || ap_c == NULL || ap_tx == NULL \
      || ap_ty == NULL) goto fail;  
  x = (double *) ap_x->data;
  y = (double *) ap_y->data;
  c = (double *) ap_c->data;
  tx = (double *) ap_tx->data;
  ty = (double *) ap_ty->data;
  nx = ap_tx->dimensions[0];
  ny = ap_ty->dimensions[0];
  mx = ap_x->dimensions[0];
  my = ap_y->dimensions[0];
  mxy = mx*my;
  ap_z = (PyArrayObject *)PyArray_FromDims(1,&mxy,PyArray_DOUBLE);
  z = (double *) ap_z->data;
  if (nux || nuy) 
    lwrk = mx*(kx+1-nux)+my*(ky+1-nuy)+(nx-kx-1)*(ny-ky-1);
  else
    lwrk = mx*(kx+1)+my*(ky+1);
  kwrk = mx+my;
  lwa = lwrk+kwrk;
  if ((wa = (double *)malloc(lwa*sizeof(double)))==NULL) {
    PyErr_NoMemory();
    goto fail;
  }
  wrk = wa;
  iwrk = (int *)(wrk+lwrk);
  if (nux || nuy)
    PARDER(tx,&nx,ty,&ny,c,&kx,&ky,&nux,&nuy,x,&mx,y,&my,z,wrk,&lwrk,iwrk,&kwrk,&ier);
  else
    BISPEV(tx,&nx,ty,&ny,c,&kx,&ky,x,&mx,y,&my,z,wrk,&lwrk,iwrk,&kwrk,&ier);

  if (wa) free(wa);
  Py_DECREF(ap_x);
  Py_DECREF(ap_y);
  Py_DECREF(ap_c);
  Py_DECREF(ap_tx);
  Py_DECREF(ap_ty);
  return Py_BuildValue("Ni",PyArray_Return(ap_z),ier);
  fail:
  if (wa) free(wa);
  Py_XDECREF(ap_x);
  Py_XDECREF(ap_y);
  Py_XDECREF(ap_z);
  Py_XDECREF(ap_c);
  Py_XDECREF(ap_tx);
  Py_XDECREF(ap_ty);
  return NULL;
}

static char doc_surfit[] = " [tx,ty,c,o] = _surfit(x,y,z,w,xb,xe,yb,ye,kx,ky,iopt,s,eps,tx,ty,nxest,nyest,wrk,lwrk1,lwrk2)";
static PyObject *fitpack_surfit(PyObject *dummy, PyObject *args) {
  int iopt,m,kx,ky,nxest,nyest,nx,ny,lwrk1,lwrk2,*iwrk,kwrk,ier,lwa,nxo,nyo,\
    i,lc,lcest,nmax;
  double *x,*y,*z,*w,xb,xe,yb,ye,s,*tx,*ty,*c,fp,*wrk1,*wrk2,*wa = NULL,eps;
  PyArrayObject *ap_x = NULL,*ap_y = NULL,*ap_z,*ap_w = NULL,\
    *ap_tx = NULL,*ap_ty = NULL,*ap_c = NULL;
  PyArrayObject *ap_wrk = NULL;
  PyObject *x_py = NULL,*y_py = NULL,*z_py = NULL,*w_py = NULL,\
    *tx_py = NULL,*ty_py = NULL;
  PyObject *wrk_py=NULL;
  nx=ny=ier=nxo=nyo=0;
  if (!PyArg_ParseTuple(args, "OOOOddddiiiddOOiiOii",\
			&x_py,&y_py,&z_py,&w_py,&xb,&xe,\
			&yb,&ye,&kx,&ky,&iopt,&s,&eps,&tx_py,&ty_py,&nxest,&nyest,\
			&wrk_py,&lwrk1,&lwrk2)) return NULL;
  ap_x = (PyArrayObject *)PyArray_ContiguousFromObject(x_py, PyArray_DOUBLE, 0, 1);
  ap_y = (PyArrayObject *)PyArray_ContiguousFromObject(y_py, PyArray_DOUBLE, 0, 1);
  ap_z = (PyArrayObject *)PyArray_ContiguousFromObject(z_py, PyArray_DOUBLE, 0, 1);
  ap_w = (PyArrayObject *)PyArray_ContiguousFromObject(w_py, PyArray_DOUBLE, 0, 1);
  ap_wrk=(PyArrayObject *)PyArray_ContiguousFromObject(wrk_py, PyArray_DOUBLE, 0, 1);
  /*ap_iwrk=(PyArrayObject *)PyArray_ContiguousFromObject(iwrk_py, PyArray_INT, 0, 1);*/
  if (ap_x == NULL || ap_y == NULL || ap_z == NULL || ap_w == NULL \
      || ap_wrk == NULL) goto fail;
  x = (double *) ap_x->data;
  y = (double *) ap_y->data;
  z = (double *) ap_z->data;
  w = (double *) ap_w->data;
  m = ap_x->dimensions[0];
  nmax=nxest;
  if (nmax<nyest) nmax=nyest;
  lcest=(nxest-kx-1)*(nyest-ky-1);
  kwrk=m+(nxest-2*kx-1)*(nyest-2*ky-1);
  lwa = 2*nmax+lcest+lwrk1+lwrk2+kwrk;
  if ((wa = (double *)malloc(lwa*sizeof(double)))==NULL) {
    PyErr_NoMemory();
    goto fail;
  }
  tx = wa;
  ty = tx + nmax;
  c = ty + nmax;
  wrk1 = c + lcest;
  iwrk = (int *)(wrk1 + lwrk1);
  wrk2 = (double *)(iwrk+kwrk);
  if (iopt) {
    ap_tx=(PyArrayObject *)PyArray_ContiguousFromObject(tx_py, PyArray_DOUBLE, 0, 1);
    ap_ty=(PyArrayObject *)PyArray_ContiguousFromObject(ty_py, PyArray_DOUBLE, 0, 1);
    if (ap_tx == NULL || ap_ty == NULL) goto fail;
    nx = nxo = ap_tx->dimensions[0];
    ny = nyo = ap_ty->dimensions[0];
    memcpy(tx,ap_tx->data,nx*sizeof(double));
    memcpy(ty,ap_ty->data,ny*sizeof(double));
  }
  if (iopt==1) {
    lc = (nx-kx-1)*(ny-ky-1);
    memcpy(wrk1,ap_wrk->data,lc*sizeof(double));
    /*memcpy(iwrk,ap_iwrk->data,n*sizeof(int));*/
  }
  SURFIT(&iopt,&m,x,y,z,w,&xb,&xe,&yb,&ye,&kx,&ky,&s,&nxest,&nyest,&nmax,&eps,&nx,tx,&ny,ty,c,&fp,wrk1,&lwrk1,wrk2,&lwrk2,iwrk,&kwrk,&ier);
  i=0;
  while ((ier>10) && (i++<5)) {
    lwrk2=ier;
    if ((wrk2 = (double *)malloc(lwrk2*sizeof(double)))==NULL) {
      PyErr_NoMemory();
      goto fail;
    }
    SURFIT(&iopt,&m,x,y,z,w,&xb,&xe,&yb,&ye,&kx,&ky,&s,&nxest,&nyest,&nmax,&eps,&nx,tx,&ny,ty,c,&fp,wrk1,&lwrk1,wrk2,&lwrk2,iwrk,&kwrk,&ier);
    if (wrk2) free(wrk2);
  }
  if (ier==10) {
	  PyErr_SetString(PyExc_ValueError, "Invalid inputs.");
	  goto fail;
  }
  lc = (nx-kx-1)*(ny-ky-1);
  Py_XDECREF(ap_tx);
  Py_XDECREF(ap_ty);
  ap_tx = (PyArrayObject *)PyArray_FromDims(1,&nx,PyArray_DOUBLE);
  ap_ty = (PyArrayObject *)PyArray_FromDims(1,&ny,PyArray_DOUBLE);
  ap_c = (PyArrayObject *)PyArray_FromDims(1,&lc,PyArray_DOUBLE);
  if (ap_tx == NULL || ap_ty == NULL || ap_c == NULL) goto fail;
  if ((iopt==0)||(nx>nxo)||(ny>nyo)) {
    Py_XDECREF(ap_wrk);
    ap_wrk = (PyArrayObject *)PyArray_FromDims(1,&lc,PyArray_DOUBLE);
    if (ap_wrk == NULL) goto fail;
    /*ap_iwrk = (PyArrayObject *)PyArray_FromDims(1,&n,PyArray_INT);*/
  }
  if(ap_wrk->dimensions[0]<lc) {
    Py_XDECREF(ap_wrk);
    ap_wrk = (PyArrayObject *)PyArray_FromDims(1,&lc,PyArray_DOUBLE);
    if (ap_wrk == NULL) goto fail;
  }
  memcpy(ap_tx->data,tx,nx*sizeof(double));
  memcpy(ap_ty->data,ty,ny*sizeof(double));
  memcpy(ap_c->data,c,lc*sizeof(double));
  memcpy(ap_wrk->data,wrk1,lc*sizeof(double));
  /*memcpy(ap_iwrk->data,iwrk,n*sizeof(int));*/
  if (wa) free(wa);
  Py_DECREF(ap_x);
  Py_DECREF(ap_y);
  Py_DECREF(ap_z);
  Py_DECREF(ap_w);
  return Py_BuildValue("NNN{s:N,s:i,s:d}",PyArray_Return(ap_tx),\
		       PyArray_Return(ap_ty),PyArray_Return(ap_c),\
		       "wrk",PyArray_Return(ap_wrk),\
		       "ier",ier,"fp",fp);
  fail:
  if (wa) free(wa);
  Py_XDECREF(ap_x);
  Py_XDECREF(ap_y);
  Py_XDECREF(ap_z);
  Py_XDECREF(ap_w);
  Py_XDECREF(ap_tx);
  Py_XDECREF(ap_ty);
  Py_XDECREF(ap_wrk);
  /*Py_XDECREF(ap_iwrk);*/
  if (!PyErr_Occurred()) {
	  PyErr_SetString(PyExc_ValueError, "An error occurred.");
  }
  return NULL;
}


static char doc_parcur[] = " [t,c,o] = _parcur(x,w,u,ub,ue,k,iopt,ipar,s,t,nest,wrk,iwrk,per)";
static PyObject *fitpack_parcur(PyObject *dummy, PyObject *args) {
  int k,iopt,ipar,nest,*iwrk,idim,m,mx,n=0,no=0,nc,ier,lc,lwa,lwrk,i,per;
  double *x,*w,*u,*c,*t,*wrk,*wa=NULL,ub,ue,fp,s;
  PyObject *x_py = NULL,*u_py = NULL,*w_py = NULL,*t_py = NULL;
  PyObject *wrk_py=NULL,*iwrk_py=NULL;
  PyArrayObject *ap_x = NULL,*ap_u = NULL,*ap_w = NULL,*ap_t = NULL,*ap_c = NULL;
  PyArrayObject *ap_wrk = NULL,*ap_iwrk = NULL;
  if (!PyArg_ParseTuple(args, "OOOddiiidOiOOi",&x_py,&w_py,&u_py,&ub,&ue,\
			&k,&iopt,&ipar,&s,&t_py,&nest,&wrk_py,&iwrk_py,&per)) return NULL;
  ap_x = (PyArrayObject *)PyArray_ContiguousFromObject(x_py, PyArray_DOUBLE, 0, 1);
  ap_u = (PyArrayObject *)PyArray_ContiguousFromObject(u_py, PyArray_DOUBLE, 0, 1);
  ap_w = (PyArrayObject *)PyArray_ContiguousFromObject(w_py, PyArray_DOUBLE, 0, 1);
  ap_wrk=(PyArrayObject *)PyArray_ContiguousFromObject(wrk_py, PyArray_DOUBLE, 0, 1);
  ap_iwrk=(PyArrayObject *)PyArray_ContiguousFromObject(iwrk_py, PyArray_INT, 0, 1);
  if (ap_x == NULL || ap_u == NULL || ap_w == NULL || ap_wrk == NULL || ap_iwrk == NULL) goto fail;
  x = (double *) ap_x->data;
  u = (double *) ap_u->data;
  w = (double *) ap_w->data;
  m = ap_w->dimensions[0];
  mx = ap_x->dimensions[0];
  idim = mx/m;
  if (per)
    lwrk=m*(k+1)+nest*(7+idim+5*k);
  else
    lwrk=m*(k+1)+nest*(6+idim+3*k);
  nc=idim*nest;
  lwa = nc+2*nest+lwrk;
  if ((wa = (double *)malloc(lwa*sizeof(double)))==NULL) {
    PyErr_NoMemory();
    goto fail;
  }
  t = wa;
  c = t + nest;
  wrk = c + nc;
  iwrk = (int *)(wrk + lwrk);
  if (iopt) {
    ap_t=(PyArrayObject *)PyArray_ContiguousFromObject(t_py, PyArray_DOUBLE, 0, 1);
    if (ap_t == NULL) goto fail;
    n = no = ap_t->dimensions[0];
    memcpy(t,ap_t->data,n*sizeof(double));
  }
  if (iopt==1) {
    memcpy(wrk,ap_wrk->data,n*sizeof(double));
    memcpy(iwrk,ap_iwrk->data,n*sizeof(int));
  }
  if (per)
    CLOCUR(&iopt,&ipar,&idim,&m,u,&mx,x,w,&k,&s,&nest,&n,t,&nc,\
	   c,&fp,wrk,&lwrk,iwrk,&ier);
  else
    PARCUR(&iopt,&ipar,&idim,&m,u,&mx,x,w,&ub,&ue,&k,&s,&nest,&n,t,&nc,\
	   c,&fp,wrk,&lwrk,iwrk,&ier);
  if (ier==10) goto fail;
  if (ier>0 && n==0) n=1;
  lc = (n-k-1)*idim;
  ap_t = (PyArrayObject *)PyArray_FromDims(1,&n,PyArray_DOUBLE);
  ap_c = (PyArrayObject *)PyArray_FromDims(1,&lc,PyArray_DOUBLE);
  if (ap_t == NULL || ap_c == NULL) goto fail;
  if ((iopt==0)||(n>no)) {
    ap_wrk = (PyArrayObject *)PyArray_FromDims(1,&n,PyArray_DOUBLE);
    ap_iwrk = (PyArrayObject *)PyArray_FromDims(1,&n,PyArray_INT);
    if (ap_wrk == NULL || ap_iwrk == NULL) goto fail;
  }
  memcpy(ap_t->data,t,n*sizeof(double));
  for (i=0;i<idim;i++)
    memcpy((double *) ap_c->data+i*(n-k-1),c+i*n,(n-k-1)*sizeof(double));
  memcpy(ap_wrk->data,wrk,n*sizeof(double));
  memcpy(ap_iwrk->data,iwrk,n*sizeof(int));  
  if (wa) free(wa);
  Py_DECREF(ap_x);
  Py_DECREF(ap_w);
  return Py_BuildValue("NN{s:N,s:d,s:d,s:N,s:N,s:i,s:d}",PyArray_Return(ap_t),PyArray_Return(ap_c),"u",PyArray_Return(ap_u),"ub",ub,"ue",ue,"wrk",PyArray_Return(ap_wrk),"iwrk",PyArray_Return(ap_iwrk),"ier",ier,"fp",fp);
  fail:
  if (wa) free(wa);
  Py_XDECREF(ap_x);
  Py_XDECREF(ap_u);
  Py_XDECREF(ap_w);
  Py_XDECREF(ap_t);
  Py_XDECREF(ap_wrk);
  Py_XDECREF(ap_iwrk);
  return NULL;
}

static char doc_curfit[] = " [t,c,o] = _curfit(x,y,w,xb,xe,k,iopt,s,t,nest,wrk,iwrk,per)";
static PyObject *fitpack_curfit(PyObject *dummy, PyObject *args) {
  int iopt,m,k,nest,n,lwrk,*iwrk,ier,lwa,lc,no=0,per;
  double *x,*y,*w,xb,xe,s,*t,*c,fp,*wrk,*wa = NULL;
  PyArrayObject *ap_x = NULL,*ap_y = NULL,*ap_w = NULL,*ap_t = NULL,*ap_c = NULL;
  PyArrayObject *ap_wrk = NULL,*ap_iwrk = NULL;
  PyObject *x_py = NULL,*y_py = NULL,*w_py = NULL,*t_py = NULL;
  PyObject *wrk_py=NULL,*iwrk_py=NULL;
  if (!PyArg_ParseTuple(args, "OOOddiidOiOOi",&x_py,&y_py,&w_py,&xb,&xe,\
			&k,&iopt,&s,&t_py,&nest,&wrk_py,&iwrk_py,&per)) return NULL;
  ap_x = (PyArrayObject *)PyArray_ContiguousFromObject(x_py, PyArray_DOUBLE, 0, 1);
  ap_y = (PyArrayObject *)PyArray_ContiguousFromObject(y_py, PyArray_DOUBLE, 0, 1);
  ap_w = (PyArrayObject *)PyArray_ContiguousFromObject(w_py, PyArray_DOUBLE, 0, 1);
  ap_wrk=(PyArrayObject *)PyArray_ContiguousFromObject(wrk_py, PyArray_DOUBLE, 0, 1);
  ap_iwrk=(PyArrayObject *)PyArray_ContiguousFromObject(iwrk_py, PyArray_INT, 0, 1);
  if (ap_x == NULL || ap_y == NULL || ap_w == NULL || ap_wrk == NULL || ap_iwrk == NULL) goto fail;
  x = (double *) ap_x->data;
  y = (double *) ap_y->data;
  w = (double *) ap_w->data;
  m = ap_x->dimensions[0];
  if (per) lwrk = m*(k+1) + nest*(8+5*k);
  else lwrk = m*(k+1) + nest*(7+3*k);
  lwa = 3*nest+lwrk;
  if ((wa = (double *)malloc(lwa*sizeof(double)))==NULL) {
    PyErr_NoMemory();
    goto fail;
  }
  t = wa;
  c = t + nest;
  wrk = c + nest;
  iwrk = (int *)(wrk + lwrk);
  if (iopt) {
    ap_t=(PyArrayObject *)PyArray_ContiguousFromObject(t_py, PyArray_DOUBLE, 0, 1);
    if (ap_t == NULL) goto fail;
    n = no = ap_t->dimensions[0];
    memcpy(t,ap_t->data,n*sizeof(double));
  }
  if (iopt==1) {
    memcpy(wrk,ap_wrk->data,n*sizeof(double));
    memcpy(iwrk,ap_iwrk->data,n*sizeof(int));
  }
  if (per)
    PERCUR(&iopt,&m,x,y,w,&k,&s,&nest,&n,t,c,&fp,wrk,&lwrk,iwrk,&ier);
  else
    CURFIT(&iopt,&m,x,y,w,&xb,&xe,&k,&s,&nest,&n,t,c,&fp,wrk,&lwrk,iwrk,&ier);
  if (ier==10) {
	  PyErr_SetString(PyExc_ValueError, "Invalid inputs.");
	  goto fail;
  }
  lc = n-k-1;
  if (!iopt) {
	  ap_t = (PyArrayObject *)PyArray_FromDims(1,&n,PyArray_DOUBLE);
	  if (ap_t == NULL) goto fail;
  }
  ap_c = (PyArrayObject *)PyArray_FromDims(1,&lc,PyArray_DOUBLE);
  if (ap_c == NULL) goto fail;
  if ((iopt==0)||(n>no)) {
    Py_XDECREF(ap_wrk);
    Py_XDECREF(ap_iwrk);
    ap_wrk = (PyArrayObject *)PyArray_FromDims(1,&n,PyArray_DOUBLE);
    ap_iwrk = (PyArrayObject *)PyArray_FromDims(1,&n,PyArray_INT);
    if (ap_wrk == NULL || ap_iwrk == NULL) goto fail;
  }
  memcpy(ap_t->data,t,n*sizeof(double));
  memcpy(ap_c->data,c,lc*sizeof(double));
  memcpy(ap_wrk->data,wrk,n*sizeof(double));
  memcpy(ap_iwrk->data,iwrk,n*sizeof(int));
  if (wa) free(wa);
  Py_DECREF(ap_x);
  Py_DECREF(ap_y);
  Py_DECREF(ap_w);
  return Py_BuildValue("NN{s:N,s:N,s:i,s:d}",PyArray_Return(ap_t),PyArray_Return(ap_c),"wrk",PyArray_Return(ap_wrk),"iwrk",PyArray_Return(ap_iwrk),"ier",ier,"fp",fp);
  fail:
  if (wa) free(wa);
  Py_XDECREF(ap_x);
  Py_XDECREF(ap_y);
  Py_XDECREF(ap_w);
  Py_XDECREF(ap_t);
  Py_XDECREF(ap_wrk);
  Py_XDECREF(ap_iwrk);
  return NULL;
}

static char doc_spl_[] = " [y,ier] = _spl_(x,nu,t,c,k )";
static PyObject *fitpack_spl_(PyObject *dummy, PyObject *args) {
  int n,nu,m,ier,k;
  double *x,*y,*t,*c,*wrk = NULL;
  PyArrayObject *ap_x = NULL,*ap_y = NULL,*ap_t = NULL,*ap_c = NULL;
  PyObject *x_py = NULL,*t_py = NULL,*c_py = NULL;
  if (!PyArg_ParseTuple(args, "OiOOi",&x_py,&nu,&t_py,&c_py,&k)) return NULL;
  ap_x = (PyArrayObject *)PyArray_ContiguousFromObject(x_py, PyArray_DOUBLE, 0, 1);
  ap_t = (PyArrayObject *)PyArray_ContiguousFromObject(t_py, PyArray_DOUBLE, 0, 1);
  ap_c = (PyArrayObject *)PyArray_ContiguousFromObject(c_py, PyArray_DOUBLE, 0, 1);
  if ((ap_x == NULL || ap_t == NULL || ap_c == NULL)) goto fail;
  x = (double *) ap_x->data;
  m = ap_x->dimensions[0];
  t = (double *) ap_t->data;
  c = (double *) ap_c->data;
  n = ap_t->dimensions[0];
  ap_y = (PyArrayObject *)PyArray_FromDims(1,&m,PyArray_DOUBLE);
  if (ap_y == NULL) goto fail;
  y = (double *) ap_y->data;
  if ((wrk = (double *)malloc(n*sizeof(double)))==NULL) {
    PyErr_NoMemory();
    goto fail;
  }
  if (nu)
    SPLDER(t,&n,c,&k,&nu,x,y,&m,wrk,&ier);
  else
    SPLEV(t,&n,c,&k,x,y,&m,&ier);
  if (wrk) free(wrk);
  Py_DECREF(ap_x);
  Py_DECREF(ap_c);
  Py_DECREF(ap_t);
  return Py_BuildValue("Ni",PyArray_Return(ap_y),ier);
 fail:
  if (wrk) free(wrk);
  Py_XDECREF(ap_x);
  Py_XDECREF(ap_c);
  Py_XDECREF(ap_t);
  return NULL;
}

static char doc_splint[] = " [aint,wrk] = _splint(t,c,k,a,b)";
static PyObject *fitpack_splint(PyObject *dummy, PyObject *args) {
  int n,k;
  double *t,*c,*wrk = NULL,a,b,aint;
  PyArrayObject *ap_t = NULL,*ap_c = NULL;
  PyArrayObject *ap_wrk = NULL;
  PyObject *t_py = NULL,*c_py = NULL;
  if (!PyArg_ParseTuple(args, "OOidd",&t_py,&c_py,&k,&a,&b)) return NULL;
  ap_t = (PyArrayObject *)PyArray_ContiguousFromObject(t_py, PyArray_DOUBLE, 0, 1);
  ap_c = (PyArrayObject *)PyArray_ContiguousFromObject(c_py, PyArray_DOUBLE, 0, 1);
  if ((ap_t == NULL || ap_c == NULL)) goto fail;
  t = (double *) ap_t->data;
  c = (double *) ap_c->data;
  n = ap_t->dimensions[0];
  ap_wrk = (PyArrayObject *)PyArray_FromDims(1,&n,PyArray_DOUBLE);
  if (ap_wrk == NULL) goto fail;
  wrk = (double *) ap_wrk->data;
  aint = SPLINT(t,&n,c,&k,&a,&b,wrk);
  Py_DECREF(ap_c);
  Py_DECREF(ap_t);
  return Py_BuildValue("dN",aint,PyArray_Return(ap_wrk));
 fail:
  Py_XDECREF(ap_c);
  Py_XDECREF(ap_t);
  return NULL;
}

static char doc_sproot[] = " [z,ier] = _sproot(t,c,k,mest)";
static PyObject *fitpack_sproot(PyObject *dummy, PyObject *args) {
  int n,k,mest,ier,m;
  double *t,*c,*z=NULL;
  PyArrayObject *ap_t = NULL,*ap_c = NULL;
  PyArrayObject *ap_z = NULL;
  PyObject *t_py = NULL,*c_py = NULL;
  if (!PyArg_ParseTuple(args, "OOii",&t_py,&c_py,&k,&mest)) return NULL;
  ap_t = (PyArrayObject *)PyArray_ContiguousFromObject(t_py, PyArray_DOUBLE, 0, 1);
  ap_c = (PyArrayObject *)PyArray_ContiguousFromObject(c_py, PyArray_DOUBLE, 0, 1);
  if ((ap_t == NULL || ap_c == NULL)) goto fail;
  t = (double *) ap_t->data;
  c = (double *) ap_c->data;
  n = ap_t->dimensions[0];
  if ((z = (double *)malloc(mest*sizeof(double)))==NULL) {
    PyErr_NoMemory();
    goto fail;
  }
  SPROOT(t,&n,c,z,&mest,&m,&ier);
  if (ier==10) m=0;
  ap_z = (PyArrayObject *)PyArray_FromDims(1,&m,PyArray_DOUBLE);
  if (ap_z == NULL) goto fail;
  memcpy(ap_z->data,z,m*sizeof(double));
  if (z) free(z);
  Py_DECREF(ap_c);
  Py_DECREF(ap_t);
  return Py_BuildValue("Ni",PyArray_Return(ap_z),ier);
 fail:
  if (z) free(z);
  Py_XDECREF(ap_c);
  Py_XDECREF(ap_t);
  return NULL;
}

static char doc_spalde[] = " [d,ier] = _spalde(t,c,k,x)";
static PyObject *fitpack_spalde(PyObject *dummy, PyObject *args) {
  int n,k,k1,ier;
  double *t,*c,*d=NULL,x;
  PyArrayObject *ap_t = NULL,*ap_c = NULL,*ap_d = NULL;
  PyObject *t_py = NULL,*c_py = NULL;
  if (!PyArg_ParseTuple(args, "OOid",&t_py,&c_py,&k,&x)) return NULL;
  ap_t = (PyArrayObject *)PyArray_ContiguousFromObject(t_py, PyArray_DOUBLE, 0, 1);
  ap_c = (PyArrayObject *)PyArray_ContiguousFromObject(c_py, PyArray_DOUBLE, 0, 1);
  if ((ap_t == NULL || ap_c == NULL)) goto fail;
  t = (double *) ap_t->data;
  c = (double *) ap_c->data;
  n = ap_t->dimensions[0];
  k1=k+1;
  ap_d = (PyArrayObject *)PyArray_FromDims(1,&k1,PyArray_DOUBLE);
  if (ap_d == NULL) goto fail;
  d = (double *) ap_d->data;
  SPALDE(t,&n,c,&k1,&x,d,&ier);
  Py_DECREF(ap_c);
  Py_DECREF(ap_t);
  return Py_BuildValue("Ni",PyArray_Return(ap_d),ier);
 fail:
  Py_XDECREF(ap_c);
  Py_XDECREF(ap_t);
  return NULL;
}

static char doc_insert[] = " [tt,cc,ier] = _insert(iopt,t,c,k,x,m)";
static PyObject *fitpack_insert(PyObject *dummy, PyObject*args) {
  int iopt, n, nn, k, nest, ier, m;
  double x;
  double *t, *c, *tt, *cc;
  PyArrayObject *ap_t = NULL, *ap_c = NULL, *ap_tt = NULL, *ap_cc = NULL;
  PyObject *t_py = NULL, *c_py = NULL;
  PyObject *ret = NULL;
  if (!PyArg_ParseTuple(args, "iOOidi",&iopt,&t_py,&c_py,&k, &x, &m)) return NULL;
  ap_t = (PyArrayObject *)PyArray_ContiguousFromObject(t_py, PyArray_DOUBLE, 0, 1);
  ap_c = (PyArrayObject *)PyArray_ContiguousFromObject(c_py, PyArray_DOUBLE, 0, 1);
  if (ap_t == NULL || ap_c == NULL) goto fail;
  t = (double *) ap_t->data;
  c = (double *) ap_c->data;
  n = ap_t->dimensions[0];
  nest = n + m;
  ap_tt = (PyArrayObject *)PyArray_FromDims(1,&nest,PyArray_DOUBLE);
  ap_cc = (PyArrayObject *)PyArray_FromDims(1,&nest,PyArray_DOUBLE);
  if (ap_tt == NULL || ap_cc == NULL) goto fail;
  tt = (double *) ap_tt->data;
  cc = (double *) ap_cc->data;
  for ( ; n < nest; n++) {
    INSERT(&iopt, t, &n, c, &k, &x, tt, &nn, cc, &nest, &ier);
    if (ier) break;
    t = tt;
    c = cc;
  }
  Py_DECREF(ap_c);
  Py_DECREF(ap_t);
  ret = Py_BuildValue("NNi",PyArray_Return(ap_tt),PyArray_Return(ap_cc),ier);
  return ret;
  
  fail:
  Py_XDECREF(ap_c);
  Py_XDECREF(ap_t);
  return NULL;
  }


static void
_deBoor(double *t, double x, int k, int ell, double *result) {
    /* On completion the result array stores 
       the k+1 non-zero values of beta_i,k(x):  for i=ell, ell-1, ell-2, ell-k.
       Where t[ell] <= x < t[ell+1]. 
    */
    /* Implements the recursive algorithm of deBoor and Cox
       Equivalent to what is done in fpbspl.f in fitpack except
       their is no upper bound on the order of the spline. 
    */

    double *hh = result + k + 1;
    double *h = result;
    double xb, xa, w;
    int ind, j, m;
    
    result[0] = 1.0;
    for (j=1; j<=k; j++) {
        memcpy(hh, h, j*sizeof(double));
        h[0] = 0.0;
        for (m=1; m<=j; m++) {
            ind = ell + m;
            xb = t[ind];
            xa = t[ind-j];
            if (xb == xa) {
                h[m] = 0.0;
                continue;
            }
            w = hh[m-1]/(xb-xa);
            h[m-1] += w*(xb-x);
            h[m] = w*(x-xa);
        }
    }
}


/* Given a set of (N+1) samples:  A default set of knots is constructed
   using the samples xk plus 2*(K-1) additional knots where 
   K = max(order,1) and the knots are chosen so that distances
   are symmetric around the first and last samples: x_0 and x_N.

   There should be a vector of N+K coefficients for the spline
   curve in coef.  These coefficients form the curve as
   
   s(x) = sum(c_j B_{j,K}(x), j=-K..N-1)

   The spline function is evaluated at all points xx. 
   The approximation interval is from xk[0] to xk[-1]
   Any xx outside that interval is set automatically to 0.0
 */
static char doc_bspleval[] = "y = _bspleval(xx,xk,coef,order)";
static PyObject *_bspleval(PyObject *dummy, PyObject *args) {
    int k,kk,N,i,ell;
    PyObject *xx_py=NULL, *coef_py=NULL, *x_i_py=NULL;
    PyArrayObject *xx=NULL, *coef=NULL, *x_i=NULL, *yy=NULL;
    PyArrayIterObject *xx_iter;
    double *t=NULL, *h=NULL, *ptr;
    double x0, xN, xN1, arg, sp, cval;
    if (!PyArg_ParseTuple(args, "OOOi", &xx_py, &x_i_py, &coef_py, &k)) 
        return NULL;
    if (k < 0) {
        PyErr_Format(PyExc_ValueError, "order (%d) must be >=0", k);
        return NULL;
    }
    kk = k;
    if (k==0) kk = 1;        
    x_i = (PyArrayObject *)PyArray_FROMANY(x_i_py, NPY_DOUBLE, 1, 1, NPY_ALIGNED);
    coef = (PyArrayObject *)PyArray_FROMANY(coef_py, NPY_DOUBLE, 1, 1, NPY_ALIGNED);
    xx = (PyArrayObject *)PyArray_FROMANY(xx_py, NPY_DOUBLE, 0, 0, NPY_ALIGNED);
    if (x_i == NULL || coef == NULL || xx == NULL) goto fail;

    N = PyArray_DIM(x_i,0)-1;

    if (PyArray_DIM(coef,0) < (N+kk)) {
        PyErr_Format(PyExc_ValueError, "too few coefficients (have %d need at least %d)", 
                     PyArray_DIM(coef,0), N+kk);
        goto fail;
    }
    
    /* create output values */
    yy = (PyArrayObject *)PyArray_EMPTY(xx->nd, xx->dimensions, NPY_DOUBLE, 0);
    if (yy == NULL) goto fail;
    /* create dummy knot array with new knots inserted at the end
       selected as mirror symmetric versions of the old knots
     */
    t = (double *)malloc(sizeof(double)*(N+2*kk-1));
    if (t==NULL) {
        PyErr_NoMemory();
        goto fail;
    }
    x0 = *((double *)PyArray_DATA(x_i));
    xN = *((double *)PyArray_DATA(x_i) + N);
    for (i=0; i<kk-1; i++) { /* fill in ends if kk > 1*/
        t[i] = 2*x0 - *((double *)(PyArray_GETPTR1(x_i,kk-1-i)));
        t[kk+N+i] = 2*xN - *((double *)(PyArray_GETPTR1(x_i,N-1-i)));
    }
    ptr = t + (kk-1);
    for (i=0; i<=N; i++) {
        *ptr++ = *((double *)(PyArray_GETPTR1(x_i, i)));
    }
   
    /* Create work array to hold computed non-zero values for
       the spline for a value of x. 
    */
    h = (double *)malloc(sizeof(double)*(2*kk+1));
    if (h==NULL) {
        PyErr_NoMemory();
        goto fail;
    }

    /* Determine the spline for each value of x */ 
    xx_iter = (PyArrayIterObject *)PyArray_IterNew((PyObject *)xx);
    if (xx_iter == NULL) goto fail;
    ptr = PyArray_DATA(yy);
    while(PyArray_ITER_NOTDONE(xx_iter)) {
        arg = *((double *)PyArray_ITER_DATA(xx_iter));
        if ((arg < x0) || (arg > xN)) {
            /* If we are outside the interpolation region, 
               fill with zeros
            */
            *ptr++ = 0.0;
        }
        else {
            /* Find the interval that arg lies between in the set of knots 
               t[ell] <= arg < t[ell+1] (last-knot use the previous interval) */
            xN1 = *((double *)PyArray_DATA(x_i) + N-1);
            if (arg >= xN1) {
                ell = N + kk - 2;
            }
            else {
                ell = kk-1;
                while ((arg > t[ell])) ell++;
                ell -= 1;
            }
            
            _deBoor(t, arg, k, ell, h);
           
            sp = 0.0;
            for (i=0; i<=k; i++) {
                cval = *((double *)(PyArray_GETPTR1(coef, ell-i+1)));
                sp += cval*h[k-i];
            }
            *ptr++ = sp;
        }
        PyArray_ITER_NEXT(xx_iter);
    }
    Py_DECREF(xx_iter);
    Py_DECREF(x_i);
    Py_DECREF(coef);
    Py_DECREF(xx);
    free(t);
    free(h);
    return PyArray_Return(yy);
  
 fail:
    Py_XDECREF(xx);
    Py_XDECREF(coef);
    Py_XDECREF(x_i);
    Py_XDECREF(yy);
    if (t != NULL) free(t);
    if (h != NULL) free(h);
    return NULL;
}


/* Given a set of (N+1) sample positions:
   Construct the diagonals of the (N+1) x (N+K) matrix that is needed to find
   the coefficients of a spline fit of order K.  
       Note that K>=2 because for K=0,1, the coefficients are just the 
       sample values themselves.

   The equation that expresses the constraints is
   
     s(x_i) = sum(c_j B_{j,K}(x_i), j=-K..N-1) = w_i   for i=0..N

   This is equivalent to 

     w = B*c   where c.T = [c_{-K}, c{-K+1}, ..., c_{N-1}] and
                     w.T = [w_{0}, w_{1}, ..., w_{N}]

   Therefore B is an (N+1) times (N+K) matrix with entries

   B_{j,K}(x_i)  for column j=-K..N-1 
                 and row i=0..N

   This routine takes the N+1 sample positions and the order k and 
      constructs the banded constraint matrix B (with k+1 non-zero diagonals)

   The returned array is (N+1) times (N+K) ready to be either used
   to compute a minimally Kth-order derivative discontinuous spline
   or to be expanded with an additional K-1 constraints to be used in 
   an exact reconstruction approach.  
 */
static char doc_bsplmat[] = "B = _bsplmat(order,xk)\n"
"Construct the constraint matrix for spline fitting of order k\n"
"given sample positions in xk.\n"
"\n"
"If xk is an integer (N+1), then the result is equivalent to\n"
"xk=arange(N+1)+x0 for any value of x0.   This produces the\n"
"integer-spaced, or cardinal spline matrix a bit faster.";
static PyObject *_bsplmat(PyObject *dummy, PyObject *args) {
    int k,N,i,numbytes,j, equal;
    int dims[2];
    PyObject *x_i_py=NULL;
    PyArrayObject *x_i=NULL, *BB=NULL;
    double *t=NULL, *h=NULL, *ptr;
    double x0, xN, arg;
    if (!PyArg_ParseTuple(args, "iO", &k, &x_i_py)) 
        return NULL;
    if (k < 2) {
        PyErr_Format(PyExc_ValueError, "order (%d) must be >=2", k);
        return NULL;
    }

    equal = 0;
    N = PySequence_Length(x_i_py);
    if (N == -1 && PyErr_Occurred()) {
        PyErr_Clear();
        N = PyInt_AsLong(x_i_py);
        if (N==-1 && PyErr_Occurred()) goto fail;
        equal = 1;
    }
    N -= 1;
    
    /* create output matrix */
    dims[0] = N+1;
    dims[1] = N+k;
    BB = (PyArrayObject *)PyArray_ZEROS(2, dims, NPY_DOUBLE, 0);
    if (BB == NULL) goto fail;

    t = (double *)malloc(sizeof(double)*(N+2*k-1));
    if (t==NULL) {
        PyErr_NoMemory();
        goto fail;
    }

    /* Create work array to hold computed non-zero values for
       the spline for a value of x. 
    */
    h = (double *)malloc(sizeof(double)*(2*k+1));
    if (h==NULL) {
        PyErr_NoMemory();
        goto fail;
    }

    numbytes = k*sizeof(double);

    if (equal) { /* points equally spaced by 1 */
        /* we run deBoor's algorithm one time with artificially created knots
           Then, we keep copying the result to every row */

        /* Create knots at equally-spaced locations from -(K-1) to N+K-1 */
        ptr = t;
        for (i=-k+1; i<N+k; i++) *ptr++ = i;
        j = k-1;
        _deBoor(t, 0, k, k-1, h);
        ptr = PyArray_DATA(BB);
        N = N+1;
        for (i=0; i<N; i++) {
            memcpy(ptr, h, numbytes);
            ptr += (N+k);
        }
        goto finish;
    }

    /* Not-equally spaced */
    x_i = (PyArrayObject *)PyArray_FROMANY(x_i_py, NPY_DOUBLE, 1, 1, NPY_ALIGNED);
    if (x_i == NULL) return NULL;

    /* create dummy knot array with new knots inserted at the end
       selected as mirror symmetric versions of the old knots
     */
    x0 = *((double *)PyArray_DATA(x_i));
    xN = *((double *)PyArray_DATA(x_i) + N);
    for (i=0; i<k-1; i++) { /* fill in ends if k > 1*/
        t[i] = 2*x0 - *((double *)(PyArray_GETPTR1(x_i,k-1-i)));
        t[k+N+i] = 2*xN - *((double *)(PyArray_GETPTR1(x_i,N-1-i)));
    }
    ptr = t + (k-1);
    for (i=0; i<=N; i++) {
        *ptr++ = *((double *)(PyArray_GETPTR1(x_i, i)));
    }
   

    /* Determine the K+1 non-zero values of the spline and place them in the 
       correct location in the matrix for each row (along the diagonals). 
       In fact, the last member is always zero so only K non-zero values
       are present. 
    */
    ptr = PyArray_DATA(BB);
    for (i=0,j=k-1; i<N; i++,j++) {
        arg = *((double *)PyArray_DATA(x_i) + i);
        _deBoor(t, arg, k, j, h);
        memcpy(ptr, h, numbytes);
        ptr += (N+k+1);  /* advance to next row shifted over one */
    }
    /* Last row is different the first coefficient is zero.*/
    _deBoor(t, xN, k, j-1, h);
    memcpy(ptr, h+1, numbytes);

 finish:
    Py_XDECREF(x_i);
    free(t);
    free(h);
    return (PyObject *)BB;
  
 fail:
    Py_XDECREF(x_i);
    Py_XDECREF(BB);
    if (t != NULL) free(t);
    if (h != NULL) free(h);
    return NULL;
}



