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
    ap_wrk = (PyArrayObject *)PyArray_FromDims(1,&lc,PyArray_DOUBLE);
    if (ap_wrk == NULL) goto fail;
    /*ap_iwrk = (PyArrayObject *)PyArray_FromDims(1,&n,PyArray_INT);*/
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




