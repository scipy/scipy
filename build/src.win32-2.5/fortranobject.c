#define FORTRANOBJECT_C
#include "fortranobject.h"

#ifdef __cplusplus
extern "C" {
#endif

/*
  This file implements: FortranObject, array_from_pyobj, copy_ND_array

  Author: Pearu Peterson <pearu@cens.ioc.ee>
  $Revision: 1.52 $
  $Date: 2005/07/11 07:44:20 $
*/

int
F2PyDict_SetItemString(PyObject *dict, char *name, PyObject *obj)
{
    if (obj==NULL) {
        fprintf(stderr, "Error loading %s\n", name);
        if (PyErr_Occurred()) {
            PyErr_Print();
            PyErr_Clear();
        }
        return -1;
    }
    return PyDict_SetItemString(dict, name, obj);
}

/************************* FortranObject *******************************/

typedef PyObject *(*fortranfunc)(PyObject *,PyObject *,PyObject *,void *);

PyObject *
PyFortranObject_New(FortranDataDef* defs, f2py_void_func init) {
    int i;
    PyFortranObject *fp = NULL;
    PyObject *v = NULL;
    if (init!=NULL)                           /* Initialize F90 module objects */
        (*(init))();
    if ((fp = PyObject_New(PyFortranObject, &PyFortran_Type))==NULL) return NULL;
    if ((fp->dict = PyDict_New())==NULL) return NULL;
    fp->len = 0;
    while (defs[fp->len].name != NULL) fp->len++;
    if (fp->len == 0) goto fail;
    fp->defs = defs;
    for (i=0;i<fp->len;i++)
        if (fp->defs[i].rank == -1) {                      /* Is Fortran routine */
            v = PyFortranObject_NewAsAttr(&(fp->defs[i]));
            if (v==NULL) return NULL;
            PyDict_SetItemString(fp->dict,fp->defs[i].name,v);
        } else
            if ((fp->defs[i].data)!=NULL) { /* Is Fortran variable or array (not allocatable) */
                if (fp->defs[i].type == PyArray_STRING) {
                    int n = fp->defs[i].rank-1;
                    v = PyArray_New(&PyArray_Type, n, fp->defs[i].dims.d,
                                    PyArray_STRING, NULL, fp->defs[i].data, fp->defs[i].dims.d[n],
                                    NPY_FARRAY, NULL);
                }
                else {
                    v = PyArray_New(&PyArray_Type, fp->defs[i].rank, fp->defs[i].dims.d,
                                    fp->defs[i].type, NULL, fp->defs[i].data, 0, NPY_FARRAY,
                                    NULL);
                }
                if (v==NULL) return NULL;
                PyDict_SetItemString(fp->dict,fp->defs[i].name,v);
            }
    Py_XDECREF(v);
    return (PyObject *)fp;
 fail:
    Py_XDECREF(v);
    return NULL;
}

PyObject *
PyFortranObject_NewAsAttr(FortranDataDef* defs) { /* used for calling F90 module routines */
    PyFortranObject *fp = NULL;
    fp = PyObject_New(PyFortranObject, &PyFortran_Type);
    if (fp == NULL) return NULL;
    if ((fp->dict = PyDict_New())==NULL) return NULL;
    fp->len = 1;
    fp->defs = defs;
    return (PyObject *)fp;
}

/* Fortran methods */

static void
fortran_dealloc(PyFortranObject *fp) {
    Py_XDECREF(fp->dict);
    PyMem_Del(fp);
}


static PyMethodDef fortran_methods[] = {
    {NULL,          NULL}           /* sentinel */
};


static PyObject *
fortran_doc (FortranDataDef def) {
    char *p;
    PyObject *s = NULL;
    int i;
    unsigned size=100;
    if (def.doc!=NULL)
        size += strlen(def.doc);
    p = (char*)malloc (size);
    if (sprintf(p,"%s - ",def.name)==0) goto fail;
    if (def.rank==-1) {
        if (def.doc==NULL) {
            if (sprintf(p,"%sno docs available",p)==0)
                goto fail;
        } else {
            if (sprintf(p,"%s%s",p,def.doc)==0)
                goto fail;
        }
    } else {
        PyArray_Descr *d = PyArray_DescrFromType(def.type);
        if (sprintf(p,"%s'%c'-",p,d->type)==0) {Py_DECREF(d); goto fail;}
        Py_DECREF(d);
        if (def.data==NULL) {
            if (sprintf(p,"%sarray(%" NPY_INTP_FMT,p,def.dims.d[0])==0) goto fail;
            for(i=1;i<def.rank;++i)
                if (sprintf(p,"%s,%" NPY_INTP_FMT,p,def.dims.d[i])==0) goto fail;
            if (sprintf(p,"%s), not allocated",p)==0) goto fail;
        } else {
            if (def.rank>0) {
                if (sprintf(p,"%sarray(%"NPY_INTP_FMT,p,def.dims.d[0])==0) goto fail;
                for(i=1;i<def.rank;i++)
                    if (sprintf(p,"%s,%" NPY_INTP_FMT,p,def.dims.d[i])==0) goto fail;
                if (sprintf(p,"%s)",p)==0) goto fail;
            } else {
                if (sprintf(p,"%sscalar",p)==0) goto fail;
            }
        }
    }
    if (sprintf(p,"%s\n",p)==0) goto fail;
    if (strlen(p)>size) {
        fprintf(stderr,"fortranobject.c:fortran_doc:len(p)=%zd>%d(size): too long doc string required, increase size\n",strlen(p),size);
        goto fail;
    }
    s = PyString_FromString(p);
 fail:
    free(p);
    return s;
}

static FortranDataDef *save_def; /* save pointer of an allocatable array */
static void set_data(char *d,npy_intp *f) {  /* callback from Fortran */
    if (*f)                               /* In fortran f=allocated(d) */
        save_def->data = d;
    else
        save_def->data = NULL;
    /* printf("set_data: d=%p,f=%d\n",d,*f); */
}

static PyObject *
fortran_getattr(PyFortranObject *fp, char *name) {
    int i,j,k,flag;
    if (fp->dict != NULL) {
        PyObject *v = PyDict_GetItemString(fp->dict, name);
        if (v != NULL) {
            Py_INCREF(v);
            return v;
        }
    }
    for (i=0,j=1;i<fp->len && (j=strcmp(name,fp->defs[i].name));i++);
    if (j==0)
        if (fp->defs[i].rank!=-1) {                   /* F90 allocatable array */
            if (fp->defs[i].func==NULL) return NULL;
            for(k=0;k<fp->defs[i].rank;++k)
                fp->defs[i].dims.d[k]=-1;
            save_def = &fp->defs[i];
            (*(fp->defs[i].func))(&fp->defs[i].rank,fp->defs[i].dims.d,set_data,&flag);
            if (flag==2)
                k = fp->defs[i].rank + 1;
            else
                k = fp->defs[i].rank;
            if (fp->defs[i].data !=NULL) {              /* array is allocated */
                PyObject *v = PyArray_New(&PyArray_Type, k, fp->defs[i].dims.d,
                                          fp->defs[i].type, NULL, fp->defs[i].data, 0, NPY_FARRAY,
                                          NULL);
                if (v==NULL) return NULL;
                /* Py_INCREF(v); */
                return v;
            } else {                                    /* array is not allocated */
                Py_INCREF(Py_None);
                return Py_None;
            }
        }
    if (strcmp(name,"__dict__")==0) {
        Py_INCREF(fp->dict);
        return fp->dict;
    }
    if (strcmp(name,"__doc__")==0) {
        PyObject *s = PyString_FromString("");
        for (i=0;i<fp->len;i++)
            PyString_ConcatAndDel(&s,fortran_doc(fp->defs[i]));
        if (PyDict_SetItemString(fp->dict, name, s))
            return NULL;
        return s;
    }
    if ((strcmp(name,"_cpointer")==0) && (fp->len==1)) {
        PyObject *cobj = PyCObject_FromVoidPtr((void *)(fp->defs[0].data),NULL);
        if (PyDict_SetItemString(fp->dict, name, cobj))
            return NULL;
        return cobj;
    }
    return Py_FindMethod(fortran_methods, (PyObject *)fp, name);
}

static int
fortran_setattr(PyFortranObject *fp, char *name, PyObject *v) {
    int i,j,flag;
    PyArrayObject *arr = NULL;
    for (i=0,j=1;i<fp->len && (j=strcmp(name,fp->defs[i].name));i++);
    if (j==0) {
        if (fp->defs[i].rank==-1) {
            PyErr_SetString(PyExc_AttributeError,"over-writing fortran routine");
            return -1;
        }
        if (fp->defs[i].func!=NULL) { /* is allocatable array */
            npy_intp dims[F2PY_MAX_DIMS];
            int k;
            save_def = &fp->defs[i];
            if (v!=Py_None) {     /* set new value (reallocate if needed --
                                     see f2py generated code for more
                                     details ) */
                for(k=0;k<fp->defs[i].rank;k++) dims[k]=-1;
                if ((arr = array_from_pyobj(fp->defs[i].type,dims,fp->defs[i].rank,F2PY_INTENT_IN,v))==NULL)
                    return -1;
                (*(fp->defs[i].func))(&fp->defs[i].rank,arr->dimensions,set_data,&flag);
            } else {             /* deallocate */
                for(k=0;k<fp->defs[i].rank;k++) dims[k]=0;
                (*(fp->defs[i].func))(&fp->defs[i].rank,dims,set_data,&flag);
                for(k=0;k<fp->defs[i].rank;k++) dims[k]=-1;
            }
            memcpy(fp->defs[i].dims.d,dims,fp->defs[i].rank*sizeof(npy_intp));
        } else {                     /* not allocatable array */
            if ((arr = array_from_pyobj(fp->defs[i].type,fp->defs[i].dims.d,fp->defs[i].rank,F2PY_INTENT_IN,v))==NULL)
                return -1;
        }
        if (fp->defs[i].data!=NULL) { /* copy Python object to Fortran array */
            npy_intp s = PyArray_MultiplyList(fp->defs[i].dims.d,arr->nd);
            if (s==-1)
                s = PyArray_MultiplyList(arr->dimensions,arr->nd);
            if (s<0 ||
                (memcpy(fp->defs[i].data,arr->data,s*PyArray_ITEMSIZE(arr)))==NULL) {
                if ((PyObject*)arr!=v) {
                    Py_DECREF(arr);
                }
                return -1;
            }
            if ((PyObject*)arr!=v) {
                Py_DECREF(arr);
            }
        } else return (fp->defs[i].func==NULL?-1:0);
        return 0; /* succesful */
    }
    if (fp->dict == NULL) {
        fp->dict = PyDict_New();
        if (fp->dict == NULL)
            return -1;
    }
    if (v == NULL) {
        int rv = PyDict_DelItemString(fp->dict, name);
        if (rv < 0)
            PyErr_SetString(PyExc_AttributeError,"delete non-existing fortran attribute");
        return rv;
    }
    else
        return PyDict_SetItemString(fp->dict, name, v);
}

static PyObject*
fortran_call(PyFortranObject *fp, PyObject *arg, PyObject *kw) {
    int i = 0;
    /*  printf("fortran call
        name=%s,func=%p,data=%p,%p\n",fp->defs[i].name,
        fp->defs[i].func,fp->defs[i].data,&fp->defs[i].data); */
    if (fp->defs[i].rank==-1) {/* is Fortran routine */
        if ((fp->defs[i].func==NULL)) {
            PyErr_Format(PyExc_RuntimeError, "no function to call");
            return NULL;
        }
        else if (fp->defs[i].data==NULL)
            /* dummy routine */
            return (*((fortranfunc)(fp->defs[i].func)))((PyObject *)fp,arg,kw,NULL);
        else
            return (*((fortranfunc)(fp->defs[i].func)))((PyObject *)fp,arg,kw,
                                                        (void *)fp->defs[i].data);
    }
    PyErr_Format(PyExc_TypeError, "this fortran object is not callable");
    return NULL;
}


PyTypeObject PyFortran_Type = {
    PyObject_HEAD_INIT(0)
    0,                    /*ob_size*/
    "fortran",                    /*tp_name*/
    sizeof(PyFortranObject),      /*tp_basicsize*/
    0,                    /*tp_itemsize*/
    /* methods */
    (destructor)fortran_dealloc, /*tp_dealloc*/
    0,                    /*tp_print*/
    (getattrfunc)fortran_getattr, /*tp_getattr*/
    (setattrfunc)fortran_setattr, /*tp_setattr*/
    0,                    /*tp_compare*/
    0,                    /*tp_repr*/
    0,                    /*tp_as_number*/
    0,                    /*tp_as_sequence*/
    0,                    /*tp_as_mapping*/
    0,                    /*tp_hash*/
    (ternaryfunc)fortran_call,                    /*tp_call*/
};

/************************* f2py_report_atexit *******************************/

#ifdef F2PY_REPORT_ATEXIT
static int passed_time = 0;
static int passed_counter = 0;
static int passed_call_time = 0;
static struct timeb start_time;
static struct timeb stop_time;
static struct timeb start_call_time;
static struct timeb stop_call_time;
static int cb_passed_time = 0;
static int cb_passed_counter = 0;
static int cb_passed_call_time = 0;
static struct timeb cb_start_time;
static struct timeb cb_stop_time;
static struct timeb cb_start_call_time;
static struct timeb cb_stop_call_time;

extern void f2py_start_clock(void) { ftime(&start_time); }
extern
void f2py_start_call_clock(void) {
    f2py_stop_clock();
    ftime(&start_call_time);
}
extern
void f2py_stop_clock(void) {
    ftime(&stop_time);
    passed_time += 1000*(stop_time.time - start_time.time);
    passed_time += stop_time.millitm - start_time.millitm;
}
extern
void f2py_stop_call_clock(void) {
    ftime(&stop_call_time);
    passed_call_time += 1000*(stop_call_time.time - start_call_time.time);
    passed_call_time += stop_call_time.millitm - start_call_time.millitm;
    passed_counter += 1;
    f2py_start_clock();
}

extern void f2py_cb_start_clock(void) { ftime(&cb_start_time); }
extern
void f2py_cb_start_call_clock(void) {
    f2py_cb_stop_clock();
    ftime(&cb_start_call_time);
}
extern
void f2py_cb_stop_clock(void) {
    ftime(&cb_stop_time);
    cb_passed_time += 1000*(cb_stop_time.time - cb_start_time.time);
    cb_passed_time += cb_stop_time.millitm - cb_start_time.millitm;
}
extern
void f2py_cb_stop_call_clock(void) {
    ftime(&cb_stop_call_time);
    cb_passed_call_time += 1000*(cb_stop_call_time.time - cb_start_call_time.time);
    cb_passed_call_time += cb_stop_call_time.millitm - cb_start_call_time.millitm;
    cb_passed_counter += 1;
    f2py_cb_start_clock();
}

static int f2py_report_on_exit_been_here = 0;
extern
void f2py_report_on_exit(int exit_flag,void *name) {
    if (f2py_report_on_exit_been_here) {
        fprintf(stderr,"             %s\n",(char*)name);
        return;
    }
    f2py_report_on_exit_been_here = 1;
    fprintf(stderr,"                      /-----------------------\\\n");
    fprintf(stderr,"                     < F2PY performance report >\n");
    fprintf(stderr,"                      \\-----------------------/\n");
    fprintf(stderr,"Overall time spent in ...\n");
    fprintf(stderr,"(a) wrapped (Fortran/C) functions           : %8d msec\n",
            passed_call_time);
    fprintf(stderr,"(b) f2py interface,           %6d calls  : %8d msec\n",
            passed_counter,passed_time);
    fprintf(stderr,"(c) call-back (Python) functions            : %8d msec\n",
            cb_passed_call_time);
    fprintf(stderr,"(d) f2py call-back interface, %6d calls  : %8d msec\n",
            cb_passed_counter,cb_passed_time);

    fprintf(stderr,"(e) wrapped (Fortran/C) functions (acctual) : %8d msec\n\n",
            passed_call_time-cb_passed_call_time-cb_passed_time);
    fprintf(stderr,"Use -DF2PY_REPORT_ATEXIT_DISABLE to disable this message.\n");
    fprintf(stderr,"Exit status: %d\n",exit_flag);
    fprintf(stderr,"Modules    : %s\n",(char*)name);
}
#endif

/********************** report on array copy ****************************/

#ifdef F2PY_REPORT_ON_ARRAY_COPY
static void f2py_report_on_array_copy(PyArrayObject* arr) {
    const long arr_size = PyArray_Size((PyObject *)arr);
    if (arr_size>F2PY_REPORT_ON_ARRAY_COPY) {
        fprintf(stderr,"copied an array: size=%ld, elsize=%d\n",
                arr_size, PyArray_ITEMSIZE(arr));
    }
}
static void f2py_report_on_array_copy_fromany(void) {
    fprintf(stderr,"created an array from object\n");
}

#define F2PY_REPORT_ON_ARRAY_COPY_FROMARR f2py_report_on_array_copy((PyArrayObject *)arr)
#define F2PY_REPORT_ON_ARRAY_COPY_FROMANY f2py_report_on_array_copy_fromany()
#else
#define F2PY_REPORT_ON_ARRAY_COPY_FROMARR
#define F2PY_REPORT_ON_ARRAY_COPY_FROMANY
#endif


/************************* array_from_obj *******************************/

/*
 * File: array_from_pyobj.c
 *
 * Description:
 * ------------
 * Provides array_from_pyobj function that returns a contigious array
 * object with the given dimensions and required storage order, either
 * in row-major (C) or column-major (Fortran) order. The function
 * array_from_pyobj is very flexible about its Python object argument
 * that can be any number, list, tuple, or array.
 *
 * array_from_pyobj is used in f2py generated Python extension
 * modules.
 *
 * Author: Pearu Peterson <pearu@cens.ioc.ee>
 * Created: 13-16 January 2002
 * $Id: fortranobject.c,v 1.52 2005/07/11 07:44:20 pearu Exp $
 */

static int
count_nonpos(const int rank,
             const npy_intp *dims) {
    int i=0,r=0;
    while (i<rank) {
        if (dims[i] <= 0) ++r;
        ++i;
    }
    return r;
}

static int check_and_fix_dimensions(const PyArrayObject* arr,
                                    const int rank,
                                    npy_intp *dims);

#ifdef DEBUG_COPY_ND_ARRAY
void dump_dims(int rank, npy_intp* dims) {
    int i;
    printf("[");
    for(i=0;i<rank;++i) {
        printf("%3" NPY_INTP_FMT, dims[i]);
    }
    printf("]\n");
}
void dump_attrs(const PyArrayObject* arr) {
    int rank = arr->nd;
    npy_intp size = PyArray_Size((PyObject *)arr);
    printf("\trank = %d, flags = %d, size = %" NPY_INTP_FMT  "\n",
           rank,arr->flags,size);
    printf("\tstrides = ");
    dump_dims(rank,arr->strides);
    printf("\tdimensions = ");
    dump_dims(rank,arr->dimensions);
}
#endif

#define SWAPTYPE(a,b,t) {t c; c = (a); (a) = (b); (b) = c; }

static int swap_arrays(PyArrayObject* arr1, PyArrayObject* arr2) {
    SWAPTYPE(arr1->data,arr2->data,char*);
    SWAPTYPE(arr1->nd,arr2->nd,int);
    SWAPTYPE(arr1->dimensions,arr2->dimensions,npy_intp*);
    SWAPTYPE(arr1->strides,arr2->strides,npy_intp*);
    SWAPTYPE(arr1->base,arr2->base,PyObject*);
    SWAPTYPE(arr1->descr,arr2->descr,PyArray_Descr*);
    SWAPTYPE(arr1->flags,arr2->flags,int);
    /* SWAPTYPE(arr1->weakreflist,arr2->weakreflist,PyObject*); */
    return 0;
}

#define ARRAY_ISCOMPATIBLE(arr,type_num)                                \
    (  (PyArray_ISINTEGER(arr) && PyTypeNum_ISINTEGER(type_num))        \
       ||(PyArray_ISFLOAT(arr) && PyTypeNum_ISFLOAT(type_num))          \
       ||(PyArray_ISCOMPLEX(arr) && PyTypeNum_ISCOMPLEX(type_num))      \
       )

extern
PyArrayObject* array_from_pyobj(const int type_num,
                                npy_intp *dims,
                                const int rank,
                                const int intent,
                                PyObject *obj) {
    /* Note about reference counting
       -----------------------------
       If the caller returns the array to Python, it must be done with
       Py_BuildValue("N",arr).
       Otherwise, if obj!=arr then the caller must call Py_DECREF(arr).

       Note on intent(cache,out,..)
       ---------------------
       Don't expect correct data when returning intent(cache) array.

    */
    char mess[200];
    PyArrayObject *arr = NULL;
    PyArray_Descr *descr;
    char typechar;
    int elsize;

    if ((intent & F2PY_INTENT_HIDE)
        || ((intent & F2PY_INTENT_CACHE) && (obj==Py_None))
        || ((intent & F2PY_OPTIONAL) && (obj==Py_None))
        ) {
        /* intent(cache), optional, intent(hide) */
        if (count_nonpos(rank,dims)) {
            int i;
            sprintf(mess,"failed to create intent(cache|hide)|optional array"
                    "-- must have defined dimensions but got (");
            for(i=0;i<rank;++i)
                sprintf(mess+strlen(mess),"%" NPY_INTP_FMT ",",dims[i]);
            sprintf(mess+strlen(mess),")");
            PyErr_SetString(PyExc_ValueError,mess);
            return NULL;
        }
        arr = (PyArrayObject *)
            PyArray_New(&PyArray_Type, rank, dims, type_num,
                        NULL,NULL,0,
                        !(intent&F2PY_INTENT_C),
                        NULL);
        if (arr==NULL) return NULL;
        if (!(intent & F2PY_INTENT_CACHE))
            PyArray_FILLWBYTE(arr, 0);
        return arr;
    }

    descr = PyArray_DescrFromType(type_num);
    elsize = descr->elsize;
    typechar = descr->type;
    Py_DECREF(descr);
    if (PyArray_Check(obj)) {
        arr = (PyArrayObject *)obj;

        if (intent & F2PY_INTENT_CACHE) {
            /* intent(cache) */
            if (PyArray_ISONESEGMENT(obj)
                && PyArray_ITEMSIZE((PyArrayObject *)obj)>=elsize) {
                if (check_and_fix_dimensions((PyArrayObject *)obj,rank,dims)) {
                    return NULL; /*XXX: set exception */
                }
                if (intent & F2PY_INTENT_OUT)
                    Py_INCREF(obj);
                return (PyArrayObject *)obj;
            }
            sprintf(mess,"failed to initialize intent(cache) array");
            if (!PyArray_ISONESEGMENT(obj))
                sprintf(mess+strlen(mess)," -- input must be in one segment");
            if (PyArray_ITEMSIZE(arr)<elsize)
                sprintf(mess+strlen(mess)," -- expected at least elsize=%d but got %d",
                        elsize,PyArray_ITEMSIZE(arr)
                        );
            PyErr_SetString(PyExc_ValueError,mess);
            return NULL;
        }

        /* here we have always intent(in) or intent(inout) or intent(inplace) */

        if (check_and_fix_dimensions(arr,rank,dims)) {
            return NULL; /*XXX: set exception */
        }

        if ((! (intent & F2PY_INTENT_COPY))
            && PyArray_ITEMSIZE(arr)==elsize
            && ARRAY_ISCOMPATIBLE(arr,type_num)
            ) {
            if ((intent & F2PY_INTENT_C)?PyArray_ISCARRAY(arr):PyArray_ISFARRAY(arr)) {
                if ((intent & F2PY_INTENT_OUT)) {
                    Py_INCREF(arr);
                }
                /* Returning input array */
                return arr;
            }
        }

        if (intent & F2PY_INTENT_INOUT) {
            sprintf(mess,"failed to initialize intent(inout) array");
            if ((intent & F2PY_INTENT_C) && !PyArray_ISCARRAY(arr))
                sprintf(mess+strlen(mess)," -- input not contiguous");
            if (!(intent & F2PY_INTENT_C) && !PyArray_ISFARRAY(arr))
                sprintf(mess+strlen(mess)," -- input not fortran contiguous");
            if (PyArray_ITEMSIZE(arr)!=elsize)
                sprintf(mess+strlen(mess)," -- expected elsize=%d but got %d",
                        elsize,
                        PyArray_ITEMSIZE(arr)
                        );
            if (!(ARRAY_ISCOMPATIBLE(arr,type_num)))
                sprintf(mess+strlen(mess)," -- input '%c' not compatible to '%c'",
                        arr->descr->type,typechar);
            PyErr_SetString(PyExc_ValueError,mess);
            return NULL;
        }

        /* here we have always intent(in) or intent(inplace) */

        {
            PyArrayObject *retarr = (PyArrayObject *) \
                PyArray_New(&PyArray_Type, arr->nd, arr->dimensions, type_num,
                            NULL,NULL,0,
                            !(intent&F2PY_INTENT_C),
                            NULL);
            if (retarr==NULL)
                return NULL;
            F2PY_REPORT_ON_ARRAY_COPY_FROMARR;
            if (PyArray_CopyInto(retarr, arr)) {
                Py_DECREF(retarr);
                return NULL;
            }
            if (intent & F2PY_INTENT_INPLACE) {
                if (swap_arrays(arr,retarr))
                    return NULL; /* XXX: set exception */
                Py_XDECREF(retarr);
                if (intent & F2PY_INTENT_OUT)
                    Py_INCREF(arr);
            } else {
                arr = retarr;
            }
        }
        return arr;
    }

    if ((intent & F2PY_INTENT_INOUT)
        || (intent & F2PY_INTENT_INPLACE)
        || (intent & F2PY_INTENT_CACHE)) {
        sprintf(mess,"failed to initialize intent(inout|inplace|cache) array"
                " -- input must be array but got %s",
                PyString_AsString(PyObject_Str(PyObject_Type(obj)))
                );
        PyErr_SetString(PyExc_TypeError,mess);
        return NULL;
    }

    {
        F2PY_REPORT_ON_ARRAY_COPY_FROMANY;
        arr = (PyArrayObject *) \
            PyArray_FromAny(obj,PyArray_DescrFromType(type_num), 0,0,
                            ((intent & F2PY_INTENT_C)?NPY_CARRAY:NPY_FARRAY) \
                            | NPY_FORCECAST, NULL);
        if (arr==NULL)
            return NULL;
        if (check_and_fix_dimensions(arr,rank,dims))
            return NULL; /*XXX: set exception */
        return arr;
    }

}

/*****************************************/
/* Helper functions for array_from_pyobj */
/*****************************************/

static
int check_and_fix_dimensions(const PyArrayObject* arr,const int rank,npy_intp *dims) {
    /*
      This function fills in blanks (that are -1\'s) in dims list using
      the dimensions from arr. It also checks that non-blank dims will
      match with the corresponding values in arr dimensions.
    */
    const npy_intp arr_size = (arr->nd)?PyArray_Size((PyObject *)arr):1;
#ifdef DEBUG_COPY_ND_ARRAY
    dump_attrs(arr);
    printf("check_and_fix_dimensions:init: dims=");
    dump_dims(rank,dims);
#endif
    if (rank > arr->nd) { /* [1,2] -> [[1],[2]]; 1 -> [[1]]  */
        npy_intp new_size = 1;
        int free_axe = -1;
        int i;
        /* Fill dims where -1 or 0; check dimensions; calc new_size; */
        for(i=0;i<arr->nd;++i) {
            if (dims[i] >= 0) {
                if (dims[i]!=arr->dimensions[i]) {
                    fprintf(stderr,"%d-th dimension must be fixed to %" NPY_INTP_FMT
                            " but got %" NPY_INTP_FMT "\n",
                            i,dims[i], arr->dimensions[i]);
                    return 1;
                }
                if (!dims[i]) dims[i] = 1;
            } else {
                dims[i] = arr->dimensions[i] ? arr->dimensions[i] : 1;
            }
            new_size *= dims[i];
        }
        for(i=arr->nd;i<rank;++i)
            if (dims[i]>1) {
                fprintf(stderr,"%d-th dimension must be %" NPY_INTP_FMT
                        " but got 0 (not defined).\n",
                        i,dims[i]);
                return 1;
            } else if (free_axe<0)
                free_axe = i;
            else
                dims[i] = 1;
        if (free_axe>=0) {
            dims[free_axe] = arr_size/new_size;
            new_size *= dims[free_axe];
        }
        if (new_size != arr_size) {
            fprintf(stderr,"confused: new_size=%" NPY_INTP_FMT
                    ", arr_size=%" NPY_INTP_FMT " (maybe too many free"
                    " indices)\n", new_size,arr_size);
            return 1;
        }
    } else if (rank==arr->nd) {
        int i;
        npy_intp d;
        for (i=0; i<rank; ++i) {
            d = arr->dimensions[i];
            if (dims[i]>=0) {
                if (d > 1 && d!=dims[i]) {
                    fprintf(stderr,"%d-th dimension must be fixed to %" NPY_INTP_FMT
                            " but got %" NPY_INTP_FMT "\n",
                            i,dims[i],d);
                    return 1;
                }
                if (!dims[i]) dims[i] = 1;
            } else dims[i] = d;
        }
    } else { /* [[1,2]] -> [[1],[2]] */
        int i,j;
        npy_intp d;
        int effrank;
        npy_intp size;
        for (i=0,effrank=0;i<arr->nd;++i)
            if (arr->dimensions[i]>1) ++effrank;
        if (dims[rank-1]>=0)
            if (effrank>rank) {
                fprintf(stderr,"too many axes: %d (effrank=%d), expected rank=%d\n",
                        arr->nd,effrank,rank);
                return 1;
            }

        for (i=0,j=0;i<rank;++i) {
            while (j<arr->nd && arr->dimensions[j]<2) ++j;
            if (j>=arr->nd) d = 1;
            else d = arr->dimensions[j++];
            if (dims[i]>=0) {
                if (d>1 && d!=dims[i]) {
                    fprintf(stderr,"%d-th dimension must be fixed to %" NPY_INTP_FMT
                            " but got %" NPY_INTP_FMT " (real index=%d)\n",
                            i,dims[i],d,j-1);
                    return 1;
                }
                if (!dims[i]) dims[i] = 1;
            } else
                dims[i] = d;
        }

        for (i=rank;i<arr->nd;++i) { /* [[1,2],[3,4]] -> [1,2,3,4] */
            while (j<arr->nd && arr->dimensions[j]<2) ++j;
            if (j>=arr->nd) d = 1;
            else d = arr->dimensions[j++];
            dims[rank-1] *= d;
        }
        for (i=0,size=1;i<rank;++i) size *= dims[i];
        if (size != arr_size) {
            fprintf(stderr,"confused: size=%" NPY_INTP_FMT ", arr_size=%" NPY_INTP_FMT
                    ", rank=%d, effrank=%d, arr.nd=%d, dims=[",
                    size,arr_size,rank,effrank,arr->nd);
            for (i=0;i<rank;++i) fprintf(stderr," %" NPY_INTP_FMT,dims[i]);
            fprintf(stderr," ], arr.dims=[");
            for (i=0;i<arr->nd;++i) fprintf(stderr," %" NPY_INTP_FMT,arr->dimensions[i]);
            fprintf(stderr," ]\n");
            return 1;
        }
    }
#ifdef DEBUG_COPY_ND_ARRAY
    printf("check_and_fix_dimensions:end: dims=");
    dump_dims(rank,dims);
#endif
    return 0;
}

/* End of file: array_from_pyobj.c */

/************************* copy_ND_array *******************************/

extern
int copy_ND_array(const PyArrayObject *arr, PyArrayObject *out)
{
    F2PY_REPORT_ON_ARRAY_COPY_FROMARR;
    return PyArray_CopyInto(out, (PyArrayObject *)arr);
}

#ifdef __cplusplus
}
#endif
/************************* EOF fortranobject.c *******************************/
