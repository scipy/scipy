/* -*- C -*-  (not really, but good for syntax highlighting) */
%module sparsetools

 /* why does SWIG complain about int arrays? a typecheck is provided */
#pragma SWIG nowarn=467

%{
#define SWIG_FILE_WITH_INIT
#include "Python.h"
#include "numpy/arrayobject.h"
#include "complex_ops.h"
#include "sparsetools.h"
%}

%feature("autodoc", "1");

%include "numpy.i"

%init %{
    import_array();
%}


%{
/*!
  Appends @a what to @a where. On input, @a where need not to be a tuple, but on
  return it always is.

  @par Revision history:
  - 17.02.2005, c
*/
PyObject *helper_appendToTuple( PyObject *where, PyObject *what ) {
  PyObject *o2, *o3;

  if ((!where) || (where == Py_None)) {
    where = what;
  } else {
    if (!PyTuple_Check( where )) {
      o2 = where;
      where = PyTuple_New( 1 );
      PyTuple_SetItem( where, 0, o2 );
    }
    o3 = PyTuple_New( 1 );
    PyTuple_SetItem( o3, 0, what );
    o2 = where;
    where = PySequence_Concat( o2, o3 );
    Py_DECREF( o2 );
    Py_DECREF( o3 );
  }
  return where;
}




%}



/*
 * Use STL vectors for ARGOUTs
 */
%define VEC_ARRAY_ARGOUT( ctype, atype  ) 
%typemap( in, numinputs=0 ) std::vector<ctype>* array_argout( std::vector<ctype>* tmp ) {
  tmp = new std::vector<ctype>(); 
  $1 = tmp; 
}; 
%typemap( argout ) std::vector<ctype>* array_argout { 
  int length = ($1)->size(); 
  PyObject *obj = PyArray_FromDims(1, &length, PyArray_##atype); 
  memcpy(PyArray_DATA(obj),&((*($1))[0]),sizeof(ctype)*length);	 
  delete $1; 
  $result = helper_appendToTuple( $result, (PyObject *)obj ); 
}; 
%enddef




%include "typemaps.i"

%typemap(typecheck) int *, const int *, int [], const int []
{
  $1 = (is_array($input) && PyArray_CanCastSafely(PyArray_TYPE($input),PyArray_INT)) ? 1 : 0;
};

%typemap(typecheck) float *, const float *, float [], const float []
{
  $1 = (is_array($input) && PyArray_CanCastSafely(PyArray_TYPE($input),PyArray_FLOAT)) ? 1 : 0;
};

%typemap(typecheck) double *, const double *, double [], const double []
{
  $1 = (is_array($input) && PyArray_CanCastSafely(PyArray_TYPE($input),PyArray_DOUBLE)) ? 1 : 0;
};

%typemap(typecheck) long double *, const long double *, long double [], const long double []
{
  $1 = (is_array($input) && PyArray_CanCastSafely(PyArray_TYPE($input),PyArray_LONGDOUBLE)) ? 1 : 0;
};

%typemap(typecheck) npy_cfloat *, const npy_cfloat *, npy_cfloat [], const npy_cfloat []
{
  $1 = (is_array($input) && PyArray_CanCastSafely(PyArray_TYPE($input),PyArray_CFLOAT)) ? 1 : 0;
};

%typemap(typecheck) npy_cdouble *, const npy_cdouble *, npy_cdouble [], const npy_cdouble []
{
  $1 = (is_array($input) && PyArray_CanCastSafely(PyArray_TYPE($input),PyArray_CDOUBLE)) ? 1 : 0;
};

%typemap(typecheck) npy_clongdouble *, const npy_clongdouble *, npy_clongdouble [], const npy_clongdouble []
{
  $1 = (is_array($input) && PyArray_CanCastSafely(PyArray_TYPE($input),PyArray_CLONGDOUBLE)) ? 1 : 0;
};






/*
 * IN types
 */
%apply int * IN_ARRAY1 {
    const int Ap [ ],
    const int Ai [ ],
    const int Aj [ ],
    const int Bp [ ],
    const int Bi [ ],	
    const int Bj [ ],
    const int offsets [ ]
};

%apply float * IN_ARRAY1 {
    const float Ax [ ],
    const float Bx [ ],
    const float Xx [ ],
    const float Yx [ ]
};

%apply double * IN_ARRAY1 {
    const double Ax [ ],
    const double Bx [ ],
    const double Xx [ ],
    const double Yx [ ]
};

%apply long double * IN_ARRAY1 {
    const long double Ax [ ],
    const long double Bx [ ],
    const long double Xx [ ],
    const long double Yx [ ]
};

%apply npy_cfloat * IN_ARRAY1 {
    const npy_cfloat Ax [ ],
    const npy_cfloat Bx [ ],
    const npy_cfloat Xx [ ],
    const npy_cfloat Yx [ ]
};

%apply npy_cdouble * IN_ARRAY1 {
    const npy_cdouble Ax [ ],
    const npy_cdouble Bx [ ],
    const npy_cdouble Xx [ ],
    const npy_cdouble Yx [ ]
};

%apply npy_clongdouble * IN_ARRAY1 {
    const npy_clongdouble Ax [ ],
    const npy_clongdouble Bx [ ],
    const npy_clongdouble Xx [ ],
    const npy_clongdouble Yx [ ]
};


%apply float * IN_ARRAY2 { const float Mx[] }
%apply double * IN_ARRAY2 { const double Mx[] }
%apply long double * IN_ARRAY2 { const long double Mx[] }
%apply npy_cfloat * IN_ARRAY2 { const npy_cfloat Mx[] }
%apply npy_cdouble * IN_ARRAY2 { const npy_cdouble Mx[] }
%apply npy_clongdouble * IN_ARRAY2 { const npy_longdouble Mx[] }


%apply float * IN_ARRAY2 { const float diags[] }
%apply double * IN_ARRAY2 { const double diags[] }
%apply long double * IN_ARRAY2 { const long double diags[] }
%apply npy_cfloat * IN_ARRAY2 { const npy_cfloat diags[] }
%apply npy_cdouble * IN_ARRAY2 { const npy_cdouble diags[] }
%apply npy_clongdouble * IN_ARRAY2 { const npy_longdouble diags[] }



/*
 * OUT types
 */
VEC_ARRAY_ARGOUT( int, INT )
%apply std::vector<int>* array_argout {
    std::vector<int>* Ap,
    std::vector<int>* Ai,
    std::vector<int>* Aj,
    std::vector<int>* Bp,
    std::vector<int>* Bi,
    std::vector<int>* Bj,
    std::vector<int>* Cp,
    std::vector<int>* Ci,
    std::vector<int>* Cj
};

VEC_ARRAY_ARGOUT( float, FLOAT )
%apply std::vector<float>* array_argout {
    std::vector<float>* Ax,
    std::vector<float>* Bx,
    std::vector<float>* Cx,
    std::vector<float>* Xx,
    std::vector<float>* Yx
};

VEC_ARRAY_ARGOUT( double, DOUBLE )
%apply std::vector<double>* array_argout {
    std::vector<double>* Ax,
    std::vector<double>* Bx,
    std::vector<double>* Cx,
    std::vector<double>* Xx,
    std::vector<double>* Yx
};


VEC_ARRAY_ARGOUT( long double, LONGDOUBLE )
%apply std::vector<long double>* array_argout {
    std::vector<long double>* Ax,
    std::vector<long double>* Bx,
    std::vector<long double>* Cx,
    std::vector<long double>* Xx,
    std::vector<long double>* Yx
};

VEC_ARRAY_ARGOUT( npy_cfloat, CFLOAT )
%apply std::vector<npy_cfloat>* array_argout {
    std::vector<npy_cfloat>* Ax,
    std::vector<npy_cfloat>* Bx,
    std::vector<npy_cfloat>* Cx,
    std::vector<npy_cfloat>* Xx,
    std::vector<npy_cfloat>* Yx
};


VEC_ARRAY_ARGOUT( npy_cdouble, CDOUBLE )
%apply std::vector<npy_cdouble>* array_argout {
    std::vector<npy_cdouble>* Ax,
    std::vector<npy_cdouble>* Bx,
    std::vector<npy_cdouble>* Cx,
    std::vector<npy_cdouble>* Xx,
    std::vector<npy_cdouble>* Yx
};


VEC_ARRAY_ARGOUT( npy_clongdouble, CLONGDOUBLE )
%apply std::vector<npy_clongdouble>* array_argout {
    std::vector<npy_clongdouble>* Ax,
    std::vector<npy_clongdouble>* Bx,
    std::vector<npy_clongdouble>* Cx,
    std::vector<npy_clongdouble>* Xx,
    std::vector<npy_clongdouble>* Yx
};



/*
 * INOUT types
 */
%apply float * INPLACE_ARRAY2 { float Mx [] }
%apply double * INPLACE_ARRAY2 { double Mx [] }
%apply long double * INPLACE_ARRAY2 { long double Mx[] }
%apply npy_cfloat * INPLACE_ARRAY2 { npy_cfloat Mx[] }
%apply npy_cdouble * INPLACE_ARRAY2 { npy_cdouble Mx[] }
%apply npy_clongdouble * INPLACE_ARRAY2 { npy_longdouble Mx[] }






%include "sparsetools.h"



 /*
  * Order may be important here, list float before double
  */

/*
 *  CSR->CSC or CSC->CSR or CSR = CSR^T or CSC = CSC^T
 */
%template(csrtocsc)   csrtocsc<float>; 
%template(csrtocsc)   csrtocsc<double>; 
%template(csrtocsc)   csrtocsc<long double>; 
%template(csrtocsc)   csrtocsc<npy_cfloat>; 
%template(csrtocsc)   csrtocsc<npy_cdouble>; 
%template(csrtocsc)   csrtocsc<npy_clongdouble>; 

%template(csctocsr)   csctocsr<float>; 
%template(csctocsr)   csctocsr<double>; 
%template(csctocsr)   csctocsr<long double>; 
%template(csctocsr)   csctocsr<npy_cfloat>; 
%template(csctocsr)   csctocsr<npy_cdouble>; 
%template(csctocsr)   csctocsr<npy_clongdouble>; 


/*
 * CSR<->COO and CSC<->COO
 */
%template(csrtocoo)   csrtocoo<float>; 
%template(csrtocoo)   csrtocoo<double>; 
%template(csrtocoo)   csrtocoo<long double>; 
%template(csrtocoo)   csrtocoo<npy_cfloat>; 
%template(csrtocoo)   csrtocoo<npy_cdouble>; 
%template(csrtocoo)   csrtocoo<npy_clongdouble>; 

%template(cootocsr)   cootocsr<float>; 
%template(cootocsr)   cootocsr<double>; 
%template(cootocsr)   cootocsr<long double>; 
%template(cootocsr)   cootocsr<npy_cfloat>; 
%template(cootocsr)   cootocsr<npy_cdouble>; 
%template(cootocsr)   cootocsr<npy_clongdouble>; 

%template(csctocoo)   csctocoo<float>; 
%template(csctocoo)   csctocoo<double>; 
%template(csctocoo)   csctocoo<long double>; 
%template(csctocoo)   csctocoo<npy_cfloat>; 
%template(csctocoo)   csctocoo<npy_cdouble>; 
%template(csctocoo)   csctocoo<npy_clongdouble>; 

%template(cootocsc)   cootocsc<float>; 
%template(cootocsc)   cootocsc<double>; 
%template(cootocsc)   cootocsc<long double>; 
%template(cootocsc)   cootocsc<npy_cfloat>; 
%template(cootocsc)   cootocsc<npy_cdouble>; 
%template(cootocsc)   cootocsc<npy_clongdouble>; 


/*
 * CSR+CSR and CSC+CSC
 */
%template(csrplcsr)   csrplcsr<float>; 
%template(csrplcsr)   csrplcsr<double>; 
%template(csrplcsr)   csrplcsr<long double>; 
%template(csrplcsr)   csrplcsr<npy_cfloat>; 
%template(csrplcsr)   csrplcsr<npy_cdouble>; 
%template(csrplcsr)   csrplcsr<npy_clongdouble>; 

%template(cscplcsc)   cscplcsc<float>; 
%template(cscplcsc)   cscplcsc<double>; 
%template(cscplcsc)   cscplcsc<long double>; 
%template(cscplcsc)   cscplcsc<npy_cfloat>; 
%template(cscplcsc)   cscplcsc<npy_cdouble>; 
%template(cscplcsc)   cscplcsc<npy_clongdouble>; 



/*
 * CSR*CSR and CSC*CSC
 */
%template(csrmucsr)   csrmucsr<float>; 
%template(csrmucsr)   csrmucsr<double>; 
%template(csrmucsr)   csrmucsr<long double>; 
%template(csrmucsr)   csrmucsr<npy_cfloat>; 
%template(csrmucsr)   csrmucsr<npy_cdouble>; 
%template(csrmucsr)   csrmucsr<npy_clongdouble>; 

%template(cscmucsc)   cscmucsc<float>; 
%template(cscmucsc)   cscmucsc<double>; 
%template(cscmucsc)   cscmucsc<long double>; 
%template(cscmucsc)   cscmucsc<npy_cfloat>; 
%template(cscmucsc)   cscmucsc<npy_cdouble>; 
%template(cscmucsc)   cscmucsc<npy_clongdouble>; 

/*
 * CSR*x and CSC*x
 */
%template(csrmux)   csrmux<float>; 
%template(csrmux)   csrmux<double>; 
%template(csrmux)   csrmux<long double>; 
%template(csrmux)   csrmux<npy_cfloat>; 
%template(csrmux)   csrmux<npy_cdouble>; 
%template(csrmux)   csrmux<npy_clongdouble>; 

%template(cscmux)   cscmux<float>; 
%template(cscmux)   cscmux<double>; 
%template(cscmux)   cscmux<long double>; 
%template(cscmux)   cscmux<npy_cfloat>; 
%template(cscmux)   cscmux<npy_cdouble>; 
%template(cscmux)   cscmux<npy_clongdouble>; 


/*
 * CSR(elmul)CSR and CSC(elmul)CSC
 */
%template(csrelmulcsr)   csrelmulcsr<float>; 
%template(csrelmulcsr)   csrelmulcsr<double>; 
%template(csrelmulcsr)   csrelmulcsr<long double>; 
%template(csrelmulcsr)   csrelmulcsr<npy_cfloat>; 
%template(csrelmulcsr)   csrelmulcsr<npy_cdouble>; 
%template(csrelmulcsr)   csrelmulcsr<npy_clongdouble>; 

%template(cscelmulcsc)   cscelmulcsc<float>; 
%template(cscelmulcsc)   cscelmulcsc<double>; 
%template(cscelmulcsc)   cscelmulcsc<long double>; 
%template(cscelmulcsc)   cscelmulcsc<npy_cfloat>; 
%template(cscelmulcsc)   cscelmulcsc<npy_cdouble>; 
%template(cscelmulcsc)   cscelmulcsc<npy_clongdouble>; 


/*
 * spdiags->CSC
 */
%template(spdiags)   spdiags<float>; 
%template(spdiags)   spdiags<double>; 
%template(spdiags)   spdiags<long double>; 
%template(spdiags)   spdiags<npy_cfloat>; 
%template(spdiags)   spdiags<npy_cdouble>; 
%template(spdiags)   spdiags<npy_clongdouble>; 


/*
 * CSR<->Dense
 */
%template(csrtodense)   csrtodense<float>; 
%template(csrtodense)   csrtodense<double>; 
%template(csrtodense)   csrtodense<long double>; 
%template(csrtodense)   csrtodense<npy_cfloat>; 
%template(csrtodense)   csrtodense<npy_cdouble>; 
%template(csrtodense)   csrtodense<npy_clongdouble>; 

%template(densetocsr)   densetocsr<float>; 
%template(densetocsr)   densetocsr<double>; 
%template(densetocsr)   densetocsr<long double>; 
%template(densetocsr)   densetocsr<npy_cfloat>; 
%template(densetocsr)   densetocsr<npy_cdouble>; 
%template(densetocsr)   densetocsr<npy_clongdouble>; 
