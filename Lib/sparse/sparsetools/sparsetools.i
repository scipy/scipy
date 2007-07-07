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
%} //end inline code



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



 /*
  * make typechecks - used for overloading
  */
%include "typemaps.i"

%define NPY_TYPECHECK( ctype, atype )
%typemap(typecheck) ctype *, const ctype *, ctype [], const ctype []
{
  $1 = (is_array($input) && PyArray_CanCastSafely(PyArray_TYPE($input),PyArray_##atype)) ? 1 : 0;
};
%enddef

NPY_TYPECHECK(         int,      INT )
NPY_TYPECHECK(        long,     LONG )
NPY_TYPECHECK(       float,    FLOAT )
NPY_TYPECHECK(      double,   DOUBLE )
NPY_TYPECHECK(  npy_cfloat,   CFLOAT )
NPY_TYPECHECK(  npy_cdouble, CDOUBLE )



 /*
  * IN types
  */
%define I_IN_ARRAY1( ctype )
%apply ctype * IN_ARRAY1 {
    const ctype Ap [ ],
    const ctype Ai [ ],
    const ctype Aj [ ],
    const ctype Bp [ ],
    const ctype Bi [ ],	
    const ctype Bj [ ],
    const ctype offsets [ ]
};
%enddef

%define T_IN_ARRAY1( ctype )
%apply ctype * IN_ARRAY1 {
    const ctype Ax [ ],
    const ctype Bx [ ],
    const ctype Xx [ ],
    const ctype Yx [ ]
};
%enddef

%define T_IN_ARRAY2( ctype )
%apply ctype * IN_ARRAY2 {
  const ctype Mx    [ ],
  const ctype diags [ ]
};
%enddef


I_IN_ARRAY1( int  )
I_IN_ARRAY1( long )

T_IN_ARRAY1( int         )
T_IN_ARRAY1( long        )
T_IN_ARRAY1( float       )
T_IN_ARRAY1( double      )
T_IN_ARRAY1( npy_cfloat  )
T_IN_ARRAY1( npy_cdouble )

T_IN_ARRAY2( int         )
T_IN_ARRAY2( long        )
T_IN_ARRAY2( float       )
T_IN_ARRAY2( double      )
T_IN_ARRAY2( npy_cfloat  )
T_IN_ARRAY2( npy_cdouble )



 /*
  * OUT types
  */
%define I_ARRAY_ARGOUT( ctype, atype )
VEC_ARRAY_ARGOUT( ctype, atype )
%apply std::vector<ctype>* array_argout {
    std::vector<ctype>* Ap,
    std::vector<ctype>* Ai,
    std::vector<ctype>* Aj,
    std::vector<ctype>* Bp,
    std::vector<ctype>* Bi,
    std::vector<ctype>* Bj,
    std::vector<ctype>* Cp,
    std::vector<ctype>* Ci,
    std::vector<ctype>* Cj
};
%enddef

%define T_ARRAY_ARGOUT( ctype, atype )
VEC_ARRAY_ARGOUT( ctype, atype )
%apply std::vector<ctype>* array_argout {
    std::vector<ctype>* Ax, 
    std::vector<ctype>* Bx,
    std::vector<ctype>* Cx, 
    std::vector<ctype>* Xx,
    std::vector<ctype>* Yx 
};
%enddef



I_ARRAY_ARGOUT( int,   INT)
I_ARRAY_ARGOUT( long, LONG)

T_ARRAY_ARGOUT( int,         INT     )
T_ARRAY_ARGOUT( long,        LONG    )
T_ARRAY_ARGOUT( float,       FLOAT   )
T_ARRAY_ARGOUT( double,      DOUBLE  )
T_ARRAY_ARGOUT( npy_cfloat,  CFLOAT  )
T_ARRAY_ARGOUT( npy_cdouble, CDOUBLE )



 /*
  * INOUT types
  */
%define T_INPLACE_ARRAY2( ctype )
%apply ctype * INPLACE_ARRAY2 {
  ctype Mx [ ]
};
%enddef

T_INPLACE_ARRAY2( int         )
T_INPLACE_ARRAY2( long        )
T_INPLACE_ARRAY2( float       )
T_INPLACE_ARRAY2( double      )
T_INPLACE_ARRAY2( npy_cfloat  )
T_INPLACE_ARRAY2( npy_cdouble )



%define I_INPLACE_ARRAY1( ctype )
%apply ctype * INPLACE_ARRAY {
  ctype Aj [ ]
};
%enddef

I_INPLACE_ARRAY1( int         )
I_INPLACE_ARRAY1( long        )


%define T_INPLACE_ARRAY1( ctype )
%apply ctype * INPLACE_ARRAY {
  ctype Ax [ ]
};
%enddef

T_INPLACE_ARRAY1( long        )
T_INPLACE_ARRAY1( float       )
T_INPLACE_ARRAY1( double      )
T_INPLACE_ARRAY1( npy_cfloat  )
T_INPLACE_ARRAY1( npy_cdouble )






%include "sparsetools.h"
 /*
  * Order may be important here, list float before double, scalar before complex
  */

%define INSTANTIATE_ALL( f_name )		     
%template(f_name)   f_name<int,int>;
%template(f_name)   f_name<int,long>;
%template(f_name)   f_name<int,float>;
%template(f_name)   f_name<int,double>;
%template(f_name)   f_name<int,npy_cfloat>;
%template(f_name)   f_name<int,npy_cdouble>;
/* 64-bit indices would go here */
%enddef



/*
 *  CSR->CSC or CSC->CSR or CSR = CSR^T or CSC = CSC^T
 */
INSTANTIATE_ALL(csrtocsc)
INSTANTIATE_ALL(csctocsr)

/*
 * CSR<->COO and CSC<->COO
 */
INSTANTIATE_ALL(csrtocoo)
INSTANTIATE_ALL(csctocoo)
INSTANTIATE_ALL(cootocsr)
INSTANTIATE_ALL(cootocsc)

/*
 * CSR+CSR and CSC+CSC
 */
INSTANTIATE_ALL(csrplcsr)
INSTANTIATE_ALL(cscplcsc)

/*
 * CSR*CSR and CSC*CSC
 */
INSTANTIATE_ALL(csrmucsr)
INSTANTIATE_ALL(cscmucsc)

/*
 * CSR*x and CSC*x
 */
INSTANTIATE_ALL(csrmux)
INSTANTIATE_ALL(cscmux)

/*
 * CSR(elmul)CSR and CSC(elmul)CSC
 */
INSTANTIATE_ALL(csrelmulcsr)
INSTANTIATE_ALL(cscelmulcsc)


/*
 * spdiags->CSC
 */
INSTANTIATE_ALL(spdiags)

/*
 * CSR<->Dense
 */
INSTANTIATE_ALL(csrtodense)
INSTANTIATE_ALL(densetocsr)

/*
 * Ensure sorted CSR/CSC indices.
 */
INSTANTIATE_ALL(sort_csr_indices)
INSTANTIATE_ALL(sort_csc_indices)

