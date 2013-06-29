/* -*- C++ -*- */
/*%module sparsetools*/

 /* why does SWIG complain about int arrays? a typecheck is provided */
#pragma SWIG nowarn=467

%{
#include "py3k.h"
#define SWIG_FILE_WITH_INIT
#include "Python.h"
#include "numpy/arrayobject.h"
#include "complex_ops.h"
#include "bool_ops.h"
/*#include "sparsetools.h"*/
%}

%feature("autodoc", "1");

%include "numpy.i"

%init %{
    import_array();
%}



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
    const ctype Cp [ ],
    const ctype Ci [ ],	
    const ctype Cj [ ],
    const ctype offsets [ ]
};
%enddef

%define T_IN_ARRAY1( ctype )
%apply ctype * IN_ARRAY1 {
    const ctype Ax [ ],
    const ctype Bx [ ],
    const ctype Cx [ ],
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


 /*
  * OUT types
  */
%define I_ARRAY_ARGOUT( ctype )
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

%define T_ARRAY_ARGOUT( ctype )
%apply std::vector<ctype>* array_argout {
    std::vector<ctype>* Ax, 
    std::vector<ctype>* Bx,
    std::vector<ctype>* Cx, 
    std::vector<ctype>* Xx,
    std::vector<ctype>* Yx 
};
%enddef




 /*
  * INOUT types
  */
%define I_INPLACE_ARRAY1( ctype )
%apply ctype * INPLACE_ARRAY {
  ctype Ap [ ],
  ctype Ai [ ],
  ctype Aj [ ],
  ctype Bp [ ],
  ctype Bi [ ],
  ctype Bj [ ],
  ctype Cp [ ],
  ctype Ci [ ],
  ctype Cj [ ],
  ctype flag [ ]
};
%enddef

%define T_INPLACE_ARRAY1( ctype )
%apply ctype * INPLACE_ARRAY {
  ctype Ax [ ],
  ctype Bx [ ],
  ctype Cx [ ],
  ctype Yx [ ]
};
%enddef

%define T_INPLACE_ARRAY2( ctype )
%apply ctype * INPLACE_ARRAY2 {
  ctype Mx [ ]
};
%enddef


/*
 * Macros to instantiate index types and data types
 */
%define DECLARE_INDEX_TYPE( ctype )
I_IN_ARRAY1( ctype )
I_ARRAY_ARGOUT( ctype )
I_INPLACE_ARRAY1( ctype )
%enddef

%define DECLARE_DATA_TYPE( ctype )
T_IN_ARRAY1( ctype )
T_IN_ARRAY2( ctype )
T_ARRAY_ARGOUT( ctype )
T_INPLACE_ARRAY1( ctype )
T_INPLACE_ARRAY2( ctype )
%enddef


/*
 * Create all desired index and data types here
 */
DECLARE_INDEX_TYPE( int       )

DECLARE_DATA_TYPE( npy_bool_wrapper        )
DECLARE_DATA_TYPE( signed char             )
DECLARE_DATA_TYPE( unsigned char           )
DECLARE_DATA_TYPE( short                   )
DECLARE_DATA_TYPE( unsigned short          )
DECLARE_DATA_TYPE( int                     )
DECLARE_DATA_TYPE( unsigned int            )
DECLARE_DATA_TYPE( long long               )
DECLARE_DATA_TYPE( unsigned long long      )
DECLARE_DATA_TYPE( float                   )
DECLARE_DATA_TYPE( double                  )
DECLARE_DATA_TYPE( long double             )
DECLARE_DATA_TYPE( npy_cfloat_wrapper      )
DECLARE_DATA_TYPE( npy_cdouble_wrapper     )
DECLARE_DATA_TYPE( npy_clongdouble_wrapper )




/*%include "sparsetools.h"*/

 /*
  * Order is important here, list int before float, float before 
  * double, scalar before complex, etc.
  */

%define INSTANTIATE_ALL( f_name )
/* 32-bit indices */
%template(f_name)   f_name<int,npy_bool_wrapper>;
%template(f_name)   f_name<int,signed char>;
%template(f_name)   f_name<int,unsigned char>;
%template(f_name)   f_name<int,short>;
%template(f_name)   f_name<int,unsigned short>;
%template(f_name)   f_name<int,int>;
%template(f_name)   f_name<int,unsigned int>;
%template(f_name)   f_name<int,long long>;
%template(f_name)   f_name<int,unsigned long long>;
%template(f_name)   f_name<int,float>;
%template(f_name)   f_name<int,double>;
%template(f_name)   f_name<int,long double>;
%template(f_name)   f_name<int,npy_cfloat_wrapper>;
%template(f_name)   f_name<int,npy_cdouble_wrapper>;
%template(f_name)   f_name<int,npy_clongdouble_wrapper>;
/* 64-bit indices would go here */
%enddef


%define INSTANTIATE_INDEX( f_name )
/* 32-bit indices */
%template(f_name)   f_name<int>;
/* 64-bit indices would go here */
%enddef

%define INSTANTIATE_BOOL_OUT( f_name )
/* 32-bit indices */
%template(f_name)   f_name<int,npy_bool_wrapper,        npy_bool_wrapper>;
%template(f_name)   f_name<int,signed char,             npy_bool_wrapper>;
%template(f_name)   f_name<int,unsigned char,           npy_bool_wrapper>;
%template(f_name)   f_name<int,short,                   npy_bool_wrapper>;
%template(f_name)   f_name<int,unsigned short,          npy_bool_wrapper>;
%template(f_name)   f_name<int,int,                     npy_bool_wrapper>;
%template(f_name)   f_name<int,unsigned int,            npy_bool_wrapper>;
%template(f_name)   f_name<int,long long,               npy_bool_wrapper>;
%template(f_name)   f_name<int,unsigned long long,      npy_bool_wrapper>;
%template(f_name)   f_name<int,float,                   npy_bool_wrapper>;
%template(f_name)   f_name<int,double,                  npy_bool_wrapper>;
%template(f_name)   f_name<int,long double,             npy_bool_wrapper>;
%template(f_name)   f_name<int,npy_cfloat_wrapper,      npy_bool_wrapper>;
%template(f_name)   f_name<int,npy_cdouble_wrapper,     npy_bool_wrapper>;
%template(f_name)   f_name<int,npy_clongdouble_wrapper, npy_bool_wrapper>;
/* 64-bit indices would go here */
%enddef
