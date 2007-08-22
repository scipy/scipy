/* -*- C -*-  (not really, but good for syntax highlighting) */
%module multigridtools

 /* why does SWIG complain about int arrays? a typecheck is provided */
#pragma SWIG nowarn=467

%{
#define SWIG_FILE_WITH_INIT
#include "numpy/arrayobject.h"

#include "ruge_stuben.h"
#include "smoothed_aggregation.h"
#include "relaxation.h"

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
    const ctype Sp [ ],
    const ctype Si [ ],	
    const ctype Sj [ ],
    const ctype Tp [ ],
    const ctype Ti [ ],	
    const ctype Tj [ ]
};
%enddef

%define T_IN_ARRAY1( ctype )
%apply ctype * IN_ARRAY1 {
    const ctype Ax [ ],
    const ctype Bx [ ],
    const ctype Sx [ ],
    const ctype Tx [ ],
    const ctype Xx [ ],
    const ctype Yx [ ],
    const ctype  x [ ],
    const ctype  y [ ],
    const ctype  b [ ]    
};
%enddef


I_IN_ARRAY1( int         )
T_IN_ARRAY1( float       )
T_IN_ARRAY1( double      )



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
    std::vector<ctype>* Cj,
    std::vector<ctype>* Sp,
    std::vector<ctype>* Si,
    std::vector<ctype>* Sj,
    std::vector<ctype>* Tp,
    std::vector<ctype>* Ti,
    std::vector<ctype>* Tj
};
%enddef

%define T_ARRAY_ARGOUT( ctype )
%apply std::vector<ctype>* array_argout {
    std::vector<ctype>* Ax, 
    std::vector<ctype>* Bx,
    std::vector<ctype>* Cx, 
    std::vector<ctype>* Sx,
    std::vector<ctype>* Tx, 
    std::vector<ctype>* Xx,
    std::vector<ctype>* Yx 
};
%enddef

I_ARRAY_ARGOUT( int    )
T_ARRAY_ARGOUT( float  )
T_ARRAY_ARGOUT( double )



 /*
  * INPLACE types
  */
%define I_INPLACE_ARRAY1( ctype )
%apply ctype * INPLACE_ARRAY {
  ctype Aj [ ]
};
%enddef

%define T_INPLACE_ARRAY1( ctype )
%apply ctype * INPLACE_ARRAY {
  ctype    x [ ],
  ctype temp [ ]
};
%enddef

I_INPLACE_ARRAY1( int    )
T_INPLACE_ARRAY1( float  )
T_INPLACE_ARRAY1( double )


%include "ruge_stuben.h"
%include "smoothed_aggregation.h"
%include "relaxation.h"

 /*
  * Order may be important here, list float before double
  */

%define INSTANTIATE_BOTH( f_name )
%template(f_name)   f_name<int,float>;
%template(f_name)   f_name<int,double>;
/* 64-bit indices would go here */
%enddef
 
%define INSTANTIATE_INDEX( f_name )
%template(f_name)   f_name<int>;
/* 64-bit indices would go here */
%enddef

%define INSTANTIATE_DATA( f_name )
%template(f_name)   f_name<float>;
%template(f_name)   f_name<double>;
%enddef
 
 

INSTANTIATE_DATA(rs_strong_connections)
INSTANTIATE_DATA(rs_interpolation)

INSTANTIATE_DATA(sa_strong_connections)
INSTANTIATE_DATA(sa_smoother)
/*INSTANTIATE_INDEX(sa_get_aggregates)*/

INSTANTIATE_BOTH(gauss_seidel)
INSTANTIATE_BOTH(jacobi)

