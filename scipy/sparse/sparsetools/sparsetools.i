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
  ctype Aj [ ]
};
%enddef

%define T_INPLACE_ARRAY1( ctype )
%apply ctype * INPLACE_ARRAY {
  ctype Ax [ ]
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
DECLARE_INDEX_TYPE( long long )

DECLARE_DATA_TYPE( signed char         )
DECLARE_DATA_TYPE( unsigned char       )
DECLARE_DATA_TYPE( short               )
DECLARE_DATA_TYPE( int                 )
DECLARE_DATA_TYPE( long long           )
DECLARE_DATA_TYPE( float               )
DECLARE_DATA_TYPE( double              )
DECLARE_DATA_TYPE( npy_cfloat_wrapper  )
DECLARE_DATA_TYPE( npy_cdouble_wrapper )




%include "sparsetools.h"
 /*
  * Order may be important here, list float before double, scalar before complex
  * 
  * Should we permit unsigned types as array indices?  
  * Do any functions require signedness? -- Nathan (Aug 2007)
  */

%define INSTANTIATE_ALL( f_name )
/* 64-bit indices */
%template(f_name)   f_name<int,signed char>;
%template(f_name)   f_name<int,unsigned char>;
%template(f_name)   f_name<int,short>;
%template(f_name)   f_name<int,int>;
%template(f_name)   f_name<int,long long>;
%template(f_name)   f_name<int,float>;
%template(f_name)   f_name<int,double>;
%template(f_name)   f_name<int,npy_cfloat_wrapper>;
%template(f_name)   f_name<int,npy_cdouble_wrapper>;
/* 64-bit indices */
%template(f_name)   f_name<long long,signed char>;
%template(f_name)   f_name<long long,unsigned char>;
%template(f_name)   f_name<long long,short>;
%template(f_name)   f_name<long long,int>;
%template(f_name)   f_name<long long,long long>;
%template(f_name)   f_name<long long,float>;
%template(f_name)   f_name<long long,double>;
%template(f_name)   f_name<long long,npy_cfloat_wrapper>;
%template(f_name)   f_name<long long,npy_cdouble_wrapper>;
%enddef


/*
 *  diag(CSR) and diag(CSC)
 */
INSTANTIATE_ALL(extract_csr_diagonal)
INSTANTIATE_ALL(extract_csc_diagonal)


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
 * CSR (binary op) CSR and CSC (binary op) CSC
 */
INSTANTIATE_ALL(csr_elmul_csr)
INSTANTIATE_ALL(csr_eldiv_csr)
INSTANTIATE_ALL(csr_plus_csr)
INSTANTIATE_ALL(csr_minus_csr)

INSTANTIATE_ALL(csc_elmul_csc)
INSTANTIATE_ALL(csc_eldiv_csc)
INSTANTIATE_ALL(csc_plus_csc)
INSTANTIATE_ALL(csc_minus_csc)



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
 * Sort CSR/CSC indices.
 */
INSTANTIATE_ALL(sort_csr_indices)
INSTANTIATE_ALL(sort_csc_indices)


/*
 * Sum duplicate CSR/CSC entries.
 */
INSTANTIATE_ALL(sum_csr_duplicates)
INSTANTIATE_ALL(sum_csc_duplicates)

INSTANTIATE_ALL(get_csr_submatrix)
