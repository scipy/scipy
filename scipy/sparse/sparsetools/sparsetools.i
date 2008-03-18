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
  ctype Cj [ ]
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
%include "csr.h" 
%include "csc.h"
%include "coo.h"
%include "bsr.h"
%include "dia.h"

 /*
  * Order is important here, list int before float, float before 
  * double, scalar before complex, etc.
  */

%define INSTANTIATE_ALL( f_name )
/* 32-bit indices */
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




/*
 *  diag(A)
 */
INSTANTIATE_ALL(csr_diagonal)
INSTANTIATE_ALL(csc_diagonal)
INSTANTIATE_ALL(bsr_diagonal)


/*
 *  scale columns
 */
INSTANTIATE_ALL(csr_scale_rows)
INSTANTIATE_ALL(csr_scale_columns)
INSTANTIATE_ALL(bsr_scale_rows)
INSTANTIATE_ALL(bsr_scale_columns)


/*
 *  CSR->CSC or CSC->CSR or CSR = CSR^T or CSC = CSC^T or CSR->BSR
 */
INSTANTIATE_ALL(csr_tocsc)
INSTANTIATE_ALL(csc_tocsr)
INSTANTIATE_ALL(csr_tobsr)
INSTANTIATE_ALL(bsr_transpose)

/*
 * CSR<->COO and CSC<->COO
 */
%template(expandptr)   expandptr<int>;
INSTANTIATE_ALL(coo_tocsr)
INSTANTIATE_ALL(coo_tocsc)

/*
 * CSR<->BSR
 */
%template(csr_count_blocks)   csr_count_blocks<int>;


/*
 * CSR*CSR and CSC*CSC
 */
%template(csr_matmat_pass1)   csr_matmat_pass1<int>;
%template(csc_matmat_pass1)   csc_matmat_pass1<int>;
INSTANTIATE_ALL(csr_matmat_pass2)
INSTANTIATE_ALL(csc_matmat_pass2)
INSTANTIATE_ALL(bsr_matmat_pass2)


/*
 * A*x
 */
INSTANTIATE_ALL(csr_matvec)
INSTANTIATE_ALL(csc_matvec)
INSTANTIATE_ALL(bsr_matvec)
INSTANTIATE_ALL(dia_matvec)

/*
 * A (binary op) B 
 */
INSTANTIATE_ALL(csr_elmul_csr)
INSTANTIATE_ALL(csr_eldiv_csr)
INSTANTIATE_ALL(csr_plus_csr)
INSTANTIATE_ALL(csr_minus_csr)

INSTANTIATE_ALL(csc_elmul_csc)
INSTANTIATE_ALL(csc_eldiv_csc)
INSTANTIATE_ALL(csc_plus_csc)
INSTANTIATE_ALL(csc_minus_csc)

INSTANTIATE_ALL(bsr_elmul_bsr)
INSTANTIATE_ALL(bsr_eldiv_bsr)
INSTANTIATE_ALL(bsr_plus_bsr)
INSTANTIATE_ALL(bsr_minus_bsr)

/*
 * Sort indices
 */
%template(csr_has_sorted_indices)   csr_has_sorted_indices<int>;
INSTANTIATE_ALL(csr_sort_indices)
INSTANTIATE_ALL(bsr_sort_indices)


/*
 * Remove zeros
 */
INSTANTIATE_ALL(csr_eliminate_zeros)

/*
 * Sum duplicate entries
 */
INSTANTIATE_ALL(csr_sum_duplicates)

/*
 * Extract submatrices
 */
INSTANTIATE_ALL(get_csr_submatrix)

/*
 * To dense matrix
 */
INSTANTIATE_ALL(coo_todense)

