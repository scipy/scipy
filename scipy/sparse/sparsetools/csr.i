%module csr

%include "sparsetools.i"

%{
#include "csr.h"
%}

%include "csr.h" 


%template(expandptr)   expandptr<int>;
%template(csr_count_blocks)   csr_count_blocks<int>;
%template(csr_matmat_pass1)   csr_matmat_pass1<int>;
%template(csr_has_sorted_indices)   csr_has_sorted_indices<int>;


INSTANTIATE_ALL(csr_diagonal)
INSTANTIATE_ALL(csr_scale_rows)
INSTANTIATE_ALL(csr_scale_columns)
INSTANTIATE_ALL(csr_tocsc)
INSTANTIATE_ALL(csr_tobsr)
INSTANTIATE_ALL(csr_matmat_pass2)
INSTANTIATE_ALL(csr_matvec)
INSTANTIATE_ALL(csr_matvecs)
INSTANTIATE_ALL(csr_elmul_csr)
INSTANTIATE_ALL(csr_eldiv_csr)
INSTANTIATE_ALL(csr_plus_csr)
INSTANTIATE_ALL(csr_minus_csr)
INSTANTIATE_ALL(csr_sort_indices)
INSTANTIATE_ALL(csr_eliminate_zeros)
INSTANTIATE_ALL(csr_sum_duplicates)
INSTANTIATE_ALL(get_csr_submatrix)

