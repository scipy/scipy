/* -*- C++ -*- */
%module csr

%include "sparsetools.i"

%{
#include "csr.h"
%}

%include "csr.h" 

INSTANTIATE_INDEX(expandptr)
INSTANTIATE_INDEX(csr_matmat_pass1)
INSTANTIATE_INDEX(csr_count_blocks)
INSTANTIATE_INDEX(csr_has_sorted_indices)

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
INSTANTIATE_ALL(csr_sample_values)

INSTANTIATE_BOOL_OUT(csr_ne_csr)
INSTANTIATE_BOOL_OUT(csr_lt_csr)
INSTANTIATE_BOOL_OUT(csr_gt_csr)
INSTANTIATE_BOOL_OUT(csr_le_csr)
INSTANTIATE_BOOL_OUT(csr_ge_csr)
