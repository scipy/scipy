/* -*- C++ -*- */
%module bsr

%include "sparsetools.i"

%{
#include "bsr.h"
%}

%include "bsr.h" 

INSTANTIATE_ALL(bsr_diagonal)
INSTANTIATE_ALL(bsr_scale_rows)
INSTANTIATE_ALL(bsr_scale_columns)
INSTANTIATE_ALL(bsr_transpose)
INSTANTIATE_ALL(bsr_matmat_pass2)
INSTANTIATE_ALL(bsr_matvec)
INSTANTIATE_ALL(bsr_matvecs)
INSTANTIATE_ALL(bsr_elmul_bsr)
INSTANTIATE_ALL(bsr_eldiv_bsr)
INSTANTIATE_ALL(bsr_plus_bsr)
INSTANTIATE_ALL(bsr_minus_bsr)
INSTANTIATE_ALL(bsr_sort_indices)

INSTANTIATE_BOOL_OUT(bsr_ne_bsr)
INSTANTIATE_BOOL_OUT(bsr_lt_bsr)
INSTANTIATE_BOOL_OUT(bsr_gt_bsr)
INSTANTIATE_BOOL_OUT(bsr_le_bsr)
INSTANTIATE_BOOL_OUT(bsr_ge_bsr)
