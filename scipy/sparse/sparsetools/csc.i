/* -*- C++ -*- */
%module csc

%include "sparsetools.i"

%{
#include "csc.h"
%}

%include "csc.h" 

INSTANTIATE_INDEX(csc_matmat_pass1);

INSTANTIATE_ALL(csc_diagonal)
INSTANTIATE_ALL(csc_tocsr)
INSTANTIATE_ALL(csc_matmat_pass2)
INSTANTIATE_ALL(csc_matvec)
INSTANTIATE_ALL(csc_matvecs)
INSTANTIATE_ALL(csc_elmul_csc)
INSTANTIATE_ALL(csc_eldiv_csc)
INSTANTIATE_ALL(csc_plus_csc)
INSTANTIATE_ALL(csc_minus_csc)

INSTANTIATE_BOOL_OUT(csc_ne_csc)
INSTANTIATE_BOOL_OUT(csc_lt_csc)
INSTANTIATE_BOOL_OUT(csc_gt_csc)
INSTANTIATE_BOOL_OUT(csc_le_csc)
INSTANTIATE_BOOL_OUT(csc_ge_csc)
