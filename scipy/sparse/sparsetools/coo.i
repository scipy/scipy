/* -*- C++ -*- */
%module coo

%include "sparsetools.i"

%{
#include "coo.h"
%}

%include "coo.h" 

INSTANTIATE_ALL(coo_tocsr)
INSTANTIATE_ALL(coo_tocsc)
INSTANTIATE_ALL(coo_todense)

INSTANTIATE_ALL(coo_matvec)

INSTANTIATE_INDEX(coo_count_diagonals)


