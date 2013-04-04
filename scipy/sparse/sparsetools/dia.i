/* -*- C++ -*- */
%module dia

%include "sparsetools.i"

%{
#include "dia.h"
%}

%include "dia.h" 

INSTANTIATE_ALL(dia_matvec)

