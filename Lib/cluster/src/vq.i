%module _vq
%{

#include "vq.h"

/* Wrappers for the template code */

void float_vq(float * obs,float* code_book, int Nobs, int Ncodes, 
              int Nfeatures, int* codes, float* lowest_dist)
{
    tvq<float>(obs,code_book,Nobs,Ncodes,Nfeatures,codes,lowest_dist);
}

void double_vq(double * obs,double* code_book, int Nobs, int Ncodes, 
              int Nfeatures, int* codes, double* lowest_dist)
{
    tvq<double>(obs,code_book,Nobs,Ncodes,Nfeatures,codes,lowest_dist);
}

%}

%include swig_num.i

void double_vq(double_IN_D0_D2 *obs,double_IN_D1_D2 *code_book, 
               int DIM0, int DIM1, int DIM2, 
               int_ARGOUT_TUPLE_D0 *codes, double_ARGOUT_TUPLE_D0 *lowest_dist);

void float_vq(float_IN_D0_D2 *obs,float_IN_D1_D2 *code_book, 
              int DIM0, int DIM1, int DIM2, 
              int_ARGOUT_TUPLE_D0 *codes, float_ARGOUT_TUPLE_D0 *lowest_dist);
