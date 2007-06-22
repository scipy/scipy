#ifndef _VQ_H_
#define _VQ_H

#include <Python.h>

#include <numpy/arrayobject.h>

int double_tvq(double* obs, double* code_book, int Nobs, int Ncodes, 
        int Nfeatures, npy_intp* codes, double* lowest_dist);

int float_tvq(float* obs, float* code_book, int Nobs, int Ncodes, 
        int Nfeatures, npy_intp* codes, float* lowest_dist);

#endif
