#ifndef _VQ_H_
#define _VQ_H

int double_tvq(double* obs, double* code_book, int Nobs, int Ncodes, 
        int Nfeatures, long long* codes, double* lowest_dist);

int float_tvq(float* obs, float* code_book, int Nobs, int Ncodes, 
        int Nfeatures, long long* codes, float* lowest_dist);

#endif
