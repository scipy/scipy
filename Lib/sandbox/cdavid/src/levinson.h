#ifndef LEVINSON1D_H
#define LEVINSON1D_H

#include <stddef.h>

int dbl_levinson1d(const double* in, size_t order, 
		double* acoeff, double *err, double* kcoeff, double* tmp);

int dbl_levinson1d_check(const double* in, size_t order, 
		double* acoeff, double *err, double* kcoeff, double* tmp);

int flt_levinson1d(const float* in, size_t order, 
		float* acoeff, float *err, float* kcoeff, float* tmp);

int flt_levinson1d_check(const float* in, size_t order, 
		float* acoeff, float *err, float* kcoeff, float* tmp);

int dbl_levinson2d(const double* in, size_t dim0, size_t dim1,
		double* acoeff, double *err, double* kcoeff);

int flt_levinson2d(const float* in, size_t dim0, size_t dim1,
		float* acoeff, float *err, float* kcoeff);

#endif
