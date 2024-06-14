#include "Python.h"
#include "numpy/arrayobject.h"
#include <math.h>


#define PYERR(message) do {PyErr_SetString(PyExc_ValueError, message); goto fail;} while(0)

static void convert_strides(npy_intp*,npy_intp*,int,int);

extern int S_cubic_spline2D(float*,float*,int,int,double,npy_intp*,npy_intp*,float);
extern int S_quadratic_spline2D(float*,float*,int,int,double,npy_intp*,npy_intp*,float);
extern int S_IIR_forback1(float,float,float*,float*,int,int,int,float);
extern int S_IIR_forback2(double,double,float*,float*,int,int,int,float);
extern int S_separable_2Dconvolve_mirror(float*,float*,int,int,float*,float*,int,int,npy_intp*,npy_intp*);

extern int D_cubic_spline2D(double*,double*,int,int,double,npy_intp*,npy_intp*,double);
extern int D_quadratic_spline2D(double*,double*,int,int,double,npy_intp*,npy_intp*,double);
extern int D_IIR_forback1(double,double,double*,double*,int,int,int,double);
extern int D_IIR_forback2(double,double,double*,double*,int,int,int,double);
extern int D_separable_2Dconvolve_mirror(double*,double*,int,int,double*,double*,int,int,npy_intp*,npy_intp*);

#ifdef __GNUC__
extern int C_IIR_forback1(__complex__ float,__complex__ float,__complex__ float*,__complex__ float*,int,int,int,float);
extern int C_separable_2Dconvolve_mirror(__complex__ float*,__complex__ float*,int,int,__complex__ float*,__complex__ float*,int,int,npy_intp*,npy_intp*);
extern int Z_IIR_forback1(__complex__ double,__complex__ double,__complex__ double*,__complex__ double*,int,int,int,double);
extern int Z_separable_2Dconvolve_mirror(__complex__ double*,__complex__ double*,int,int,__complex__ double*,__complex__ double*,int,int,npy_intp*,npy_intp*);
#endif
