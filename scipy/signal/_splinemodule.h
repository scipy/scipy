#include "Python.h"
#include "numpy/arrayobject.h"
#include <math.h>


#define PYERR(message) do {PyErr_SetString(PyExc_ValueError, message); goto fail;} while(0)

static void convert_strides(npy_intp*,npy_intp*,int,int);

extern int S_SYM_IIR1_initial(float, float*, float*, int, int, float);
extern int S_SYM_IIR2_initial_fwd(double, double, float*, float*, int, int, float);
extern int S_SYM_IIR2_initial_bwd(double, double, float*, float*, int, int, float);
extern int S_separable_2Dconvolve_mirror(float*,float*,int,int,float*,float*,int,int,npy_intp*,npy_intp*);

extern int D_SYM_IIR1_initial(double, double*, double*, int, int, double);
extern int D_SYM_IIR2_initial_fwd(double, double, double*, double*, int, int, double);
extern int D_SYM_IIR2_initial_bwd(double, double, double*, double*, int, int, double);
extern int D_separable_2Dconvolve_mirror(double*,double*,int,int,double*,double*,int,int,npy_intp*,npy_intp*);

#ifdef __GNUC__
extern int C_SYM_IIR1_initial(__complex__ float, __complex__ float*, __complex__ float*, int, int, float);
extern int C_separable_2Dconvolve_mirror(__complex__ float*,__complex__ float*,int,int,__complex__ float*,__complex__ float*,int,int,npy_intp*,npy_intp*);
extern int Z_SYM_IIR1_initial(__complex__ double, __complex__ double*, __complex__ double*, int, int, double);
extern int Z_separable_2Dconvolve_mirror(__complex__ double*,__complex__ double*,int,int,__complex__ double*,__complex__ double*,int,int,npy_intp*,npy_intp*);
#endif
