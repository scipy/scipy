
/*--------------------------------------------------------------------*/

#include "Python.h"
#define NO_IMPORT_ARRAY
#include "numpy/ndarrayobject.h"


/* defined below */
void f_medfilt2(float*,float*,npy_intp*,npy_intp*);
void d_medfilt2(double*,double*,npy_intp*,npy_intp*);
void b_medfilt2(unsigned char*,unsigned char*,npy_intp*,npy_intp*);
extern char *check_malloc (int);


/* 2-D median filter with zero-padding on edges. */
#define MEDIAN_FILTER_2D(NAME, TYPE, SELECT)                            \
void NAME(TYPE* in, TYPE* out, npy_intp* Nwin, npy_intp* Ns)                    \
{                                                                       \
    int nx, ny, hN[2];                                                  \
    int pre_x, pre_y, pos_x, pos_y;                                     \
    int subx, suby, k, totN;                                            \
    TYPE *myvals, *fptr1, *fptr2, *ptr1, *ptr2;                         \
                                                                        \
    totN = Nwin[0] * Nwin[1];                                           \
    myvals = (TYPE *) check_malloc( totN * sizeof(TYPE));               \
                                                                        \
    Py_BEGIN_ALLOW_THREADS                                              \
                                                                        \
    hN[0] = Nwin[0] >> 1;                                               \
    hN[1] = Nwin[1] >> 1;                                               \
    ptr1 = in;                                                          \
    fptr1 = out;                                                        \
    for (ny = 0; ny < Ns[0]; ny++)                                      \
        for (nx = 0; nx < Ns[1]; nx++) {                                \
            pre_x = hN[1];                                              \
            pre_y = hN[0];                                              \
            pos_x = hN[1];                                              \
            pos_y = hN[0];                                              \
            if (nx < hN[1]) pre_x = nx;                                 \
            if (nx >= Ns[1] - hN[1]) pos_x = Ns[1] - nx - 1;            \
            if (ny < hN[0]) pre_y = ny;                                 \
            if (ny >= Ns[0] - hN[0]) pos_y = Ns[0] - ny - 1;            \
            fptr2 = myvals;                                             \
            ptr2 = ptr1 - pre_x - pre_y*Ns[1];                          \
            for (suby = -pre_y; suby <= pos_y; suby++) {                \
                for (subx = -pre_x; subx <= pos_x; subx++)              \
                    *fptr2++ = *ptr2++;                                 \
                ptr2 += Ns[1] - (pre_x + pos_x + 1);                    \
            }                                                           \
            ptr1++;                                                     \
                                                                        \
            /* Zero pad */                                              \
            for (k = (pre_x + pos_x + 1)*(pre_y + pos_y + 1); k < totN; k++) \
                *fptr2++ = 0.0;                                         \
                                                                        \
            /*      *fptr1++ = median(myvals,totN); */                  \
            *fptr1++ = SELECT(myvals,totN);                             \
        }                                                               \
                                                                        \
    Py_END_ALLOW_THREADS                                                \
                                                                        \
    free(myvals);                                                       \
}


/* Stolen from CYTHON_RESTRICT and BOOST_RESTRICT */
#ifndef SCIPY_RESTRICT
#  if defined(__GNUC__) && __GNUC__ > 3
#    define SCIPY_RESTRICT __restrict__
#  elif defined(_MSC_VER) && _MSC_VER >= 1400
#    define SCIPY_RESTRICT __restrict
#  elif defined(__STDC_VERSION__) && __STDC_VERSION__ >= 199901L
#    define SCIPY_RESTRICT restrict
#  else
#    define SCIPY_RESTRICT
#  endif
#endif

/* define quick_select for floats, doubles, and unsigned characters */
typedef int (*CompareFunction)(const void * const, const void * const);
extern void *quick_select(void * SCIPY_RESTRICT base, const size_t num_elements,
                   const size_t element_size,
		   const CompareFunction comparison_function,
		   const size_t element_to_return);
extern int FLOAT_compare(const float *const ip1, const float *const ip2);
extern int DOUBLE_compare(const double *const ip1, const double *const ip2);
extern int UBYTE_compare(const npy_byte *const ip1, const npy_byte *const ip2);

/* index = (n - 1) / 2; use lower of middle values for even-length arrays */
float f_quick_select(float arr[], const int n) {
  float *result = quick_select(arr, n, sizeof(float), FLOAT_compare, (n - 1) / 2);
  return *result;
}
double d_quick_select(double arr[], const int n) {
  double *result = quick_select(arr, n, sizeof(double), DOUBLE_compare, (n - 1) / 2);
  return *result;
}
unsigned char b_quick_select(unsigned char arr[], const int n) {
  unsigned char *result = quick_select(arr, n, sizeof(unsigned char), UBYTE_compare, (n - 1) / 2);
  return *result;
}

/* define medfilt for floats, doubles, and unsigned characters */
MEDIAN_FILTER_2D(f_medfilt2, float, f_quick_select)
MEDIAN_FILTER_2D(d_medfilt2, double, d_quick_select)
MEDIAN_FILTER_2D(b_medfilt2, unsigned char, b_quick_select)
