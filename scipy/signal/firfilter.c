#define NO_IMPORT_ARRAY
#include "numpy/ndarrayobject.h"
#include "sigtools.h"
#include <stdbool.h>
#include <stdint.h>

static int elsizes[] = {sizeof(npy_bool),
                        sizeof(npy_byte),
                        sizeof(npy_ubyte),
                        sizeof(npy_short),
                        sizeof(npy_ushort),
                        sizeof(int),
                        sizeof(npy_uint),
                        sizeof(long),
                        sizeof(npy_ulong),
                        sizeof(npy_longlong),
                        sizeof(npy_ulonglong),
                        sizeof(float),
                        sizeof(double),
                        sizeof(npy_longdouble),
                        sizeof(npy_cfloat),
                        sizeof(npy_cdouble),
                        sizeof(npy_clongdouble),
                        sizeof(void *),
			0,0,0,0};

typedef void (OneMultAddFunction) (char *, char *, int64_t, char **, int64_t);

#define MAKE_ONEMULTADD(fname, type) \
static void fname ## _onemultadd(char *sum, char *term1, int64_t str, char **pvals, int64_t n) { \
        type dsum = *(type*)sum; \
        for (int64_t k=0; k < n; k++) { \
          type tmp = *(type*)(term1 + k * str); \
          dsum += tmp * *(type*)pvals[k]; \
        } \
        *(type*)(sum) = dsum; \
}

MAKE_ONEMULTADD(UBYTE, npy_ubyte)
MAKE_ONEMULTADD(USHORT, npy_ushort)
MAKE_ONEMULTADD(UINT, npy_uint)
MAKE_ONEMULTADD(ULONG, npy_ulong)
MAKE_ONEMULTADD(ULONGLONG, npy_ulonglong)

MAKE_ONEMULTADD(BYTE, npy_byte)
MAKE_ONEMULTADD(SHORT, short)
MAKE_ONEMULTADD(INT, int)
MAKE_ONEMULTADD(LONG, long)
MAKE_ONEMULTADD(LONGLONG, npy_longlong)

MAKE_ONEMULTADD(FLOAT, float)
MAKE_ONEMULTADD(DOUBLE, double)
MAKE_ONEMULTADD(LONGDOUBLE, npy_longdouble)
 
#ifdef __GNUC__
MAKE_ONEMULTADD(CFLOAT, __complex__ float)
MAKE_ONEMULTADD(CDOUBLE, __complex__ double)
MAKE_ONEMULTADD(CLONGDOUBLE, __complex__ long double)
#else
#define MAKE_C_ONEMULTADD(fname, type) \
static void fname ## _onemultadd2(char *sum, char *term1, char *term2) { \
  ((type *) sum)[0] += ((type *) term1)[0] * ((type *) term2)[0] \
    - ((type *) term1)[1] * ((type *) term2)[1]; \
  ((type *) sum)[1] += ((type *) term1)[0] * ((type *) term2)[1] \
    + ((type *) term1)[1] * ((type *) term2)[0]; \
  return; }

#define MAKE_C_ONEMULTADD2(fname, type) \
static void fname ## _onemultadd(char *sum, char *term1, int64_t str, \
                                 char **pvals, int64_t n) { \
        for (int64_t k=0; k < n; k++) { \
          fname ## _onemultadd2(sum, term1 + k * str, pvals[k]); \
        } \
}
MAKE_C_ONEMULTADD(CFLOAT, float)
MAKE_C_ONEMULTADD(CDOUBLE, double)
MAKE_C_ONEMULTADD(CLONGDOUBLE, npy_longdouble)
MAKE_C_ONEMULTADD2(CFLOAT, float)
MAKE_C_ONEMULTADD2(CDOUBLE, double)
MAKE_C_ONEMULTADD2(CLONGDOUBLE, npy_longdouble)
#endif /* __GNUC__ */

static OneMultAddFunction *OneMultAdd[]={NULL,
					 BYTE_onemultadd,
					 UBYTE_onemultadd,
					 SHORT_onemultadd,
                                         USHORT_onemultadd,
					 INT_onemultadd,
                                         UINT_onemultadd,
					 LONG_onemultadd,
					 ULONG_onemultadd,
					 LONGLONG_onemultadd,
					 ULONGLONG_onemultadd,
					 FLOAT_onemultadd,
					 DOUBLE_onemultadd,
					 LONGDOUBLE_onemultadd,
					 CFLOAT_onemultadd,
					 CDOUBLE_onemultadd,
					 CLONGDOUBLE_onemultadd,
                                         NULL, NULL, NULL, NULL};


/* This could definitely be more optimized... */


int pylab_convolve_2d (char  *in,        /* Input data Ns[0] x Ns[1] */
		       npy_intp   *instr,     /* Input strides */
		       char  *out,       /* Output data */
		       npy_intp   *outstr,    /* Output strides */
		       char  *hvals,     /* coefficients in filter */
		       npy_intp   *hstr,      /* coefficients strides */ 
		       npy_intp   *Nwin,     /* Size of kernel Nwin[0] x Nwin[1] */
		       npy_intp   *Ns,        /* Size of image Ns[0] x Ns[1] */
		       int   flag,       /* convolution parameters */
		       char  *fillvalue) /* fill value */
{
  const int boundary = flag & BOUNDARY_MASK;  /* flag can be fill, reflecting, circular */
  const int outsize = flag & OUTSIZE_MASK;
  const int convolve = flag & FLIP_MASK;
  const int type_num = (flag & TYPE_MASK) >> TYPE_SHIFT;
  /*type_size*/

  OneMultAddFunction *mult_and_add = OneMultAdd[type_num];
  if (mult_and_add == NULL) return -5;  /* Not available for this type */

  if (type_num < 0 || type_num > MAXTYPES) return -4;  /* Invalid type */
  const int type_size = elsizes[type_num];

  int64_t Os[2];
  if (outsize == FULL) {Os[0] = Ns[0]+Nwin[0]-1; Os[1] = Ns[1]+Nwin[1]-1;}
  else if (outsize == SAME) {Os[0] = Ns[0]; Os[1] = Ns[1];}
  else if (outsize == VALID) {Os[0] = Ns[0]-Nwin[0]+1; Os[1] = Ns[1]-Nwin[1]+1;}
  else return -1; /* Invalid output flag */

  if ((boundary != PAD) && (boundary != REFLECT) && (boundary != CIRCULAR))
    return -2; /* Invalid boundary flag */

  char **indices = malloc(Nwin[1] * sizeof(indices[0]));
  if (indices == NULL) return -3; /* No memory */

  /* Speed this up by not doing any if statements in the for loop.  Need 3*3*2=18 different
     loops executed for different conditions */

  for (int64_t m=0; m < Os[0]; m++) {
    /* Reposition index into input image based on requested output size */
    int64_t new_m;
    if (outsize == FULL) new_m = convolve ? m : (m-Nwin[0]+1);
    else if (outsize == SAME) new_m = convolve ? (m+((Nwin[0]-1)>>1)) : (m-((Nwin[0]-1) >> 1));
    else new_m = convolve ? (m+Nwin[0]-1) : m; /* VALID */

    for (int64_t n=0; n < Os[1]; n++) {  /* loop over columns */
      char * sum = out+m*outstr[0]+n*outstr[1];
      memset(sum, 0, type_size); /* sum = 0.0; */

      int64_t new_n;
      if (outsize == FULL) new_n = convolve ? n : (n-Nwin[1]+1);
      else if (outsize == SAME) new_n = convolve ? (n+((Nwin[1]-1)>>1)) : (n-((Nwin[1]-1) >> 1));
      else new_n = convolve ? (n+Nwin[1]-1) : n;

      /* Sum over kernel, if index into image is out of bounds
	 handle it according to boundary flag */
      for (int64_t j=0; j < Nwin[0]; j++) {
	int64_t ind0 = convolve ? (new_m-j): (new_m+j);
	bool bounds_pad_flag = false;

	if (ind0 < 0) {
	  if (boundary == REFLECT) ind0 = -1-ind0;
	  else if (boundary == CIRCULAR) ind0 = Ns[0] + ind0;
	  else bounds_pad_flag = true;
	}
	else if (ind0 >= Ns[0]) {
	  if (boundary == REFLECT) ind0 = Ns[0]+Ns[0]-1-ind0;
	  else if (boundary == CIRCULAR) ind0 = ind0 - Ns[0];
	  else bounds_pad_flag = true;
	}

	const int64_t ind0_memory = ind0*instr[0];

	if (bounds_pad_flag) {
	  for (int64_t k=0; k < Nwin[1]; k++) {
	      indices[k] = fillvalue;
	  }
	}
	else  {
	  for (int64_t k=0; k < Nwin[1]; k++) {
	    int64_t ind1 = convolve ? (new_n-k) : (new_n+k);
	    if (ind1 < 0) {
	      if (boundary == REFLECT) ind1 = -1-ind1;
	      else if (boundary == CIRCULAR) ind1 = Ns[1] + ind1;
	      else bounds_pad_flag = true;
	    }
	    else if (ind1 >= Ns[1]) {
	      if (boundary == REFLECT) ind1 = Ns[1]+Ns[1]-1-ind1;
	      else if (boundary == CIRCULAR) ind1 = ind1 - Ns[1];
	      else bounds_pad_flag = true;
	    }

	    if (bounds_pad_flag) {
	      indices[k] = fillvalue;
	    }
	    else {
	      indices[k] = in+ind0_memory+ind1*instr[1];
	    }
	    bounds_pad_flag = false;
	  }
	}
	mult_and_add(sum, hvals+j*hstr[0], hstr[1], indices, Nwin[1]);
      }
    }
  }
  free(indices);
  return 0;
}
