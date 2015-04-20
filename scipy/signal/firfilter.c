#define NO_IMPORT_ARRAY
#include "sigtools.h"

static int elsizes[] = {sizeof(Bool),
			sizeof(byte),
                        sizeof(ubyte),
                        sizeof(short),
                        sizeof(ushort),
                        sizeof(int),
			sizeof(uint),
			sizeof(long),
                        sizeof(ulong),
                        sizeof(longlong),
			sizeof(ulonglong),
                        sizeof(float),
                        sizeof(double),
			sizeof(longdouble),
                        sizeof(cfloat),
                        sizeof(cdouble),
			sizeof(clongdouble),
                        sizeof(void *),
			0,0,0,0};

typedef void (OneMultAddFunction) (char *, char *, intp, char **, intp);

#define MAKE_ONEMULTADD(fname, type) \
static void fname ## _onemultadd(char *sum, char *term1, intp str, char **pvals, intp n) { \
        intp k; \
        type dsum = *(type*)sum; \
        for (k=0; k < n; k++) { \
          type tmp = *(type*)(term1 + k * str); \
          dsum += tmp * *(type*)pvals[k]; \
        } \
        *(type*)(sum) = dsum; \
}

MAKE_ONEMULTADD(UBYTE, ubyte)
MAKE_ONEMULTADD(USHORT, ushort)
MAKE_ONEMULTADD(UINT, uint)
MAKE_ONEMULTADD(ULONG, ulong)
MAKE_ONEMULTADD(ULONGLONG, ulonglong)

MAKE_ONEMULTADD(BYTE, byte)
MAKE_ONEMULTADD(SHORT, short)
MAKE_ONEMULTADD(INT, int)
MAKE_ONEMULTADD(LONG, long)
MAKE_ONEMULTADD(LONGLONG, longlong)

MAKE_ONEMULTADD(FLOAT, float)
MAKE_ONEMULTADD(DOUBLE, double)
MAKE_ONEMULTADD(LONGDOUBLE, longdouble)
 
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
static void fname ## _onemultadd(char *sum, char *term1, intp str, char **pvals, intp n) { \
        intp k; \
        for (k=0; k < n; k++) { \
          fname ## _onemultadd2(sum, term1 + k * str, pvals[k]); \
        } \
}
MAKE_C_ONEMULTADD(CFLOAT, float)
MAKE_C_ONEMULTADD(CDOUBLE, double)
MAKE_C_ONEMULTADD(CLONGDOUBLE, longdouble)
MAKE_C_ONEMULTADD2(CFLOAT, float)
MAKE_C_ONEMULTADD2(CDOUBLE, double)
MAKE_C_ONEMULTADD2(CLONGDOUBLE, longdouble)
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
		       intp   *instr,     /* Input strides */
		       char  *out,       /* Output data */
		       intp   *outstr,    /* Ouput strides */
		       char  *hvals,     /* coefficients in filter */
		       intp   *hstr,      /* coefficients strides */ 
		       intp   *Nwin,     /* Size of kernel Nwin[0] x Nwin[1] */
		       intp   *Ns,        /* Size of image Ns[0] x Ns[1] */
		       int   flag,       /* convolution parameters */
		       char  *fillvalue) /* fill value */
{
  int bounds_pad_flag = 0;
  int m, n, j, ind0, ind1;
  int Os[2];
  int new_m, new_n, ind0_memory=0;
  int boundary, outsize, convolve, type_num, type_size;
  char ** indices;
  OneMultAddFunction *mult_and_add;

  boundary = flag & BOUNDARY_MASK;  /* flag can be fill, reflecting, circular */
  outsize = flag & OUTSIZE_MASK;
  convolve = flag & FLIP_MASK;
  type_num = (flag & TYPE_MASK) >> TYPE_SHIFT;
  /*type_size*/

  mult_and_add = OneMultAdd[type_num];
  if (mult_and_add == NULL) return -5;  /* Not available for this type */

  if (type_num < 0 || type_num > MAXTYPES) return -4;  /* Invalid type */
  type_size = elsizes[type_num];

  indices = malloc(Nwin[1] * sizeof(indices[0]));
  if (indices == NULL) return -3; /* No memory */

  if (outsize == FULL) {Os[0] = Ns[0]+Nwin[0]-1; Os[1] = Ns[1]+Nwin[1]-1;}
  else if (outsize == SAME) {Os[0] = Ns[0]; Os[1] = Ns[1];}
  else if (outsize == VALID) {Os[0] = Ns[0]-Nwin[0]+1; Os[1] = Ns[1]-Nwin[1]+1;}
  else return -1;  /* Invalid output flag */  
  
  if ((boundary != PAD) && (boundary != REFLECT) && (boundary != CIRCULAR)) 
    return -2;   /* Invalid boundary flag */

  /* Speed this up by not doing any if statements in the for loop.  Need 3*3*2=18 different
     loops executed for different conditions */

  for (m=0; m < Os[0]; m++) {
    /* Reposition index into input image based on requested output size */
    if (outsize == FULL) new_m = convolve ? m : (m-Nwin[0]+1);
    else if (outsize == SAME) new_m = convolve ? (m+((Nwin[0]-1)>>1)) : (m-((Nwin[0]-1) >> 1));
    else new_m = convolve ? (m+Nwin[0]-1) : m; /* VALID */

    for (n=0; n < Os[1]; n++) {  /* loop over columns */
      char * sum = out+m*outstr[0]+n*outstr[1];
      memset(sum, 0, type_size); /* sum = 0.0; */

      if (outsize == FULL) new_n = convolve ? n : (n-Nwin[1]+1);
      else if (outsize == SAME) new_n = convolve ? (n+((Nwin[1]-1)>>1)) : (n-((Nwin[1]-1) >> 1));
      else new_n = convolve ? (n+Nwin[1]-1) : n;

      /* Sum over kernel, if index into image is out of bounds
	 handle it according to boundary flag */
      for (j=0; j < Nwin[0]; j++) {
	ind0 = convolve ? (new_m-j): (new_m+j);
	bounds_pad_flag = 0;

	if (ind0 < 0) {
	  if (boundary == REFLECT) ind0 = -1-ind0;
	  else if (boundary == CIRCULAR) ind0 = Ns[0] + ind0;
	  else bounds_pad_flag = 1;
	}
	else if (ind0 >= Ns[0]) {
	  if (boundary == REFLECT) ind0 = Ns[0]+Ns[0]-1-ind0;
	  else if (boundary == CIRCULAR) ind0 = ind0 - Ns[0];
	  else bounds_pad_flag = 1;
	}
	
	if (!bounds_pad_flag) ind0_memory = ind0*instr[0];

        if (bounds_pad_flag) {
          intp k;
          for (k=0; k < Nwin[1]; k++) {
              indices[k] = fillvalue;
          }
        }
        else  {
          intp k;
	  for (k=0; k < Nwin[1]; k++) {
	    ind1 = convolve ? (new_n-k) : (new_n+k);
	    if (ind1 < 0) {
	      if (boundary == REFLECT) ind1 = -1-ind1;
	      else if (boundary == CIRCULAR) ind1 = Ns[1] + ind1;
	      else bounds_pad_flag = 1;
	    }
	    else if (ind1 >= Ns[1]) {
	      if (boundary == REFLECT) ind1 = Ns[1]+Ns[1]-1-ind1;
	      else if (boundary == CIRCULAR) ind1 = ind1 - Ns[1];
	      else bounds_pad_flag = 1;
	    }

	    if (bounds_pad_flag) {
              indices[k] = fillvalue;
            }
	    else {
              indices[k] = in+ind0_memory+ind1*instr[1];
            }
	    bounds_pad_flag = 0;
	  }
        }
        mult_and_add(sum, hvals+j*hstr[0], hstr[1], indices, Nwin[1]);
      }
    }
  }
  free(indices);
  return 0;
}
