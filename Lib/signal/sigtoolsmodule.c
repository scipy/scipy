/* SIGTOOLS module by Travis Oliphant

Copyright 2001 Travis Oliphant
Permission to use, copy, modify, and distribute this software without fee
is granted, under the terms of the LGPL provided this notification remain.
*/

#include "Python.h"
#include "Numeric/arrayobject.h"
#include <setjmp.h>
#include "sigtools.h"

#define PYERR(message) {PyErr_SetString(PyExc_ValueError, message); goto fail;}

#define DATA(arr) ((arr)->data)
#define DIMS(arr) ((arr)->dimensions)
#define STRIDES(arr) ((arr)->strides)
#define ELSIZE(arr) ((arr)->descr->elsize)
#define OBJECTTYPE(arr) ((arr)->descr->type_num)
#define BASEOBJ(arr) ((PyArrayObject *)((arr)->base))
#define RANK(arr) ((arr)->nd)
#define ISCONTIGUOUS(m) ((m)->flags & CONTIGUOUS)


jmp_buf MALLOC_FAIL;

char *check_malloc (int);

char *check_malloc (size)
	int size;
{
    char *the_block;
    
    the_block = (char *)malloc(size);
    if (the_block == NULL)
	{
	    printf("\nERROR: unable to allocate %d bytes!\n", size);
	    longjmp(MALLOC_FAIL,-1);
	}
    return(the_block);
}



/************************************************************************
 * Start of portable, non-python specific routines.                     *
 ************************************************************************/

/* Some core routines are written
in a portable way so that they could be used in other applications.  The 
order filtering, however uses python-specific constructs in its guts 
and is therefore Python dependent.  This could be changed in a 
straightforward way but I haven't done it for lack of time.*/

static int index_out_of_bounds(int *indices, int *max_indices, int ndims) {
  int bad_index = 0, k = 0;

  while (!bad_index && (k++ < ndims)) {
    bad_index = ((*(indices) >= *(max_indices++)) || (*(indices) < 0));
    indices++;
  }
  return bad_index;
}

/* This maybe could be redone with stride information so it could be 
 * called with non-contiguous arrays:  I think offsets is related to 
 * the difference between the strides.  I'm not sure about init_offset 
 * just yet.  I think it needs to be calculated because of mode_dep
 * but probably with dim1 being the size of the "original, unsliced" array
 */

static long compute_offsets (unsigned long *offsets, long *offsets2, int *dim1, int *dim2, int *dim3, int *mode_dep, int nd) {
  int k,i;
  long init_offset = 0;

  for (k = 0; k < nd - 1; k++) 
    {
      init_offset += mode_dep[k];
      init_offset *= dim1[k+1];
    }
  init_offset += mode_dep[k] - 2;
  
  k = nd;
  while(k--) {
    offsets[k] = 0;
    offsets2[k] = 0;
    for (i = k + 1; i < nd - 1; i++) {
      offsets[k] += dim1[i] - dim2[i];
      offsets[k] *= dim1[i+1];

      offsets2[k] += dim1[i] - dim3[i];
      offsets2[k] *= dim1[i+1];
    }

    if (k < nd - 1) {
      offsets[k] += dim1[i] - dim2[i];
      offsets2[k] += dim1[i] - dim3[i];
    }
    offsets[k] += 1;
    offsets2[k] += 1;
  }
  return init_offset;
}

/* increment by 1 the index into an N-D array, doing the necessary
   carrying when the index reaches the dimension along that axis */ 
static int increment(int *ret_ind, int nd, int *max_ind) {    
    int k, incr = 1;
    
    k = nd - 1;
    if (++ret_ind[k] >= max_ind[k]) {
      while (k >= 0 && (ret_ind[k] >= max_ind[k]-1)) {
	incr++;
	ret_ind[k--] = 0;
      }
      if (k >= 0) ret_ind[k]++;
    }
    return incr;
}

/*
 All of these MultAdd functions loop over all the elements of the smallest
 array, incrementing an array of indices into the large N-D array at
 the same time.  The
 bounds for the other array are checked and if valid the product is
 added to the running sum.  If invalid bounds are found nothing is
 done (zero is added).  This has the effect of zero padding the array
 to handle edges.
 */

static void UBYTE_MultAdd(char *ip1, int is1, char *ip2, int is2, char *op, int *dims1, int *dims2, int ndims, int nels2, int check, int *loop_ind, int *temp_ind, unsigned long *offset) { 
  unsigned char tmp=(unsigned char)0.0;  int i, k, incr = 1;
  unsigned char *ptr1 = (unsigned char *)ip1, *ptr2 = (unsigned char *)ip2;

  i = nels2;

  temp_ind[ndims-1]--;
  while (i--) { 
    /* Adjust index array and move ptr1 to right place */
    k = ndims - 1;
    while(--incr) {
      temp_ind[k] -= dims2[k] - 1;   /* Return to start for these dimensions */
      k--;
    }
    ptr1 += offset[k];               /* Precomputed offset array */
    temp_ind[k]++;

    if (!(check && index_out_of_bounds(temp_ind,dims1,ndims))) { 
      tmp += (*ptr1) * (*ptr2); 
    } 
    incr = increment(loop_ind, ndims, dims2);  /* Returns number of N-D indices incremented. */
    ptr2++;

  }
  *((unsigned char *)op) = tmp; 
}

static void SBYTE_MultAdd(char *ip1, int is1, char *ip2, int is2, char *op, int *dims1, int *dims2, int ndims, int nels2, int check, int *loop_ind, int *temp_ind, unsigned long *offset) { 
  signed char tmp=(signed char)0.0;  int i, k, incr = 1;
  signed char *ptr1 = (signed char *)ip1, *ptr2 = (signed char *)ip2;

  i = nels2;

  temp_ind[ndims-1]--;
  while (i--) { 
    /* Adjust index array and move ptr1 to right place */
    k = ndims - 1;
    while(--incr) {
      temp_ind[k] -= dims2[k] - 1;   /* Return to start for these dimensions */
      k--;
    }
    ptr1 += offset[k];               /* Precomputed offset array */
    temp_ind[k]++;

    if (!(check && index_out_of_bounds(temp_ind,dims1,ndims))) { 
      tmp += (*ptr1) * (*ptr2); 
    } 
    incr = increment(loop_ind, ndims, dims2);  /* Returns number of N-D indices incremented. */
    ptr2++;

  }
  *((signed char *)op) = tmp; 
}

static void SHORT_MultAdd(char *ip1, int is1, char *ip2, int is2, char *op, int *dims1, int *dims2, int ndims, int nels2, int check, int *loop_ind, int *temp_ind, unsigned long *offset) { 
  short tmp=(short)0.0;  int i, k, incr = 1;
  short *ptr1 = (short *)ip1, *ptr2 = (short *)ip2;

  i = nels2;

  temp_ind[ndims-1]--;
  while (i--) { 
    /* Adjust index array and move ptr1 to right place */
    k = ndims - 1;
    while(--incr) {
      temp_ind[k] -= dims2[k] - 1;   /* Return to start for these dimensions */
      k--;
    }
    ptr1 += offset[k];               /* Precomputed offset array */
    temp_ind[k]++;

    if (!(check && index_out_of_bounds(temp_ind,dims1,ndims))) { 
      tmp += (*ptr1) * (*ptr2); 
    } 
    incr = increment(loop_ind, ndims, dims2);  /* Returns number of N-D indices incremented. */
    ptr2++;

  }
  *((short *)op) = tmp; 
}

static void INT_MultAdd(char *ip1, int is1, char *ip2, int is2, char *op, int *dims1, int *dims2, int ndims, int nels2, int check, int *loop_ind, int *temp_ind, unsigned long *offset) { 
  int tmp=(int)0.0;  int i, k, incr = 1;
  int *ptr1 = (int *)ip1, *ptr2 = (int *)ip2;

  i = nels2;

  temp_ind[ndims-1]--;
  while (i--) { 
    /* Adjust index array and move ptr1 to right place */
    k = ndims - 1;
    while(--incr) {
      temp_ind[k] -= dims2[k] - 1;   /* Return to start for these dimensions */
      k--;
    }
    ptr1 += offset[k];               /* Precomputed offset array */
    temp_ind[k]++;

    if (!(check && index_out_of_bounds(temp_ind,dims1,ndims))) { 
      tmp += (*ptr1) * (*ptr2); 
    } 
    incr = increment(loop_ind, ndims, dims2);  /* Returns number of N-D indices incremented. */
    ptr2++;

  }
  *((int *)op) = tmp; 
}

static void LONG_MultAdd(char *ip1, int is1, char *ip2, int is2, char *op, int *dims1, int *dims2, int ndims, int nels2, int check, int *loop_ind, int *temp_ind, unsigned long *offset) { 
  long tmp=(long)0.0;  int i, k, incr = 1;
  long *ptr1 = (long *)ip1, *ptr2 = (long *)ip2;

  i = nels2;

  temp_ind[ndims-1]--;
  while (i--) { 
    /* Adjust index array and move ptr1 to right place */
    k = ndims - 1;
    while(--incr) {
      temp_ind[k] -= dims2[k] - 1;   /* Return to start for these dimensions */
      k--;
    }
    ptr1 += offset[k];               /* Precomputed offset array */
    temp_ind[k]++;

    if (!(check && index_out_of_bounds(temp_ind,dims1,ndims))) { 
      tmp += (*ptr1) * (*ptr2); 
    } 
    incr = increment(loop_ind, ndims, dims2);  /* Returns number of N-D indices incremented. */
    ptr2++;

  }
  *((long *)op) = tmp; 
}

static void FLOAT_MultAdd(char *ip1, int is1, char *ip2, int is2, char *op, int *dims1, int *dims2, int ndims, int nels2, int check, int *loop_ind, int *temp_ind, unsigned long *offset) { 
  float tmp=(float)0.0;  int i, k, incr = 1;
  float *ptr1 = (float *)ip1, *ptr2 = (float *)ip2;

  i = nels2;

  temp_ind[ndims-1]--;
  while (i--) { 
    /* Adjust index array and move ptr1 to right place */
    k = ndims - 1;
    while(--incr) {
      temp_ind[k] -= dims2[k] - 1;   /* Return to start for these dimensions */
      k--;
    }
    ptr1 += offset[k];               /* Precomputed offset array */
    temp_ind[k]++;

    if (!(check && index_out_of_bounds(temp_ind,dims1,ndims))) { 
      tmp += (*ptr1) * (*ptr2); 
    } 
    incr = increment(loop_ind, ndims, dims2);  /* Returns number of N-D indices incremented. */
    ptr2++;

  }
  *((float *)op) = tmp; 
}

static void DOUBLE_MultAdd(char *ip1, int is1, char *ip2, int is2, char *op, int *dims1, int *dims2, int ndims, int nels2, int check, int *loop_ind, int *temp_ind, unsigned long *offset) { 
  double tmp=(double)0.0;  int i, k, incr = 1;
  double *ptr1 = (double *)ip1, *ptr2 = (double *)ip2;

  i = nels2;

  temp_ind[ndims-1]--;
  while (i--) { 
    /* Adjust index array and move ptr1 to right place */
    k = ndims - 1;
    while(--incr) {
      temp_ind[k] -= dims2[k] - 1;   /* Return to start for these dimensions */
      k--;
    }
    ptr1 += offset[k];               /* Precomputed offset array */
    temp_ind[k]++;

    if (!(check && index_out_of_bounds(temp_ind,dims1,ndims))) { 
      tmp += (*ptr1) * (*ptr2); 
    } 
    incr = increment(loop_ind, ndims, dims2);  /* Returns number of N-D indices incremented. */
    ptr2++;

  }
  *((double *)op) = tmp; 
}

static void CFLOAT_MultAdd(char *ip1, int is1, char *ip2, int is2, char *op, int *dims1, int *dims2, int ndims, int nels2, int check, int *loop_ind, int *temp_ind, unsigned long *offset) { 
  float tmpr= 0.0, tmpi = 0.0; int i, k, incr = 1;
  float *ptr1 = (float *)ip1, *ptr2 = (float *)ip2;

  i = nels2;

  temp_ind[ndims-1]--;
  while (i--) { 
    /* Adjust index array and move ptr1 to right place */
    k = ndims - 1;
    while(--incr) {
      temp_ind[k] -= dims2[k] - 1;   /* Return to start for these dimensions */
      k--;
    }
    ptr1 += 2*offset[k];               /* Precomputed offset array */
    temp_ind[k]++;

    if (!(check && index_out_of_bounds(temp_ind,dims1,ndims))) { 
      tmpr += ptr1[0] * ptr2[0] - ptr1[1] * ptr2[1]; 
      tmpi += ptr1[1] * ptr2[0] + ptr1[0] * ptr2[1]; 
    } 
    incr = increment(loop_ind, ndims, dims2);  /* Returns number of N-D indices incremented. */
    ptr2 += 2;

  }
  ((float *)op)[0] = tmpr; ((float *)op)[1] = tmpi; 
}

static void CDOUBLE_MultAdd(char *ip1, int is1, char *ip2, int is2, char *op, int *dims1, int *dims2, int ndims, int nels2, int check, int *loop_ind, int *temp_ind, unsigned long *offset) { 
  double tmpr= 0.0, tmpi = 0.0; int i, k, incr = 1;
  double *ptr1 = (double *)ip1, *ptr2 = (double *)ip2;

  i = nels2;

  temp_ind[ndims-1]--;
  while (i--) { 
    /* Adjust index array and move ptr1 to right place */
    k = ndims - 1;
    while(--incr) {
      temp_ind[k] -= dims2[k] - 1;   /* Return to start for these dimensions */
      k--;
    }
    ptr1 += 2*offset[k];               /* Precomputed offset array */
    temp_ind[k]++;

    if (!(check && index_out_of_bounds(temp_ind,dims1,ndims))) { 
      tmpr += ptr1[0] * ptr2[0] - ptr1[1] * ptr2[1]; 
      tmpi += ptr1[1] * ptr2[0] + ptr1[0] * ptr2[1]; 
    } 
    incr = increment(loop_ind, ndims, dims2);  /* Returns number of N-D indices incremented. */
    ptr2 += 2;

  }
  ((double *)op)[0] = tmpr; ((double *)op)[1] = tmpi; 
}


static void correlateND(Generic_Array *ap1, Generic_Array *ap2, Generic_Array *ret, MultAddFunction *multiply_and_add_ND, int mode) {
	int *a_ind, *b_ind, *temp_ind, *mode_dep, *check_ind;
	unsigned long *offsets, offset1;
	long *offsets2;
	int i, k, check, incr = 1;
	int bytes_in_array, num_els_ret, num_els_ap2;
	int is1, is2, os;
	char *ip1, *ip2, *op, *ap1_ptr;
	int *ret_ind;

        num_els_ret = 1;
	for (i = 0; i < ret->nd; i++) num_els_ret *= ret->dimensions[i];
	num_els_ap2 = 1;
	for (i = 0; i < ret->nd; i++) num_els_ap2 *= ap2->dimensions[i];
	bytes_in_array = ap1->nd * sizeof(int);
	mode_dep = (int *)malloc(bytes_in_array);
	switch(mode) {
	case 0:
	        for (i = 0; i < ap1->nd; i++) mode_dep[i] = 0;
		break;
	case 1:
	        for (i = 0; i < ap1->nd; i++) mode_dep[i] = -((ap2->dimensions[i]) >> 1);
		break;
	case 2:
	        for (i = 0; i < ap1->nd; i++) mode_dep[i] = 1 - ap2->dimensions[i];
	}

	is1 = ap1->elsize; is2 = ap2->elsize;
	op = ret->data; os = ret->elsize;
	ip1 = ap1->data; ip2 = ap2->data;
	op = ret->data;

	b_ind = (int *)malloc(bytes_in_array);  /* loop variables */
	memset(b_ind,0,bytes_in_array);
	a_ind = (int *)malloc(bytes_in_array);
	ret_ind = (int *)malloc(bytes_in_array);
	memset(ret_ind,0,bytes_in_array);
	temp_ind = (int *)malloc(bytes_in_array);
	check_ind = (int*)malloc(bytes_in_array);
	offsets = (unsigned long *)malloc(ap1->nd*sizeof(unsigned long));
	offsets2 = (long *)malloc(ap1->nd*sizeof(long));
	offset1 = compute_offsets(offsets,offsets2,ap1->dimensions,ap2->dimensions,ret->dimensions,mode_dep,ap1->nd);
	/* The convolution proceeds by looping through the output array
	   and for each value summing all contributions from the summed
	   element-by-element product of the two input arrays.  Index
	   counters are used for book-keeping in the area so that we 
	   can tell where we are in all of the arrays and be sure that 
	   we are not trying to access areas outside the arrays definition.

	   The inner loop is implemented separately but equivalently for each
	   datatype. The outer loop is similar in structure and form to
	   to the inner loop.
	*/
	/* Need to keep track of a ptr to place in big (first) input
	   array where we start the multiplication (we pass over it in the
	   inner loop (and not dereferenced) 
	   if it is pointing outside dataspace)
	*/
	/* Calculate it once and the just move it around appropriately */
	ap1_ptr = ip1 + offset1*is1;
	for (k=0; k < ap1->nd; k++) {a_ind[k] = mode_dep[k]; check_ind[k] = ap1->dimensions[k] - ap2->dimensions[k] - mode_dep[k] - 1;}
	a_ind[ap1->nd-1]--;
	i = num_els_ret;
	while (i--) {
	  k = ap1->nd - 1;
	  while(--incr) {
	    a_ind[k] -= ret->dimensions[k] - 1;   /* Return to start */
	    k--;
	  }
	  ap1_ptr += offsets2[k]*is1;
	  a_ind[k]++;
	  memcpy(temp_ind, a_ind, bytes_in_array);

	  check = 0; k = -1;
	  while(!check && (++k < ap1->nd))
	    check = check || (ret_ind[k] < -mode_dep[k]) || (ret_ind[k] > check_ind[k]);

	  multiply_and_add_ND(ap1_ptr,is1,ip2,is2,op,ap1->dimensions,ap2->dimensions,ap1->nd,num_els_ap2,check,b_ind,temp_ind,offsets);

	  incr = increment(ret_ind,ret->nd,ret->dimensions); /* increment index counter */
	  op += os;   /* increment to next output index */

	}
	free(b_ind); free(a_ind); free(ret_ind);
	free(offsets); free(offsets2); free(temp_ind);
	free(check_ind); free(mode_dep);
}

/*****************************************************************
 *   This is code for a 1-D linear-filter along an arbitrary     *
 *   dimension of an N-D array.                                  *
 *****************************************************************/

static void FLOAT_filt(char *b, char *a, char *x, char *y, char *Z, int len_b, unsigned int len_x, int stride_X, int stride_Y ) {
  char *ptr_x = x, *ptr_y = y;
  float *ptr_Z, *ptr_b;
  float *ptr_a;
  float *xn, *yn;
  const float a0 = *((float *)a);
  int k, n;

  for (k = 0; k < len_x; k++) {
    ptr_b = (float *)b;          /* Reset a and b pointers */
    ptr_a = (float *)a;
    xn = (float *)ptr_x;
    yn = (float *)ptr_y;
    if (len_b > 1) {
      ptr_Z = ((float *)Z);
      *yn = *ptr_Z + *ptr_b / a0 * *xn;    /* Calculate first delay (output) */
      ptr_b++; ptr_a++;
    /* Fill in middle delays */
      for (n = 0; n < len_b - 2; n++) {
	*ptr_Z = ptr_Z[1] + *xn * (*ptr_b / a0) - *yn * (*ptr_a / a0);
	ptr_b++; ptr_a++; ptr_Z++;
      }
    /* Calculate last delay */
      *ptr_Z = *xn * (*ptr_b / a0) - *yn * (*ptr_a / a0);
    }
    else {
      *yn = *xn * (*ptr_b / a0);
    }

    ptr_y += stride_Y;             /* Move to next input/output point */
    ptr_x += stride_X;     
  }
}


static void DOUBLE_filt(char *b, char *a, char *x, char *y, char *Z, int len_b, unsigned int len_x, int stride_X, int stride_Y ) {
  char *ptr_x = x, *ptr_y = y;
  double *ptr_Z, *ptr_b;
  double *ptr_a;
  double *xn, *yn;
  double a0;
  int k, n;

  a0 = *((double *)a);
  for (k = 0; k < len_x; k++) {
    ptr_b = (double *)b;          /* Reset a and b pointers */
    ptr_a = (double *)a;
    xn = (double *)ptr_x;
    yn = (double *)ptr_y;
    if (len_b > 1) {
      ptr_Z = ((double *)Z);
      *yn = *ptr_Z + *ptr_b / a0 * *xn;    /* Calculate first delay (output) */
      ptr_b++; ptr_a++;
      /* Fill in middle delays */
      for (n = 0; n < len_b - 2; n++) {
	*ptr_Z = ptr_Z[1] + *xn * (*ptr_b / a0) - *yn * (*ptr_a / a0);
	ptr_b++; ptr_a++; ptr_Z++;
      }
    /* Calculate last delay */
      *ptr_Z = *xn * (*ptr_b / a0) - *yn * (*ptr_a / a0);
    }
    else {
      *yn = *xn * (*ptr_b / a0);
    }

    ptr_y += stride_Y;             /* Move to next input/output point */
    ptr_x += stride_X;     
  }
}

 
static void CFLOAT_filt(char *b, char *a, char *x, char *y, char *Z, int len_b, unsigned int len_x, int stride_X, int stride_Y ) {
  char *ptr_x = x, *ptr_y = y;
  float *ptr_Z, *ptr_b;
  float *ptr_a;
  float *xn, *yn;
  float a0r = ((float *)a)[0];
  float a0i = ((float *)a)[1];
  float a0_mag, tmpr, tmpi;
  int k, n;

  a0_mag = a0r*a0r + a0i*a0i;
  for (k = 0; k < len_x; k++) {
    ptr_b = (float *)b;          /* Reset a and b pointers */
    ptr_a = (float *)a;
    xn = (float *)ptr_x;
    yn = (float *)ptr_y;
    if (len_b > 1) {
      ptr_Z = ((float *)Z);
      tmpr = ptr_b[0]*a0r + ptr_b[1]*a0i;
      tmpi = ptr_b[1]*a0r - ptr_b[0]*a0i;
      /* Calculate first delay (output) */
      yn[0] = ptr_Z[0] + (tmpr * xn[0] - tmpi * xn[1])/a0_mag;    
      yn[1] = ptr_Z[1] + (tmpi * xn[0] + tmpr * xn[1])/a0_mag;
      ptr_b += 2; ptr_a += 2;
      /* Fill in middle delays */
      for (n = 0; n < len_b - 2; n++) {
	tmpr = ptr_b[0]*a0r + ptr_b[1]*a0i;
	tmpi = ptr_b[1]*a0r - ptr_b[0]*a0i;
	ptr_Z[0] = ptr_Z[2] + (tmpr * xn[0] - tmpi * xn[1])/a0_mag;
	ptr_Z[1] = ptr_Z[3] + (tmpi * xn[0] + tmpr * xn[1])/a0_mag;
	tmpr = ptr_a[0]*a0r + ptr_a[1]*a0i;
	tmpi = ptr_a[1]*a0r - ptr_a[0]*a0i;
	ptr_Z[0] -= (tmpr * yn[0] - tmpi * yn[1])/a0_mag;
	ptr_Z[1] -= (tmpi * yn[0] + tmpr * yn[1])/a0_mag;
	ptr_b += 2; ptr_a += 2; ptr_Z += 2;
      }
      /* Calculate last delay */

      tmpr = ptr_b[0]*a0r + ptr_b[1]*a0i;
      tmpi = ptr_b[1]*a0r - ptr_b[0]*a0i;
      ptr_Z[0] = (tmpr * xn[0] - tmpi * xn[1])/a0_mag;
      ptr_Z[1] = (tmpi * xn[0] + tmpr * xn[1])/a0_mag;
      tmpr = ptr_a[0]*a0r + ptr_a[1]*a0i;
      tmpi = ptr_a[1]*a0r - ptr_a[0]*a0i;
      ptr_Z[0] -= (tmpr * yn[0] - tmpi * yn[1])/a0_mag;
      ptr_Z[1] -= (tmpi * yn[0] + tmpr * yn[1])/a0_mag;
    }
    else {
      tmpr = ptr_b[0]*a0r + ptr_b[1]*a0i;
      tmpi = ptr_b[1]*a0r - ptr_b[0]*a0i;
      yn[0] = (tmpr * xn[0] - tmpi * xn[1])/a0_mag;
      yn[1] = (tmpi * xn[0] + tmpr * xn[1])/a0_mag;      
    }

    ptr_y += stride_Y;             /* Move to next input/output point */
    ptr_x += stride_X;     

  }
}


static void CDOUBLE_filt(char *b, char *a, char *x, char *y, char *Z, int len_b, unsigned int len_x, int stride_X, int stride_Y ) {
  char *ptr_x = x, *ptr_y = y;
  double *ptr_Z, *ptr_b;
  double *ptr_a;
  double *xn, *yn;
  double a0r = ((double *)a)[0];
  double a0i = ((double *)a)[1];
  double a0_mag, tmpr, tmpi;
  int k, n;

  a0_mag = a0r*a0r + a0i*a0i;
  for (k = 0; k < len_x; k++) {
    ptr_b = (double *)b;          /* Reset a and b pointers */
    ptr_a = (double *)a;
    xn = (double *)ptr_x;
    yn = (double *)ptr_y;
    if (len_b > 1) {
      ptr_Z = ((double *)Z);
      tmpr = ptr_b[0]*a0r + ptr_b[1]*a0i;
      tmpi = ptr_b[1]*a0r - ptr_b[0]*a0i;
      /* Calculate first delay (output) */
      yn[0] = ptr_Z[0] + (tmpr * xn[0] - tmpi * xn[1])/a0_mag;    
      yn[1] = ptr_Z[1] + (tmpi * xn[0] + tmpr * xn[1])/a0_mag;
      ptr_b += 2; ptr_a += 2;
      /* Fill in middle delays */
      for (n = 0; n < len_b - 2; n++) {
	tmpr = ptr_b[0]*a0r + ptr_b[1]*a0i;
	tmpi = ptr_b[1]*a0r - ptr_b[0]*a0i;
	ptr_Z[0] = ptr_Z[2] + (tmpr * xn[0] - tmpi * xn[1])/a0_mag;
	ptr_Z[1] = ptr_Z[3] + (tmpi * xn[0] + tmpr * xn[1])/a0_mag;
	tmpr = ptr_a[0]*a0r + ptr_a[1]*a0i;
	tmpi = ptr_a[1]*a0r - ptr_a[0]*a0i;
	ptr_Z[0] -= (tmpr * yn[0] - tmpi * yn[1])/a0_mag;
	ptr_Z[1] -= (tmpi * yn[0] + tmpr * yn[1])/a0_mag;
	ptr_b += 2; ptr_a += 2; ptr_Z += 2;
      }
      /* Calculate last delay */
      tmpr = ptr_b[0]*a0r + ptr_b[1]*a0i;
      tmpi = ptr_b[1]*a0r - ptr_b[0]*a0i;
      ptr_Z[0] = (tmpr * xn[0] - tmpi * xn[1])/a0_mag;
      ptr_Z[1] = (tmpi * xn[0] + tmpr * xn[1])/a0_mag;
      tmpr = ptr_a[0]*a0r + ptr_a[1]*a0i;
      tmpi = ptr_a[1]*a0r - ptr_a[0]*a0i;
      ptr_Z[0] -= (tmpr * yn[0] - tmpi * yn[1])/a0_mag;
      ptr_Z[1] -= (tmpi * yn[0] + tmpr * yn[1])/a0_mag;
    }
    else {
      tmpr = ptr_b[0]*a0r + ptr_b[1]*a0i;
      tmpi = ptr_b[1]*a0r - ptr_b[0]*a0i;
      yn[0] = (tmpr * xn[0] - tmpi * xn[1])/a0_mag;
      yn[1] = (tmpi * xn[0] + tmpr * xn[1])/a0_mag;      
    } 
    ptr_y += stride_Y;             /* Move to next input/output point */
    ptr_x += stride_X;     
  }
}




/* This routine expects the vector to a to begin with a non-zero value */

static void RawFilter(Generic_Vector Vb, Generic_Vector Va, Generic_Array X, Generic_Array Y, Generic_Array *Vi, Generic_Array *Vf, BasicFilterFunction *filt, int along_dimen) {

  int k, n, byte_len;
  int ndims, num_loops, count, stride_X, stride_Y, incr = 1, filt_size;
  unsigned int len;
  int *loop_index, *loop_strides_X, *loop_strides_Y, *loop_strides_Vi;
  int *loop_strides_Vf, *max_index;
  char *ptrX, *ptrY, *ptrVi=NULL, *ptrVf, *ptra, *ptrb;
  char *pa1, *pa2, *pb1, *pb2;

  /* Make dimension array for the looping index that has 
   * along_dimen removed.  Also remove along_dimension from strides
   * vectors.
   */
  ndims = X.nd - 1; 
  if (ndims < 1)
    ndims = 1;                 /* make sure at least 1 for allocation */
  byte_len = ndims*sizeof(int);
  max_index = (int *)malloc(byte_len);
  loop_index = (int *)malloc(byte_len);
  memset(loop_index, 0, byte_len);
  loop_strides_X = (int *)malloc(byte_len);
  loop_strides_Y = (int *)malloc(byte_len);
  loop_strides_Vi = (int *)malloc(byte_len);
  loop_strides_Vf = (int *)malloc(byte_len);
  num_loops = 1;
  count = 0;
  for (k=0; k < X.nd; k++) {
    if (k != along_dimen) {
      loop_strides_X[count] = X.strides[k];
      loop_strides_Y[count] = Y.strides[k];
      max_index[count++] = X.dimensions[k];
      num_loops *= X.dimensions[k];
      if (Vi != NULL) {
	loop_strides_Vi[k] = Vi->strides[k];
	loop_strides_Vf[k] = Vf->strides[k];
      }
    }
  }

  /* Make same length buffers for the a and b filter coefficients.
     and fill with received values and zero for others.
   */
#define max(x,y) ((x) > (y) ? (x) : (y))
  filt_size = max(Vb.numels, Va.numels);
#undef max

  ptra = (char *)malloc(filt_size*Va.elsize);
  ptrb = (char *)malloc(filt_size*Vb.elsize);
  pa1 = ptra; pa2 = Va.data;
  pb1 = ptrb; pb2 = Vb.data;
  for (k = 0; k < filt_size; k++) {
    if (k < Va.numels)
      memcpy(pa1, pa2, Va.elsize);
    else
      memcpy(pa1, Va.zero, Va.elsize);
    if (k < Vb.numels)
      memcpy(pb1, pb2, Vb.elsize);
    else
      memcpy(pb1, Vb.zero, Vb.elsize);
    pa1 += Va.elsize; 
    pa2 += Va.elsize;
    pb1 += Vb.elsize;
    pb2 += Vb.elsize;
  }    
    
  len = X.dimensions[along_dimen];
  stride_X = X.strides[along_dimen]; 
  stride_Y = Y.strides[along_dimen];

  if (Vi != NULL) {          /* Copy initial data to final data space */
    ptrVi = Vi->data;         /*  Where it will be used in calcluations */
    ptrVf = Vf->data;
    pa1 = ptrVi; pa2 = ptrVf;
    for (k=0; k < filt_size - 1; k++) {
      memcpy(pa2, pa1, Vi->elsize);
      pa2 += Vf->elsize;              /* Contiguous */
      pa1 += Vi->strides[along_dimen]; /* Possibly non-contiguous */
    }
  }
  else { /* Initial conditions are zero, and final conditions will be ignored 
	    Create some area to use for the filter delays and fill with zeros.
          */
    ptrVf = (char *)malloc((filt_size-1)*X.elsize);
    pa2 = ptrVf;
    for (k = 0; k < filt_size - 1; k++) {
      memcpy(pa2, X.zero, X.elsize);
      pa2 += X.elsize;
    }
  }

  ptrX = X.data;
  ptrY = Y.data;

  while (num_loops--) { 
    filt(ptrb, ptra, ptrX, ptrY, ptrVf, filt_size, len, stride_X, stride_Y);
    incr = increment(loop_index, ndims, max_index);  /* Returns number of N-D indices incremented. */
    k = ndims - incr;

    if (num_loops > 0) {
    /* Adjust index array and move pointers to right place for next iteration */

      ptrX += loop_strides_X[k];         /* Stride information */
      ptrY += loop_strides_Y[k];

      if (Vi != NULL) {
	ptrVi += loop_strides_Vi[k];
	ptrVf += loop_strides_Vf[k];            /* Move to new storage area */
	pa1 = ptrVi; pa2 = ptrVf;
	for (n=0; n < filt_size - 1; n++) {   /* Copy new initial conditions */
	  memcpy(pa2, pa1, Vi->elsize);
	  pa2 += Vf->elsize;               /* Contiguous */
	  pa1 += Vi->strides[along_dimen]; /* Possibly non-contiguous */
	}
      }
      else {             /* Re-initialize delays to zero */
	pa2 = ptrVf;
	for (n = 0; n < filt_size - 1; n++) {
	  memcpy(pa2, X.zero, X.elsize);
	  pa2 += X.elsize;
	  
	}
      }    
    }
  }
  if (Vi == NULL) {
     free(ptrVf);
  }

  free(ptra); free(ptrb);
  free(loop_index); free(max_index);
  free(loop_strides_X); free(loop_strides_Y); 
  free(loop_strides_Vi); free(loop_strides_Vf);
}



/********************************************************
 *
 *  Code taken from remez.c by Erik Kvaleberg which was 
 *    converted from an original FORTRAN by
 *
 * AUTHORS: JAMES H. MCCLELLAN
 *
 *         DEPARTMENT OF ELECTRICAL ENGINEERING AND COMPUTER SCIENCE
 *         MASSACHUSETTS INSTITUTE OF TECHNOLOGY
 *         CAMBRIDGE, MASS. 02139
 *
 *         THOMAS W. PARKS
 *         DEPARTMENT OF ELECTRICAL ENGINEERING
 *         RICE UNIVERSITY
 *         HOUSTON, TEXAS 77001
 *
 *         LAWRENCE R. RABINER
 *         BELL LABORATORIES
 *         MURRAY HILL, NEW JERSEY 07974
 *
 *  
 *  Adaptation to C by 
 *      egil kvaleberg
 *      husebybakken 14a
 *      0379 oslo, norway
 *  Email:
 *      egil@kvaleberg.no
 *  Web:
 *      http://www.kvaleberg.com/
 * 
 *
 *********************************************************/


#define BANDPASS       1
#define DIFFERENTIATOR 2
#define HILBERT        3

#define GOBACK goto
#define DOloop(a,from,to) for ( (a) = (from); (a) <= (to); ++(a))
#define PI    3.14159265358979323846
#define TWOPI (PI+PI)

/*
 *-----------------------------------------------------------------------
 * FUNCTION: lagrange_interp (d)
 *  FUNCTION TO CALCULATE THE LAGRANGE INTERPOLATION
 *  COEFFICIENTS FOR USE IN THE FUNCTION gee.
 *-----------------------------------------------------------------------
 */
static double lagrange_interp(int k, int n, int m, double *x)
{
    int j, l;
    double q, retval;

    retval = 1.0;
    q = x[k];
    DOloop(l,1,m) {
	for (j = l; j <= n; j += m) {
	    if (j != k)
		retval *= 2.0 * (q - x[j]);
	}
    }
    return 1.0 / retval;
}

/*
 *-----------------------------------------------------------------------
 * FUNCTION: freq_eval (gee)
 *  FUNCTION TO EVALUATE THE FREQUENCY RESPONSE USING THE
 *  LAGRANGE INTERPOLATION FORMULA IN THE BARYCENTRIC FORM
 *-----------------------------------------------------------------------
 */
static double freq_eval(int k, int n, double *grid, double *x, double *y, double *ad)
{
    int j;
    double p,c,d,xf;

    d = 0.0;
    p = 0.0;
    xf = cos(TWOPI * grid[k]);

    DOloop(j,1,n) {
	c = ad[j] / (xf - x[j]);
	d += c;
	p += c * y[j];
    }

    return p/d;
}


/*
 *-----------------------------------------------------------------------
 * SUBROUTINE: remez
 *  THIS SUBROUTINE IMPLEMENTS THE REMEZ EXCHANGE ALGORITHM
 *  FOR THE WEIGHTED CHEBYSHEV APPROXIMATION OF A CONTINUOUS
 *  FUNCTION WITH A SUM OF COSINES.  INPUTS TO THE SUBROUTINE
 *  ARE A DENSE GRID WHICH REPLACES THE FREQUENCY AXIS, THE
 *  DESIRED FUNCTION ON THIS GRID, THE WEIGHT FUNCTION ON THE
 *  GRID, THE NUMBER OF COSINES, AND AN INITIAL GUESS OF THE
 *  EXTREMAL FREQUENCIES.  THE PROGRAM MINIMIZES THE CHEBYSHEV
 *  ERROR BY DETERMINING THE BSMINEST LOCATION OF THE EXTREMAL
 *  FREQUENCIES (POINTS OF MAXIMUM ERROR) AND THEN CALCULATES
 *  THE COEFFICIENTS OF THE BEST APPROXIMATION.
 *-----------------------------------------------------------------------
 */
static int remez(double *dev, double des[], double grid[], double edge[],  
	   double wt[], int ngrid, int nbands, int iext[], double alpha[],
	   int nfcns, int itrmax, double *work, int dimsize)
		/* dev, iext, alpha                         are output types */
		/* des, grid, edge, wt, ngrid, nbands, nfcns are input types */
{
    int k, k1, kkk, kn, knz, klow, kup, nz, nzz, nm1;
    int cn;
    int j, jchnge, jet, jm1, jp1;
    int l, luck=0, nu, nut, nut1=0, niter;

    double ynz=0.0, comp=0.0, devl, gtemp, fsh, y1=0.0, err, dtemp, delf, dnum, dden;
    double aa=0.0, bb=0.0, ft, xe, xt;

    static double *a, *p, *q;
    static double *ad, *x, *y;

    a = work; p = a + dimsize+1; q = p + dimsize+1; 
    ad = q + dimsize+1; x = ad + dimsize+1; y = x + dimsize+1;
    devl = -1.0;
    nz  = nfcns+1;
    nzz = nfcns+2;
    niter = 0;

    do {
    L100:
	iext[nzz] = ngrid + 1;
	++niter;

	if (niter > itrmax) break;

	/* printf("ITERATION %2d: ",niter); */

	DOloop(j,1,nz) {
	    x[j] = cos(grid[iext[j]]*TWOPI);
	}
	jet = (nfcns-1) / 15 + 1;

	DOloop(j,1,nz) {
	    ad[j] = lagrange_interp(j,nz,jet,x);
	}

	dnum = 0.0;
	dden = 0.0;
	k = 1;

	DOloop(j,1,nz) {
	    l = iext[j];
	    dnum += ad[j] * des[l];
	    dden += (double)k * ad[j] / wt[l];
	    k = -k;
	}
	*dev = dnum / dden;

	/* printf("DEVIATION = %lg\n",*dev); */

	nu = 1;
	if ( (*dev) > 0.0 ) nu = -1;
	(*dev) = -(double)nu * (*dev);
	k = nu;
	DOloop(j,1,nz) {
	    l = iext[j];
	    y[j] = des[l] + (double)k * (*dev) / wt[l];
	    k = -k;
	}
	if ( (*dev) <= devl ) {
	    /* finished */
	    return -1;
	}
	devl = (*dev);
	jchnge = 0;
	k1 = iext[1];
	knz = iext[nz];
	klow = 0;
	nut = -nu;
	j = 1;

    /*
     * SEARCH FOR THE EXTREMAL FREQUENCIES OF THE BEST APPROXIMATION
     */

    L200:
	if (j == nzz) ynz = comp;
	if (j >= nzz) goto L300;
	kup = iext[j+1];
	l = iext[j]+1;
	nut = -nut;
	if (j == 2) y1 = comp;
	comp = (*dev);
	if (l >= kup) goto L220;
	err = (freq_eval(l,nz,grid,x,y,ad)-des[l]) * wt[l];
	if (((double)nut*err-comp) <= 0.0) goto L220;
	comp = (double)nut * err;
    L210:
	if (++l >= kup) goto L215;
	err = (freq_eval(l,nz,grid,x,y,ad)-des[l]) * wt[l];
	if (((double)nut*err-comp) <= 0.0) goto L215;
	comp = (double)nut * err;
	GOBACK L210;

    L215:
	iext[j++] = l - 1;
	klow = l - 1;
	++jchnge;
	GOBACK L200;

    L220:
	--l;
    L225:
	if (--l <= klow) goto L250;
	err = (freq_eval(l,nz,grid,x,y,ad)-des[l]) * wt[l];
	if (((double)nut*err-comp) > 0.0) goto L230;
	if (jchnge <= 0) goto L225;
	goto L260;

    L230:
	comp = (double)nut * err;
    L235:
	if (--l <= klow) goto L240;
	err = (freq_eval(l,nz,grid,x,y,ad)-des[l]) * wt[l];
	if (((double)nut*err-comp) <= 0.0) goto L240;
	comp = (double)nut * err;
	GOBACK L235;
    L240:
	klow = iext[j];
	iext[j] = l+1;
	++j;
	++jchnge;
	GOBACK L200;

    L250:
	l = iext[j]+1;
	if (jchnge > 0) GOBACK L215;

    L255:
	if (++l >= kup) goto L260;
	err = (freq_eval(l,nz,grid,x,y,ad)-des[l]) * wt[l];
	if (((double)nut*err-comp) <= 0.0) GOBACK L255;
	comp = (double)nut * err;

	GOBACK L210;
    L260:
	klow = iext[j++];
	GOBACK L200;

    L300:
	if (j > nzz) goto L320;
	if (k1 > iext[1] ) k1 = iext[1];
	if (knz < iext[nz]) knz = iext[nz];
	nut1 = nut;
	nut = -nu;
	l = 0;
	kup = k1;
	comp = ynz*(1.00001);
	luck = 1;
    L310:
	if (++l >= kup) goto L315;
	err = (freq_eval(l,nz,grid,x,y,ad)-des[l]) * wt[l];
	if (((double)nut*err-comp) <= 0.0) GOBACK L310;
	comp = (double) nut * err;
	j = nzz;
	GOBACK L210;

    L315:
	luck = 6;
	goto L325;

    L320:
	if (luck > 9) goto L350;
	if (comp > y1) y1 = comp;
	k1 = iext[nzz];
    L325:
	l = ngrid+1;
	klow = knz;
	nut = -nut1;
	comp = y1*(1.00001);
    L330:
	if (--l <= klow) goto L340;
	err = (freq_eval(l,nz,grid,x,y,ad)-des[l]) * wt[l];
	if (((double)nut*err-comp) <= 0.0) GOBACK L330;
	j = nzz;
	comp = (double) nut * err;
	luck = luck + 10;
	GOBACK L235;
    L340:
	if (luck == 6) goto L370;
	DOloop(j,1,nfcns) {
	    iext[nzz-j] = iext[nz-j];
	}
	iext[1] = k1;
	GOBACK L100;
    L350:
	kn = iext[nzz];
	DOloop(j,1,nfcns) iext[j] = iext[j+1];
	iext[nz] = kn;

	GOBACK L100;
    L370:
	;
    } while (jchnge > 0);

/*
 *    CALCULATION OF THE COEFFICIENTS OF THE BEST APPROXIMATION
 *    USING THE INVERSE DISCRETE FOURIER TRANSFORM
 */
    nm1 = nfcns - 1;
    fsh = 1.0e-06;
    gtemp = grid[1];
    x[nzz] = -2.0;
    cn  = 2*nfcns - 1;
    delf = 1.0/cn;
    l = 1;
    kkk = 0;

    if (edge[1] == 0.0 && edge[2*nbands] == 0.5) kkk = 1;

    if (nfcns <= 3) kkk = 1;
    if (kkk !=     1) {
	dtemp = cos(TWOPI*grid[1]);
	dnum  = cos(TWOPI*grid[ngrid]);
	aa    = 2.0/(dtemp-dnum);
	bb    = -(dtemp+dnum)/(dtemp-dnum);
    }

    DOloop(j,1,nfcns) {
	ft = (j - 1) * delf;
	xt = cos(TWOPI*ft);
	if (kkk != 1) {
	    xt = (xt-bb)/aa;
#if 0
	    /*XX* ckeck up !! */
	    xt1 = sqrt(1.0-xt*xt);
	    ft = atan2(xt1,xt)/TWOPI;
#else
	    ft = acos(xt)/TWOPI;
#endif
	}
L410:
	xe = x[l];
	if (xt > xe) goto L420;
	if ((xe-xt) < fsh) goto L415;
	++l;
	GOBACK L410;
L415:
	a[j] = y[l];
	goto L425;
L420:
	if ((xt-xe) < fsh) GOBACK L415;
	grid[1] = ft;
	a[j] = freq_eval(1,nz,grid,x,y,ad);
L425:
	if (l > 1) l = l-1;
    }

    grid[1] = gtemp;
    dden = TWOPI / cn;
    DOloop (j,1,nfcns) {
	dtemp = 0.0;
	dnum = (j-1) * dden;
	if (nm1 >= 1) {
	    DOloop(k,1,nm1) {
		dtemp += a[k+1] * cos(dnum*k);
	    }
	}
	alpha[j] = 2.0 * dtemp + a[1];
    }

    DOloop(j,2,nfcns) alpha[j] *= 2.0 / cn;
    alpha[1] /= cn;

    if (kkk != 1) {
	p[1] = 2.0*alpha[nfcns]*bb+alpha[nm1];
	p[2] = 2.0*aa*alpha[nfcns];
	q[1] = alpha[nfcns-2]-alpha[nfcns];
	DOloop(j,2,nm1) {
	    if (j >= nm1) {
		aa *= 0.5;
		bb *= 0.5;
	    }
	    p[j+1] = 0.0;
	    DOloop(k,1,j) {
		a[k] = p[k];
		p[k] = 2.0 * bb * a[k];
	    }
	    p[2] += a[1] * 2.0 *aa;
	    jm1 = j - 1;
	    DOloop(k,1,jm1) p[k] += q[k] + aa * a[k+1];
	    jp1 = j + 1;
	    DOloop(k,3,jp1) p[k] += aa * a[k-1];

	    if (j != nm1) {
		DOloop(k,1,j) q[k] = -a[k];
		q[1] += alpha[nfcns - 1 - j];
	    }
	}
	DOloop(j,1,nfcns) alpha[j] = p[j];
    }

    if (nfcns <= 3) {
	  alpha[nfcns+1] = alpha[nfcns+2] = 0.0;
    }
    return 0;
}


/*
 *-----------------------------------------------------------------------
 * FUNCTION: eff
 *  FUNCTION TO CALCULATE THE DESIRED MAGNITUDE RESPONSE
 *  AS A FUNCTION OF FREQUENCY.
 *  AN ARBITRARY FUNCTION OF FREQUENCY CAN BE
 *  APPROXIMATED IF THE USER REPLACES THIS FUNCTION
 *  WITH THE APPROPRIATE CODE TO EVALUATE THE IDEAL
 *  MAGNITUDE.  NOTE THAT THE PARAMETER FREQ IS THE
 *  VALUE OF NORMALIZED FREQUENCY NEEDED FOR EVALUATION.
 *-----------------------------------------------------------------------
 */
static double eff(double freq, double *fx, int lband, int jtype)
{
      if (jtype != 2) return fx[lband];
      else            return fx[lband] * freq;
}

/*
 *-----------------------------------------------------------------------
 * FUNCTION: wate
 *  FUNCTION TO CALCULATE THE WEIGHT FUNCTION AS A FUNCTION
 *  OF FREQUENCY.  SIMILAR TO THE FUNCTION eff, THIS FUNCTION CAN
 *  BE REPLACED BY A USER-WRITTEN ROUTINE TO CALCULATE ANY
 *  DESIRED WEIGHTING FUNCTION.
 *-----------------------------------------------------------------------
 */
static double wate(double freq, double *fx, double *wtx, int lband, int jtype)
{
      if (jtype != 2)          return wtx[lband];
      if (fx[lband] >= 0.0001) return wtx[lband] / freq;
      return                          wtx[lband];
}

/*********************************************************/

/*  This routine accepts basic input information and puts it in 
 *  the form expected by remez.

 *  Adpated from main() by Travis Oliphant
 */

static int pre_remez(double *h2, int numtaps, int numbands, double *bands, double *response, double *weight, int type, int maxiter, int grid_density) {
  
  int jtype, nbands, nfilt, lgrid, nz;
  int neg, nodd, nm1;
  int j, k, l, lband, dimsize;
  double delf, change, fup, temp;
  double *tempstor, *edge, *h, *fx, *wtx;
  double *des, *grid, *wt, *alpha, *work;
  double dev;
  int ngrid;
  int *iext;
  int nfcns, wrksize, total_dsize, total_isize;

  lgrid = grid_density;
  dimsize = (int) ceil(numtaps/2.0 + 2);
  wrksize = grid_density * dimsize;
  nfilt = numtaps;
  jtype = type; nbands = numbands;
  /* Note:  code assumes these arrays start at 1 */
  edge = bands-1; 
  h = h2 - 1;
  fx = response - 1;
  wtx = weight - 1;

  total_dsize = (dimsize+1)*7 + 3*(wrksize+1);
  total_isize = (dimsize+1);
  /* Need space for:  (all arrays ignore the first element).

     des  (wrksize+1)
     grid (wrksize+1)
     wt   (wrksize+1)
     iext (dimsize+1)   (integer)
     alpha (dimsize+1)
     work  (dimsize+1)*6 

  */
  tempstor = malloc((total_dsize)*sizeof(double)+(total_isize)*sizeof(int));
  if (tempstor == NULL) return -2;

  des = tempstor; grid = des + wrksize+1;
  wt = grid + wrksize+1; alpha = wt + wrksize+1;
  work = alpha + dimsize+1; iext = (int *)(work + (dimsize+1)*6);

  /* Set up problem on dense_grid */

  neg = 1;
  if (jtype == 1) neg = 0;
  nodd = nfilt % 2;
  nfcns = nfilt / 2;
  if (nodd == 1 && neg == 0) nfcns = nfcns + 1;

    /*
     * SET UP THE DENSE GRID. THE NUMBER OF POINTS IN THE GRID
     * IS (FILTER LENGTH + 1)*GRID DENSITY/2
     */
    grid[1] = edge[1];
    delf = lgrid * nfcns;
    delf = 0.5 / delf;
    if (neg != 0) {
	if (edge[1] < delf) grid[1] = delf;
    }
    j = 1;
    l = 1;
    lband = 1;

    /*
     * CALCULATE THE DESIRED MAGNITUDE RESPONSE AND THE WEIGHT
     * FUNCTION ON THE GRID
     */
    for (;;) {
	fup = edge[l + 1];
	do {
	    temp = grid[j];
	    des[j] = eff(temp,fx,lband,jtype);
	    wt[j] = wate(temp,fx,wtx,lband,jtype);
	    if (++j > wrksize) { free(tempstor); return -1;} /* too many points, or too dense grid */
	    grid[j] = temp + delf;
	} while (grid[j] <= fup);

	grid[j-1] = fup;
	des[j-1] = eff(fup,fx,lband,jtype);
	wt[j-1] = wate(fup,fx,wtx,lband,jtype);
	++lband;
	l += 2;
	if (lband > nbands) break;
	grid[j] = edge[l];
    }

    ngrid = j - 1;
    if (neg == nodd) {
	if (grid[ngrid] > (0.5-delf)) --ngrid;
    }

    /*
     * SET UP A NEW APPROXIMATION PROBLEM WHICH IS EQUIVALENT
     * TO THE ORIGINAL PROBLEM
     */
    if (neg <= 0) {
	if (nodd != 1) {
	    DOloop(j,1,ngrid) {
		change = cos(PI*grid[j]);
		des[j] = des[j] / change;
		wt[j]  = wt[j] * change;
	    }
	}
    } else {
	if (nodd != 1) {
	    DOloop(j,1,ngrid) {
		change = sin(PI*grid[j]);
		des[j] = des[j] / change;
		wt[j]  = wt[j]  * change;
	    }
	} else {
	    DOloop(j,1,ngrid) {
		change = sin(TWOPI * grid[j]);
		des[j] = des[j] / change;
		wt[j]  = wt[j]  * change;
	    }
	}
    }

    /*XX*/
    temp = (double)(ngrid-1) / (double)nfcns;
    DOloop(j,1,nfcns) {
	iext[j] = (int)((j-1)*temp) + 1; /* round? !! */
    }
    iext[nfcns+1] = ngrid;
    nm1 = nfcns - 1;
    nz  = nfcns + 1;

    if (remez(&dev, des, grid, edge, wt, ngrid, numbands, iext, alpha, nfcns, maxiter, work, dimsize) < 0) { free(tempstor); return -1; }

    /*
     * CALCULATE THE IMPULSE RESPONSE.
     */
    if (neg <= 0) {

	if (nodd != 0) {
	    DOloop(j,1,nm1) {
		h[j] = 0.5 * alpha[nz-j];
	    }
	    h[nfcns] = alpha[1];
	} else {
	    h[1] = 0.25 * alpha[nfcns];
	    DOloop(j,2,nm1) {
		h[j] = 0.25 * (alpha[nz-j] + alpha[nfcns+2-j]);
	    }
	    h[nfcns] = 0.5*alpha[1] + 0.25*alpha[2];
	}
    } else {
	if (nodd != 0) {
	    h[1] = 0.25 * alpha[nfcns];
	    h[2] = 0.25 * alpha[nm1];
	    DOloop(j,3,nm1) {
		h[j] = 0.25 * (alpha[nz-j] - alpha[nfcns+3-j]);
	    }
	    h[nfcns] = 0.5 * alpha[1] - 0.25 * alpha[3];
	    h[nz] = 0.0;
	} else {
	    h[1] = 0.25 * alpha[nfcns];
	    DOloop(j,2,nm1) {
		h[j] = 0.25 * (alpha[nz-j] - alpha[nfcns+2-j]);
	    }
	    h[nfcns] = 0.5 * alpha[1] - 0.25 * alpha[2];
	}
    }

    DOloop(j,1,nfcns){
        k = nfilt + 1 - j;
        if (neg == 0)
           h[k] = h[j];
        else
           h[k] = -h[j];
    }
    if (neg == 1 && nodd == 1) h[nz] = 0.0;

  free(tempstor);
  return 0;

}

/**************************************************************
 * End of remez routines 
 **************************************************************/


/****************************************************/
/* End of python-independent routines               */
/****************************************************/

static void OBJECT_MultAdd(char *ip1, int is1, char *ip2, int is2, char *op, int *dims1, int *dims2, int ndims, int nels2, int check, int *loop_ind, int *temp_ind, unsigned long *offset) { 
  int i, k, first_time = 1, incr = 1; 
  PyObject *tmp1=NULL, *tmp2=NULL, *tmp=NULL;

  i = nels2;

  temp_ind[ndims-1]--;
  while (i--) { 
    /* Adjust index array and move ptr1 to right place */
    k = ndims - 1;
    while(--incr) {
      temp_ind[k] -= dims2[k] - 1;   /* Return to start for these dimensions */
      k--;
    }
    ip1 += offset[k]*is1;           /* Precomputed offset array */
    temp_ind[k]++;

    if (!(check && index_out_of_bounds(temp_ind,dims1,ndims))) { 
      tmp1 = PyNumber_Multiply(*((PyObject **)ip1),*((PyObject **)ip2));
      if (first_time) {
	tmp = tmp1;
	first_time = 0;
      } else {
	tmp2 = PyNumber_Add(tmp, tmp1);
	Py_XDECREF(tmp);
	tmp = tmp2;
	Py_XDECREF(tmp1);
      }
    } 
    incr = increment(loop_ind, ndims, dims2); 
    ip2 += is2; 
  }
  Py_XDECREF(*((PyObject **)op));
  *((PyObject **)op) = tmp; 
}

#ifdef PyArray_UNSIGNED_TYPES
static MultAddFunction *MultiplyAddFunctions[] = {NULL,UBYTE_MultAdd,SBYTE_MultAdd,SHORT_MultAdd,NULL,INT_MultAdd,NULL,LONG_MultAdd,
FLOAT_MultAdd, DOUBLE_MultAdd, 
CFLOAT_MultAdd, CDOUBLE_MultAdd, OBJECT_MultAdd};
#else
static MultAddFunction *MultiplyAddFunctions[] = {NULL,UBYTE_MultAdd,SBYTE_MultAdd,SHORT_MultAdd,INT_MultAdd,LONG_MultAdd,
FLOAT_MultAdd, DOUBLE_MultAdd, 
CFLOAT_MultAdd, CDOUBLE_MultAdd, OBJECT_MultAdd};
#endif


static void OBJECT_filt(char *b, char *a, char *x, char *y, char *Z, int len_b, unsigned int len_x, int stride_X, int stride_Y ) {
  char *ptr_x = x, *ptr_y = y;
  PyObject  **ptr_Z, **ptr_b;
  PyObject  **ptr_a;
  PyObject  **xn, **yn;
  PyObject  **a0 = (PyObject  **)a;
  PyObject  *tmp1, *tmp2, *tmp3;
  int k, n;

  /* My reference counting might not be right */
  for (k = 0; k < len_x; k++) {
    ptr_b = (PyObject  **)b;          /* Reset a and b pointers */
    ptr_a = (PyObject  **)a;
    xn = (PyObject  **)ptr_x;
    yn = (PyObject  **)ptr_y;
    if (len_b > 1) {
      ptr_Z = ((PyObject  **)Z);
      /* Calculate first delay (output) */
      tmp1 = PyNumber_Multiply(*ptr_b,*xn);
      tmp2 = PyNumber_Divide(tmp1,*a0);
      tmp3 = PyNumber_Add(tmp2,*ptr_Z);
      Py_XDECREF(*yn);
      *yn = tmp3; Py_DECREF(tmp1); Py_DECREF(tmp2);
      ptr_b++; ptr_a++;

      /* Fill in middle delays */
      for (n = 0; n < len_b - 2; n++) {
	tmp1 = PyNumber_Multiply(*xn, *ptr_b);
	tmp2 = PyNumber_Divide(tmp1,*a0);
	tmp3 = PyNumber_Add(tmp2,ptr_Z[1]);
	Py_DECREF(tmp1); Py_DECREF(tmp2);
	tmp1 = PyNumber_Multiply(*yn, *ptr_a);
	tmp2 = PyNumber_Divide(tmp1, *a0); Py_DECREF(tmp1);
	Py_XDECREF(*ptr_Z);
	*ptr_Z = PyNumber_Subtract(tmp3, tmp2); Py_DECREF(tmp2);
	Py_DECREF(tmp3);
	ptr_b++; ptr_a++; ptr_Z++;
      }
    /* Calculate last delay */
      tmp1 = PyNumber_Multiply(*xn,*ptr_b);
      tmp3 = PyNumber_Divide(tmp1,*a0); Py_DECREF(tmp1);
      tmp1 = PyNumber_Multiply(*yn, *ptr_a);
      tmp2 = PyNumber_Divide(tmp1, *a0); Py_DECREF(tmp1);
      Py_XDECREF(*ptr_Z);
      *ptr_Z = PyNumber_Subtract(tmp3,tmp2); Py_DECREF(tmp2);
      Py_DECREF(tmp3);
    }
    else {
      tmp1 = PyNumber_Multiply(*xn,*ptr_b);
      Py_XDECREF(*yn);
      *yn = PyNumber_Divide(tmp1,*a0); Py_DECREF(tmp1);
    }

    ptr_y += stride_Y;             /* Move to next input/output point */
    ptr_x += stride_X;     
  }
}

/************************/
/* N-D Order Filtering. */


static void fill_buffer(char *ip1, PyArrayObject *ap1, PyArrayObject *ap2, char *sort_buffer, int nels2, int check, int *loop_ind, int *temp_ind, unsigned long *offset){ 
  int i, k, incr = 1;
  int ndims = ap1->nd;
  int *dims2 = ap2->dimensions;
  int *dims1 = ap1->dimensions;
  int is1 = ap1->strides[ndims-1];
  int is2 = ap2->strides[ndims-1];
  char *ip2 = ap2->data;
  int elsize = ap1->descr->elsize;

  i = nels2;

  temp_ind[ndims-1]--;
  while (i--) { 
    /* Adjust index array and move ptr1 to right place */
    k = ndims - 1;
    while(--incr) {
      temp_ind[k] -= dims2[k] - 1;   /* Return to start for these dimensions */
      k--;
    }
    ip1 += offset[k]*is1;               /* Precomputed offset array */
    temp_ind[k]++;

    if (!(check && index_out_of_bounds(temp_ind,dims1,ndims)) && memcmp(ip2, ap2->descr->zero, ap2->descr->elsize)) { 
      memcpy(sort_buffer, ip1, elsize);
      sort_buffer += elsize;
    } 
    incr = increment(loop_ind, ndims, dims2);  /* Returns number of N-D indices incremented. */
    ip2 += is2;

  }
  return;
}

#define COMPARE(fname, type) \
int fname(type *ip1, type *ip2) { return *ip1 < *ip2 ? -1 : *ip1 == *ip2 ? 0 : 1; }

COMPARE(DOUBLE_compare, double)
COMPARE(FLOAT_compare, float)
COMPARE(SHORT_compare, short)
COMPARE(INT_compare, int)
COMPARE(LONG_compare, long)
COMPARE(BYTE_compare, char)
COMPARE(UNSIGNEDBYTE_compare, unsigned char)

#ifdef PyArray_UNSIGNED_TYPES
COMPARE(USHORT_compare, unsigned short)
COMPARE(UINT_compare, unsigned int)
#endif


int OBJECT_compare(PyObject **ip1, PyObject **ip2) {
        return PyObject_Compare(*ip1, *ip2);
}

typedef int (*CompareFunction) Py_FPROTO((const void *, const void *));


CompareFunction compare_functions[] = {NULL,
(CompareFunction)UNSIGNEDBYTE_compare,(CompareFunction)BYTE_compare,(CompareFunction)SHORT_compare,
#ifdef PyArray_UNSIGNED_TYPES
(CompareFunction)USHORT_compare,
#endif
(CompareFunction)INT_compare,
#ifdef PyArray_UNSIGNED_TYPES
(CompareFunction)UINT_compare,
#endif
(CompareFunction)LONG_compare,
(CompareFunction)FLOAT_compare,(CompareFunction)DOUBLE_compare,
NULL,NULL,(CompareFunction)OBJECT_compare};


PyObject *PyArray_OrderFilterND(PyObject *op1, PyObject *op2, int order) {
	PyArrayObject *ap1=NULL, *ap2=NULL, *ret=NULL;
	int *a_ind, *b_ind, *temp_ind, *mode_dep, *check_ind;
	unsigned long *offsets, offset1;
	long *offsets2;
	int i, n2, n2_nonzero, k, check, incr = 1;
	int typenum, bytes_in_array;
	int is1, os;
	char *op, *ap1_ptr, *ap2_ptr, *sort_buffer;
	int *ret_ind;
	CompareFunction compare_func;

	/* Get Array objects from input */
	typenum = PyArray_ObjectType(op1, 0);  
	typenum = PyArray_ObjectType(op2, typenum);
	
	ap1 = (PyArrayObject *)PyArray_ContiguousFromObject(op1, typenum, 0, 0);
	if (ap1 == NULL) return NULL;
	ap2 = (PyArrayObject *)PyArray_ContiguousFromObject(op2, typenum, 0, 0);
	if (ap2 == NULL) goto fail;

	if (ap1->nd != ap2->nd) {
	  PyErr_SetString(PyExc_ValueError, "All input arrays must have the same number of dimensions.");
	  goto fail;
	}

	n2 = PyArray_Size((PyObject *)ap2);
	n2_nonzero = 0;
	ap2_ptr = ap2->data;
	/* Find out the number of non-zero entries in domain (allows for
	 *  different shapped rank-filters to be used besides just rectangles)
	 */
	for (k=0; k < n2; k++) {
	  n2_nonzero += (memcmp(ap2_ptr,ap2->descr->zero,ap2->descr->elsize) != 0);
	  ap2_ptr += ap2->descr->elsize;
	}

	if ((order >= n2_nonzero) || (order < 0)) {
	  PyErr_SetString(PyExc_ValueError, "Order must be non-negative and less than number of nonzero elements in domain.");
	  goto fail;
	}
	
	ret = (PyArrayObject *)PyArray_FromDims(ap1->nd, ap1->dimensions, typenum);
	if (ret == NULL) goto fail;
	
	compare_func = compare_functions[ap1->descr->type_num];
	if (compare_func == NULL) {
		PyErr_SetString(PyExc_ValueError, 
			"order_filterND not available for this type");
		goto fail;
	}

	is1 = ap1->descr->elsize;
	
	if (!(sort_buffer = malloc(n2_nonzero*is1))) goto fail;

	op = ret->data; os = ret->descr->elsize;

	op = ret->data;

	bytes_in_array = ap1->nd*sizeof(int);
	mode_dep = malloc(bytes_in_array);
	for (k = 0; k < ap1->nd; k++) { 
	  mode_dep[k] = -((ap2->dimensions[k]-1) >> 1);
	}	

	b_ind = (int *)malloc(bytes_in_array);  /* loop variables */
	memset(b_ind,0,bytes_in_array);
	a_ind = (int *)malloc(bytes_in_array);
	ret_ind = (int *)malloc(bytes_in_array);
	memset(ret_ind,0,bytes_in_array);
	temp_ind = (int *)malloc(bytes_in_array);
	check_ind = (int*)malloc(bytes_in_array);
	offsets = (unsigned long *)malloc(ap1->nd*sizeof(unsigned long));
	offsets2 = (long *)malloc(ap1->nd*sizeof(long));
	offset1 = compute_offsets(offsets,offsets2,ap1->dimensions,ap2->dimensions,ret->dimensions,mode_dep,ap1->nd);
	/* The filtering proceeds by looping through the output array
	   and for each value filling a buffer from the 
	   element-by-element product of the two input arrays.  The buffer
	   is then sorted and the order_th element is kept as output. Index
	   counters are used for book-keeping in the area so that we 
	   can tell where we are in all of the arrays and be sure that 
	   we are not trying to access areas outside the arrays definition.

	   The inner loop is implemented separately but equivalently for each
	   datatype. The outer loop is similar in structure and form to
	   to the inner loop.
	*/
	/* Need to keep track of a ptr to place in big (first) input
	   array where we start the multiplication (we pass over it in the
	   inner loop (and not dereferenced) 
	   if it is pointing outside dataspace)
	*/
	/* Calculate it once and the just move it around appropriately */
	ap1_ptr = ap1->data + offset1*is1;
	for (k=0; k < ap1->nd; k++) {a_ind[k] = mode_dep[k]; check_ind[k] = ap1->dimensions[k] - ap2->dimensions[k] - mode_dep[k] - 1;}
	a_ind[ap1->nd-1]--;
	i = PyArray_Size((PyObject *)ret);
	while (i--) {
	  /* Zero out the sort_buffer (has effect of zero-padding
	     on boundaries). Treat object arrays right.*/
	  ap2_ptr = sort_buffer;
	  for (k=0; k < n2_nonzero; k++) {
	    memcpy(ap2_ptr,ap1->descr->zero,is1);
	    ap2_ptr += is1;
	  }
	    
	  k = ap1->nd - 1;
	  while(--incr) {
	    a_ind[k] -= ret->dimensions[k] - 1;   /* Return to start */
	    k--;
	  }
	  ap1_ptr += offsets2[k]*is1;
	  a_ind[k]++;
	  memcpy(temp_ind, a_ind, bytes_in_array);

	  check = 0; k = -1;
	  while(!check && (++k < ap1->nd))
	    check = check || (ret_ind[k] < -mode_dep[k]) || (ret_ind[k] > check_ind[k]);

	  fill_buffer(ap1_ptr,ap1,ap2,sort_buffer,n2,check,b_ind,temp_ind,offsets);
	  qsort(sort_buffer, n2_nonzero, is1, compare_func);
	  memcpy(op, sort_buffer + order*is1, os);

	  incr = increment(ret_ind,ret->nd,ret->dimensions); /* increment index counter */
	  op += os;   /* increment to next output index */

	}
	free(b_ind); free(a_ind); free(ret_ind);
	free(offsets); free(offsets2); free(temp_ind);
	free(check_ind); free(mode_dep);
	free(sort_buffer);
	
	Py_DECREF(ap1);
	Py_DECREF(ap2);

	return PyArray_Return(ret);

fail:
	Py_XDECREF(ap1);
	Py_XDECREF(ap2);
	Py_XDECREF(ret);
	return NULL;
}




static BasicFilterFunction *BasicFilterFunctions[] = {NULL,NULL,NULL,NULL,NULL,NULL,
#ifdef PyArray_UNSIGNED_TYPES
NULL,
NULL,
#endif
FLOAT_filt, DOUBLE_filt, 
CFLOAT_filt, CDOUBLE_filt, OBJECT_filt }; 
/* There is the start of an OBJECT_filt, but it may need work */


/* Copy data from PyArray to Generic header for use in C routines */
static void Py_copy_info(Generic_Array *gen, PyArrayObject *py_arr) {
        gen->data = py_arr->data;
	gen->nd = py_arr->nd;
	gen->dimensions = py_arr->dimensions;
	gen->elsize = py_arr->descr->elsize;
	gen->strides = py_arr->strides;
	gen->zero = py_arr->descr->zero;
	return;
}

static void Py_copy_info_vec(Generic_Vector *gen, PyArrayObject *py_arr) {
        gen->data = py_arr->data;
	gen->elsize = py_arr->descr->elsize;
	gen->numels = PyArray_Size((PyObject *)py_arr);
	gen->zero = py_arr->descr->zero;
	return;
}

/******************************************/

static char doc_correlateND[] = "out = _correlateND(a,kernel,mode) \n\n   mode = 0 - 'valid', 1 - 'same', \n  2 - 'full' (default)";

static PyObject *sigtools_correlateND(PyObject *dummy, PyObject *args) {
	PyObject *kernel, *a0;
	PyArrayObject *ap1, *ap2, *ret;
	Generic_Array in1, in2, out;
	int *ret_dimens;
	int mode=2, n1, n2, i, typenum;
	MultAddFunction *multiply_and_add_ND;
	
	if (!PyArg_ParseTuple(args, "OO|i", &a0, &kernel, &mode)) return NULL;

	typenum = PyArray_ObjectType(a0, 0);  
	typenum = PyArray_ObjectType(kernel, typenum);
	
	ret = NULL;
	ap1 = (PyArrayObject *)PyArray_ContiguousFromObject(a0, typenum, 0, 0);
	if (ap1 == NULL) return NULL;
	ap2 = (PyArrayObject *)PyArray_ContiguousFromObject(kernel, typenum, 0, 0);
	if (ap2 == NULL) goto fail;

	if (ap1->nd != ap2->nd) {
	  PyErr_SetString(PyExc_ValueError, "Arrays must have the same number of dimensions.");
	  goto fail;
	}

        if (ap1->nd == 0) {  /* Zero-dimensional arrays */
          PyErr_SetString(PyExc_ValueError, "Cannot convolve zero-dimensional arrays.");
          goto fail;
        }
	
	n1 = PyArray_Size((PyObject *)ap1);
	n2 = PyArray_Size((PyObject *)ap2);

	/* Swap if first argument is not the largest */
	if (n1 < n2) { ret = ap1; ap1 = ap2; ap2 = ret; ret = NULL; }
	ret_dimens = malloc(ap1->nd*sizeof(int));
	switch(mode) {
	case 0:
	        for (i = 0; i < ap1->nd; i++) { 
		  ret_dimens[i] = ap1->dimensions[i] - ap2->dimensions[i] + 1;
		  if (ret_dimens[i] < 0) {
		    PyErr_SetString(PyExc_ValueError, "no part of the output is valid, use option 1 (same) or 2 (full) for third argument");
		    goto fail;
		  }
		}
		break;
	case 1:
	        for (i = 0; i < ap1->nd; i++) { ret_dimens[i] = ap1->dimensions[i];}
		break;
	case 2:
	        for (i = 0; i < ap1->nd; i++) { ret_dimens[i] = ap1->dimensions[i] + ap2->dimensions[i] - 1;}
		break;
	default: 
	        PyErr_SetString(PyExc_ValueError, 
			    "mode must be 0 (valid), 1 (same), or 2 (full)");
		goto fail;
	}
	
	ret = (PyArrayObject *)PyArray_FromDims(ap1->nd, ret_dimens, typenum);
	free(ret_dimens);
	if (ret == NULL) goto fail;
	
	multiply_and_add_ND = MultiplyAddFunctions[(int)(ret->descr->type_num)];
	if (multiply_and_add_ND == NULL) {
		PyErr_SetString(PyExc_ValueError, 
			"correlateND not available for this type");
		goto fail;
	}

	/* copy header information to generic structures */
	Py_copy_info(&in1, ap1);
	Py_copy_info(&in2, ap2);
	Py_copy_info(&out, ret);
       
	correlateND(&in1, &in2, &out, multiply_and_add_ND, mode);

	Py_DECREF(ap1);
	Py_DECREF(ap2);
	return PyArray_Return(ret);

fail:
	Py_XDECREF(ap1);
	Py_XDECREF(ap2);
	Py_XDECREF(ret);
	return NULL;	
}


/*******************************************************************/

static char doc_convolve2d[] = "out = _convolve2d(in1, in2, flip, mode, boundary, fillvalue)";

extern int pylab_convolve_2d(char*,int*,char*,int*,char*,int*,int*,int*,int,char*);

static PyObject *sigtools_convolve2d(PyObject *dummy, PyObject *args) {

    PyObject *in1=NULL, *in2=NULL, *fill_value=NULL;
    int mode=2, boundary=0, typenum, flag, flip=1, ret;
    int *aout_dimens, *dims=NULL;
    char zeros[32];  /* Zeros */
    int n1, n2, i;
    PyArrayObject *ain1=NULL, *ain2=NULL, *aout=NULL;
    PyArrayObject *afill=NULL, *newfill=NULL;

    if (!PyArg_ParseTuple(args, "OO|iiiO", &in1, &in2, &flip, &mode, &boundary, &fill_value)) {
        return NULL;
    }

    typenum = PyArray_ObjectType(in1, 0);
    typenum = PyArray_ObjectType(in2, typenum);
    ain1 = (PyArrayObject *)PyArray_FromObject(in1, typenum, 2, 2);
    if (ain1 == NULL) goto fail;
    ain2 = (PyArrayObject *)PyArray_FromObject(in2, typenum, 2, 2);
    if (ain2 == NULL) goto fail;

    if ((boundary != PAD) && (boundary != REFLECT) && (boundary != CIRCULAR))
      PYERR("Incorrect boundary value.");
    if (boundary == PAD) {
	if (fill_value == NULL) {
	    newfill = (PyArrayObject *)PyArray_FromDimsAndData(0, dims, typenum, zeros);
	}
	else {
	    afill = (PyArrayObject *)PyArray_FromObject(fill_value, PyArray_CDOUBLE, 0, 0);
	    if (afill == NULL) goto fail;
	    newfill = (PyArrayObject *)PyArray_Cast(afill, typenum);
	}
	if (newfill == NULL) goto fail;
    }
    else {
	newfill = (PyArrayObject *)PyArray_FromDimsAndData(0, dims, typenum, zeros);
	if (newfill == NULL) goto fail;
    }
    
    n1 = PyArray_Size((PyObject *)ain1);
    n2 = PyArray_Size((PyObject *)ain2);
    
    /* Swap if first argument is not the largest */
    if (n1 < n2) { aout = ain1; ain1 = ain2; ain2 = aout; aout = NULL; }
    aout_dimens = malloc(ain1->nd*sizeof(int));
    switch(mode & OUTSIZE_MASK) {
    case VALID:
	for (i = 0; i < ain1->nd; i++) { 
	    aout_dimens[i] = ain1->dimensions[i] - ain2->dimensions[i] + 1;
	    if (aout_dimens[i] < 0) {
		PyErr_SetString(PyExc_ValueError, "no part of the output is valid, use option 1 (same) or 2 (full) for third argument");
		goto fail;
	    }
	}
	break;
    case SAME:
	for (i = 0; i < ain1->nd; i++) { aout_dimens[i] = ain1->dimensions[i];}
	break;
    case FULL:
	for (i = 0; i < ain1->nd; i++) { aout_dimens[i] = ain1->dimensions[i] + ain2->dimensions[i] - 1;}
	break;
    default: 
	PyErr_SetString(PyExc_ValueError, 
			"mode must be 0 (valid), 1 (same), or 2 (full)");
	goto fail;
    }
	
    aout = (PyArrayObject *)PyArray_FromDims(ain1->nd, aout_dimens, typenum);
    free(aout_dimens);
    if (aout == NULL) goto fail;

    flag = mode + boundary + (typenum << TYPE_SHIFT) + \
      (flip != 0) * FLIP_MASK;
    
    ret = pylab_convolve_2d (DATA(ain1),      /* Input data Ns[0] x Ns[1] */
		             STRIDES(ain1),   /* Input strides */
		             DATA(aout),      /* Output data */
		             STRIDES(aout),   /* Ouput strides */
		             DATA(ain2),      /* coefficients in filter */
		             STRIDES(ain2),   /* coefficients strides */ 
		             DIMS(ain2),      /* Size of kernel Nwin[2] */
			     DIMS(ain1),      /* Size of image Ns[0] x Ns[1] */
		             flag,            /* convolution parameters */
		             DATA(newfill));  /* fill value */


    switch (ret) {
    case 0:
      Py_DECREF(ain1);
      Py_DECREF(ain2);
      Py_XDECREF(afill);
      Py_XDECREF(newfill);
      return (PyObject *)aout;
      break;
    case -5:
    case -4:
      PyErr_SetString(PyExc_ValueError,
		      "convolve2d not available for this type.");
      goto fail;
    case -3:
      PyErr_NoMemory();
      goto fail;
    case -2:
      PyErr_SetString(PyExc_ValueError,
		      "Invalid boundary type.");
      goto fail;
    case -1:
      PyErr_SetString(PyExc_ValueError,
		      "Invalid output flag.");
      goto fail;
    }

fail:
    Py_XDECREF(ain1);
    Py_XDECREF(ain2);
    Py_XDECREF(aout);
    Py_XDECREF(afill);
    Py_XDECREF(newfill);
    return NULL;
}

/*******************************************************************/

static char doc_order_filterND[] = "out = _order_filterND(a,domain,order)";

static PyObject *sigtools_order_filterND(PyObject *dummy, PyObject *args) {
	PyObject *domain, *a0;
	int order=0;
	
	if (!PyArg_ParseTuple(args, "OO|i", &a0, &domain, &order)) return NULL;
	
	return PyArray_OrderFilterND(a0, domain, order);
}



static char doc_linear_filter[] = "(y,Vf) = _linear_filter(b,a,X,Dim=-1,Vi=None)  implemented using Direct Form II transposed flow diagram. If Vi is not given, Vf is not returned.";
 
static PyObject *sigtools_linear_filter(PyObject *dummy, PyObject *args) {
	PyObject *b=NULL, *a=NULL, *X=NULL, *Vi=NULL;
	PyArrayObject *arY=NULL, *arb=NULL, *ara=NULL, *arX=NULL, *arVi=NULL, *arVf=NULL;
	Generic_Array x, y, *vi=NULL, *vf=NULL;
	Generic_Vector Vb, Va;
	int dim = -1, typenum, thedim;
	char *ara_ptr, input_flag = 0;
	BasicFilterFunction *basic_filter;

	if (!PyArg_ParseTuple(args, "OOO|iO", &b, &a, &X, &dim, &Vi))
	  return NULL;


	typenum = PyArray_ObjectType(b, 0);  
	typenum = PyArray_ObjectType(a, typenum);
	typenum = PyArray_ObjectType(X, typenum);
	if (Vi != NULL) typenum = PyArray_ObjectType(Vi, typenum);

	arY = NULL; arVf = NULL; 
	ara = NULL; arb = NULL; arX = NULL; arVi = NULL;
	ara = (PyArrayObject *)PyArray_ContiguousFromObject(a, typenum, 1, 1);
	arb = (PyArrayObject *)PyArray_ContiguousFromObject(b, typenum, 1, 1);
	arX = (PyArrayObject *)PyArray_FromObject(X, typenum, 0, 0);
	if (ara == NULL || arb == NULL || arX == NULL) goto fail;

	if (dim < -arX->nd || dim > arX->nd - 1) {
	  PyErr_SetString(PyExc_ValueError, 
			  "selected axis is out of range");
	  goto fail;
	}

	if (dim < 0) 
	  thedim = arX->nd + dim;
	else
	  thedim = dim;

	if (Vi != NULL) {
	  arVi = (PyArrayObject *)PyArray_FromObject(Vi, typenum, arX->nd, arX->nd);
	  if (arVi == NULL) goto fail;
	  input_flag = (PyArray_Size((PyObject *)arVi) > 0);
	}

	arY = (PyArrayObject *)PyArray_FromDims(arX->nd, arX->dimensions, typenum);
	if (arY == NULL) goto fail;

	if (input_flag) {
	  arVf = (PyArrayObject *)PyArray_FromDims(arVi->nd, arVi->dimensions, typenum);
	}
	
       	basic_filter = BasicFilterFunctions[(int)(arX->descr->type_num)];
 	if (basic_filter == NULL) {
 		PyErr_SetString(PyExc_ValueError, 
 			"linear_filter not available for this type");
 		goto fail;
 	}

	/* copy header information to generic structures ready 
	    to transfer to non-Python-Specific code */

	Py_copy_info_vec(&Va, ara);
	Py_copy_info_vec(&Vb, arb);
	Py_copy_info(&x, arX);
	Py_copy_info(&y, arY);
	/* Skip over leading zeros in vector representing denominator (a)*/
	ara_ptr = ara->data;
	while(memcmp(ara_ptr,Va.zero,Va.elsize) == 0) {
	  ara_ptr += Va.elsize;
	  Va.data = ara_ptr;
	  Va.numels--;
	}

	if (input_flag) {
	  if (arVi->dimensions[thedim] != (Va.numels > Vb.numels ? Va.numels : Vb.numels) - 1) {
	    PyErr_SetString(PyExc_ValueError,
			    "The number of initial conditions must be max([len(a),len(b)]) - 1");
	    goto fail;
	  }
        vi = (Generic_Array *)malloc(sizeof(Generic_Array));
	vf = (Generic_Array *)malloc(sizeof(Generic_Array));
	Py_copy_info(vi, arVi);
	Py_copy_info(vf, arVf);
	}

	/* If dimension to filter along is negative, make it the
	   correct positive dimension */

        /* fprintf(stderr, "Here.\n"); */

	RawFilter(Vb, Va, x, y, vi, vf, basic_filter, thedim);
        /* fprintf(stderr, "Now, Here.\n");*/

	Py_XDECREF(ara);
	Py_XDECREF(arb);
	Py_XDECREF(arX);
	Py_XDECREF(arVi);

	if (!input_flag) {
	  return PyArray_Return(arY);
	}
	else {
	  free(vi); free(vf);
	  return Py_BuildValue("(NN)",arY,arVf);
 	}

 fail:
	Py_XDECREF(ara);
	Py_XDECREF(arb);
	Py_XDECREF(arX);
	Py_XDECREF(arVi);
	Py_XDECREF(arVf);
	Py_XDECREF(arY);
        return NULL;
}


static char doc_remez[] = "h = _remez(numtaps, bands, des, weight, type, Hz, maxiter, grid_density) \n  returns the optimal (in the Chebyshev/minimax sense) FIR filter impulse \n  response given a set of band edges, the desired response on those bands,\n  and the weight given to the error in those bands.  Bands is a monotonic\n   vector with band edges given in frequency domain where Hz is the sampling\n   frequency.";
 
static PyObject *sigtools_remez(PyObject *dummy, PyObject *args) {
        PyObject *bands, *des, *weight;
        int k, numtaps, numbands, type = BANDPASS, err; 
	PyArrayObject *a_bands=NULL, *a_des=NULL, *a_weight=NULL;
	PyArrayObject *h=NULL;
	int ret_dimens, maxiter = 25, grid_density = 16;
	double oldvalue, *dptr, Hz = 1.0;
        char mystr[255];
       


        if (!PyArg_ParseTuple(args, "iOOO|idii", &numtaps, &bands, &des, &weight, &type, &Hz, &maxiter, &grid_density))
	    return NULL;

	if (type != BANDPASS && type != DIFFERENTIATOR && type != HILBERT) {
	  PyErr_SetString(PyExc_ValueError,
			  "The type must be BANDPASS, DIFFERENTIATOR, or HILBERT.");
	  return NULL;
	}
	
	if (numtaps < 2) {
	  PyErr_SetString(PyExc_ValueError,
			  "The number of taps must be greater than 1.");
	  return NULL;
	}
	
	
	a_bands = (PyArrayObject *)PyArray_ContiguousFromObject(bands, PyArray_DOUBLE,1,1);
	if (a_bands == NULL) goto fail;
	a_des = (PyArrayObject *)PyArray_ContiguousFromObject(des, PyArray_DOUBLE,1,1);
	if (a_des == NULL) goto fail;
	a_weight = (PyArrayObject *)PyArray_ContiguousFromObject(weight, PyArray_DOUBLE,1,1);
	if (a_weight == NULL) goto fail;


	numbands = a_des->dimensions[0];
	if ((a_bands->dimensions[0] != 2*numbands) || (a_weight->dimensions[0] != numbands)) {
	  PyErr_SetString(PyExc_ValueError,
			  "The inputs desired and weight must have same length.\n  The input bands must have twice this length.");
	  goto fail;
	}

      /* Check the bands input to see if it is monotonic, divide by 
	 Hz to take from range 0 to 0.5 and check to see if in that range */ 

	dptr = (double *)a_bands->data;
	oldvalue = 0;
	for (k=0; k < 2*numbands; k++) {
	  if (*dptr < oldvalue) {
	    PyErr_SetString(PyExc_ValueError,
			  "Bands must be monotonic starting at zero.");
	    goto fail;
	  }
	  if (*dptr * 2 > Hz) {
	    PyErr_SetString(PyExc_ValueError,
			  "Band edges should be less than 1/2 the sampling frequency");
	    goto fail;
	  }
	  oldvalue = *dptr;
	  *dptr = oldvalue / Hz;  /* Change so that sampling frequency is 1.0 */
	  dptr++;
	}

	ret_dimens = numtaps;
	h = (PyArrayObject *)PyArray_FromDims(1, &ret_dimens, PyArray_DOUBLE);
	if (h == NULL) goto fail;

	err=pre_remez((double *)h->data, numtaps, numbands, (double *)a_bands->data, (double *)a_des->data, (double *)a_weight->data, type, maxiter, grid_density);
        if (err < 0) {
	  if (err == -1) {
            sprintf(mystr,"Failure to converge after %d iterations.\n      Design may still be correct.",maxiter);
	    PyErr_SetString(PyExc_ValueError, mystr);
	    goto fail;
	  }
	  else if (err == -2) {
	    PyErr_NoMemory();
            goto fail;
	  }
	}

	Py_DECREF(a_bands);
	Py_DECREF(a_des);
	Py_DECREF(a_weight);

	return PyArray_Return(h);

 fail:
	Py_XDECREF(a_bands);
	Py_XDECREF(a_des);
	Py_XDECREF(a_weight);
	Py_XDECREF(h);
	return NULL;
}
   
static char doc_median2d[] = "filt = _median2d(data, size)";

extern void f_medfilt2(float*,float*,int*,int*);
extern void d_medfilt2(double*,double*,int*,int*);
extern void b_medfilt2(unsigned char*,unsigned char*,int*,int*);

static PyObject *sigtools_median2d(PyObject *dummy, PyObject *args)
{
    PyObject *image=NULL, *size=NULL;
    int typenum;
    PyArrayObject *a_image=NULL, *a_size=NULL;
    PyArrayObject *a_out=NULL;
    int Nwin[2] = {3,3};

    if (!PyArg_ParseTuple(args, "O|O", &image, &size)) return NULL;

    typenum = PyArray_ObjectType(image, 0);
    a_image = (PyArrayObject *)PyArray_ContiguousFromObject(image, typenum, 2, 2);
    if (a_image == NULL) goto fail;

    if (size != NULL) {
	a_size = (PyArrayObject *)PyArray_ContiguousFromObject(size, PyArray_LONG, 1, 1);
	if (a_size == NULL) goto fail;
	if ((RANK(a_size) != 1) || (DIMS(a_size)[0] < 2)) 
	    PYERR("Size must be a length two sequence");
	Nwin[0] = ((long *)DATA(a_size))[0];
	Nwin[1] = ((long *)DATA(a_size))[1];
    }  

    a_out = (PyArrayObject *)PyArray_FromDims(2,DIMS(a_image),typenum);
    if (a_out == NULL) goto fail;

    if (setjmp(MALLOC_FAIL)) {
	PYERR("Memory allocation error.");
    }
    else {
	switch (typenum) {
	case PyArray_UBYTE:
	    b_medfilt2((unsigned char *)DATA(a_image), (unsigned char *)DATA(a_out), Nwin, DIMS(a_image));
	    break;
	case PyArray_FLOAT:
	    f_medfilt2((float *)DATA(a_image), (float *)DATA(a_out), Nwin, DIMS(a_image));
	    break;
	case PyArray_DOUBLE:
	    d_medfilt2((double *)DATA(a_image), (double *)DATA(a_out), Nwin, DIMS(a_image));
	    break;
	default:
	  PYERR("2D median filter only supports Int8, Float32, and Float64.");
	}
    }

    Py_DECREF(a_image);
    Py_XDECREF(a_size);

    return PyArray_Return(a_out);
 
 fail:
    Py_XDECREF(a_image);
    Py_XDECREF(a_size);
    Py_XDECREF(a_out);
    return NULL;

}



static struct PyMethodDef toolbox_module_methods[] = {
	{"_correlateND", sigtools_correlateND, METH_VARARGS, doc_correlateND},
	{"_convolve2d", sigtools_convolve2d, METH_VARARGS, doc_convolve2d},
	{"_order_filterND", sigtools_order_filterND, METH_VARARGS, doc_order_filterND},
	{"_linear_filter",sigtools_linear_filter, METH_VARARGS, doc_linear_filter},
	{"_remez",sigtools_remez, METH_VARARGS, doc_remez},
	{"_medfilt2d", sigtools_median2d, METH_VARARGS, doc_median2d},
	{NULL,		NULL, 0}		/* sentinel */
};

/* Initialization function for the module (*must* be called initsigtools) */

DL_EXPORT(void) initsigtools(void) {
        PyObject *m, *d;
	
	/* Create the module and add the functions */
	m = Py_InitModule("sigtools", toolbox_module_methods);

	/* Import the C API function pointers for the Array Object*/
	import_array();

	/* Make sure the multiarraymodule is loaded so that the zero
	   and one objects are defined */
	PyImport_ImportModule("multiarray");
	/* { PyObject *multi = PyImport_ImportModule("multiarray"); } */

	/* Add some symbolic constants to the module */
	d = PyModule_GetDict(m);

	/* PyDict_SetItemString(d,"BANDPASS", PyInt_FromLong((long) BANDPASS));
        PyDict_SetItemString(d,"DIFFERENTIATOR", PyInt_FromLong((long) DIFFERENTIATOR));
        PyDict_SetItemString(d,"HILBERT", PyInt_FromLong((long) HILBERT));
        */


	/* Check for errors */
	if (PyErr_Occurred())
		Py_FatalError("can't initialize module array");
}
