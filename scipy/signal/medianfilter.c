
/*--------------------------------------------------------------------*/

#include "Python.h"
#define NO_IMPORT_ARRAY
#include "numpy/noprefix.h"


/* defined below */
void f_medfilt2(float*,float*,intp*,intp*);
void d_medfilt2(double*,double*,intp*,intp*);
void b_medfilt2(unsigned char*,unsigned char*,intp*,intp*);
extern char *check_malloc (int);


/* The QUICK_SELECT routine is based on Hoare's Quickselect algorithm,
 * with unrolled recursion. 
 * Author: Thouis R. Jones, 2008 
 */

#define ELEM_SWAP(t, a, x, y) {register t temp = (a)[x]; (a)[x] = (a)[y]; (a)[y] = temp;}
#define FIRST_LOWEST(x, y, z) (((x) < (y)) && ((x) < (z)))
#define FIRST_HIGHEST(x, y, z) (((x) > (y)) && ((x) > (z)))
#define LOWEST_IDX(a, x, y) (((a)[x] < (a)[y]) ? (x) : (y))
#define HIGHEST_IDX(a, x, y) (((a)[x] > (a)[y]) ? (x) : (y))

/* if (l is index of lowest) {return lower of mid,hi} else if (l is index of highest) {return higher of mid,hi} else return l */
#define MEDIAN_IDX(a, l, m, h) (FIRST_LOWEST((a)[l], (a)[m], (a)[h]) ? LOWEST_IDX(a, m, h) : (FIRST_HIGHEST((a)[l], (a)[m], (a)[h]) ? HIGHEST_IDX(a, m, h) : (l)))

#define QUICK_SELECT(NAME, TYPE)                                        \
TYPE NAME(TYPE arr[], int n)                                            \
{                                                                       \
    int lo, hi, mid, md;                                                \
    int median_idx;                                                     \
    int ll, hh;                                                         \
    TYPE piv;                                                           \
                                                                        \
    lo = 0; hi = n-1;                                                   \
    median_idx = (n - 1) / 2; /* lower of middle values for even-length arrays */ \
                                                                        \
    while (1) {                                                         \
        if ((hi - lo) < 2) {                                            \
            if (arr[hi] < arr[lo]) ELEM_SWAP(TYPE, arr, lo, hi);        \
            return arr[median_idx];                                     \
        }                                                               \
                                                                        \
        mid = (hi + lo) / 2;                                            \
        /* put the median of lo,mid,hi at position lo - this will be the pivot */ \
        md = MEDIAN_IDX(arr, lo, mid, hi);                              \
        ELEM_SWAP(TYPE, arr, lo, md);                                   \
                                                                        \
        /* Nibble from each end towards middle, swapping misordered items */ \
        piv = arr[lo];                                                  \
        for (ll = lo+1, hh = hi;; ll++, hh--) {                         \
	    while (arr[ll] < piv) ll++;					\
	    while (arr[hh] > piv) hh--;					\
	    if (hh < ll) break;						\
	    ELEM_SWAP(TYPE, arr, ll, hh);				\
        }                                                               \
        /* move pivot to top of lower partition */                      \
        ELEM_SWAP(TYPE, arr, hh, lo);                                   \
        /* set lo, hi for new range to search */                        \
        if (hh < median_idx) /* search upper partition */               \
            lo = hh+1;                                                  \
        else if (hh > median_idx) /* search lower partition */          \
            hi = hh-1;                                                  \
        else                                                            \
            return piv;                                                 \
    }                                                                   \
}


/* 2-D median filter with zero-padding on edges. */
#define MEDIAN_FILTER_2D(NAME, TYPE, SELECT)                            \
void NAME(TYPE* in, TYPE* out, intp* Nwin, intp* Ns)                    \
{                                                                       \
    int nx, ny, hN[2];                                                  \
    int pre_x, pre_y, pos_x, pos_y;                                     \
    int subx, suby, k, totN;                                            \
    TYPE *myvals, *fptr1, *fptr2, *ptr1, *ptr2;                         \
                                                                        \
    totN = Nwin[0] * Nwin[1];                                           \
    myvals = (TYPE *) check_malloc( totN * sizeof(TYPE));               \
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
    free(myvals);                                                       \
}


/* define quick_select for floats, doubles, and unsigned characters */
QUICK_SELECT(f_quick_select, float)
QUICK_SELECT(d_quick_select, double)
QUICK_SELECT(b_quick_select, unsigned char)

/* define medfilt for floats, doubles, and unsigned characters */
MEDIAN_FILTER_2D(f_medfilt2, float, f_quick_select)
MEDIAN_FILTER_2D(d_medfilt2, double, d_quick_select)
MEDIAN_FILTER_2D(b_medfilt2, unsigned char, b_quick_select)
