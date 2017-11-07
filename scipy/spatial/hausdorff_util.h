// Utility C code for directed_hausdorff
// Try to improve on the performance of the original Cython directed_hausdorff

#include <stdlib.h>
#include <stdbool.h>
#include <math.h>
#include <stddef.h>
#include <stdio.h>

struct return_values 
{
    double cmax;
    int    index_1;
    int    index_2;
};

// function for inner loop of directed_hausdorff algorithm
void hausdorff_loop(const int data_dims,
                    double ar1[],
                    double ar2[],
                    int N1,
                    int N2,
                    struct return_values * ret_vals)
{
    double               d, cmin, diff;
    bool                 no_break_happened;
    double * const ar1_start_Ptr = ar1;
    double * const ar2_start_Ptr = ar2;
    int size_ar1 = data_dims * N1;
    int size_ar2 = data_dims * N2;
    const double * const ar1_end_Ptr = &ar1[size_ar1 - 1];
    const double * const ar2_end_Ptr = &ar2[size_ar2 - 1];

    ret_vals->cmax = 0;
    
    while (ar1 <= ar1_end_Ptr) {
        no_break_happened = 1;
        cmin = INFINITY;
        ar2 = ar2_start_Ptr;
        while (ar2 <= ar2_end_Ptr) {
            d = 0;
            for ( int k = 0; k < data_dims; ++k, ++ar1, ++ar2) {
                diff = *ar1 - *ar2;
                d += (diff * diff);
            }
            if (d < (ret_vals->cmax)) {
                --no_break_happened;
                ar1 -= data_dims;
                break;
            }

            if (d < cmin)
                cmin = d;
        ar1 -= data_dims;
        }
        if ( (cmin >= (ret_vals->cmax)) && (no_break_happened)) {
            ret_vals->cmax = cmin;
            ret_vals->index_1 = (ar1 + data_dims - ar1_start_Ptr + 1) / data_dims - 1;
            ret_vals->index_2 = (ar2 - ar2_start_Ptr + 1) / data_dims - 1;
        }
    ar1 += data_dims;
    }
}
