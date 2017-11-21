// Utility C code for directed_hausdorff
// Try to improve on the performance of the original Cython directed_hausdorff

#include <math.h>

struct return_values 
{
    double cmax;
    int    index_1;
    int    index_2;
};

// function for inner loop of directed_hausdorff algorithm
void hausdorff_loop(const int data_dims,
                    const double ar1[],
                    const double ar2[],
                    const int N1,
                    const int N2,
                    struct return_values * ret_vals)
{
    double               d, cmin, diff;
    const double * const ar1_start_Ptr = ar1;
    const double * const ar2_start_Ptr = ar2;
    const int size_ar1 = data_dims * N1;
    const int size_ar2 = data_dims * N2;
    const double * const ar1_end_Ptr = &ar1[size_ar1 - 1];
    const double * const ar2_end_Ptr = &ar2[size_ar2 - 1];
    int k;

    ret_vals->cmax = 0;
    
    for ( ; ar1 <= ar1_end_Ptr; ar1 += data_dims ) {
        cmin = INFINITY;
        for ( ar2 = ar2_start_Ptr; ar2 <= ar2_end_Ptr; ar2 += data_dims) {
            d = 0;
            for ( k = 0; k < data_dims; ++k ) {
                diff = *(ar1 + k) - *(ar2 + k);
                d += (diff * diff);
            }
            if (d < ret_vals->cmax)
                goto main_loop_continue;
            else if (d < cmin)
                cmin = d;
        }

        if ( cmin >= ret_vals->cmax) {
            ret_vals->cmax = cmin;
            ret_vals->index_1 = ar1 - ar1_start_Ptr;
            ret_vals->index_2 = ar2 - ar2_start_Ptr;
        }

    main_loop_continue: ;
    }
    // minimize processing of hausdorff pair
    // indices by placing operations outside
    // looping logic
    ret_vals->index_1 = (++ret_vals->index_1) / data_dims;
    ret_vals->index_2 = ret_vals->index_2 / data_dims - 1;
}
