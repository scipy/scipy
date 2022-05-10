
ctypedef fused lapack_t:
    float
    double
    (float complex)
    (double complex)

ctypedef fused lapack_cz_t:
    (float complex)
    (double complex)

ctypedef fused lapack_sd_t:
    float
    double

# ====================== swap_c_and_f_layout : s, d, c, z ====================
cdef void swap_c_and_f_layout(lapack_t *a, lapack_t *b, int r, int c, int n) nogil
