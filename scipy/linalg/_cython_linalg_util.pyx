from ._cython_linalg_util cimport lapack_t

# ====================== swap_c_and_f_layout : s, d, c, z ====================
@cython.cdivision(True)
@cython.wraparound(False)
@cython.boundscheck(False)
@cython.initializedcheck(False)
cdef void swap_c_and_f_layout(lapack_t *a, lapack_t *b, int r, int c, int n) nogil:
    """Recursive matrix transposition for square arrays"""
    cdef int i, j, ith_row, r2, c2
    cdef lapack_t *bb=b
    cdef lapack_t *aa=a
    if r < 16:
        for j in range(c):
            ith_row = 0
            for i in range(r):
            # Basically b[i*n+j] = a[j*n+i] without index math
                bb[ith_row] = aa[i]
                ith_row += n
            aa += n
            bb += 1
    else:
        # If tall
        if (r > c):
            r2 = r//2
            swap_c_and_f_layout(a, b, r2, c, n)
            swap_c_and_f_layout(a + r2, b+(r2)*n, r-r2, c, n)
        else:  # Nope
            c2 = c//2
            swap_c_and_f_layout(a, b, r, c2, n);
            swap_c_and_f_layout(a+(c2)*n, b+c2, r, c-c2, n)
# ============================================================================
