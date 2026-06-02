# -*- cython -*-
"""
List of sets of integers, low-level C implementation

Works similarly as

    setlist = [set() for j in range(n)]

but with integer values.

"""

cimport libc.stdlib
cimport numpy as np
import numpy as np

cdef struct setlist_t:
    size_t n
    size_t *sizes
    size_t *alloc_sizes
    int **sets

cdef inline int init(setlist_t *setlist, size_t n, size_t size_guess) except -1:
    """
    Initialise a list of `n` sets with a given guessed size
    """
    cdef int j

    setlist.n = n

    setlist.sets = <int**>libc.stdlib.malloc(sizeof(int*) * n)
    if setlist.sets == NULL:
        raise MemoryError("Failed to allocate memory in setlist.init()")

    setlist.sizes = <size_t*>libc.stdlib.malloc(sizeof(size_t) * n)
    if setlist.sizes == NULL:
        libc.stdlib.free(setlist.sets)
        raise MemoryError("Failed to allocate memory in setlist.init()")

    setlist.alloc_sizes = <size_t*>libc.stdlib.malloc(sizeof(size_t) * n)
    if setlist.alloc_sizes == NULL:
        libc.stdlib.free(setlist.sets)
        libc.stdlib.free(setlist.sizes)
        raise MemoryError("Failed to allocate memory in setlist.init()")

    for j in range(n):
        setlist.sizes[j] = 0
        setlist.alloc_sizes[j] = size_guess
        setlist.sets[j] = <int*>libc.stdlib.malloc(sizeof(int) * size_guess)
        if setlist.sets[j] == NULL:
            for i in range(j):
                libc.stdlib.free(setlist.sets[i])
            libc.stdlib.free(setlist.sets)
            libc.stdlib.free(setlist.sizes)
            libc.stdlib.free(setlist.alloc_sizes)
            raise MemoryError("Failed to allocate memory in setlist.init()")

    return 0

cdef inline void free(setlist_t *setlist) noexcept:
    """
    Free the set list
    """

    cdef int j
    for j in range(setlist.n):
        libc.stdlib.free(setlist.sets[j])
    libc.stdlib.free(setlist.sets)
    libc.stdlib.free(setlist.sizes)
    libc.stdlib.free(setlist.alloc_sizes)
    setlist.sets = NULL
    setlist.sizes = NULL
    setlist.alloc_sizes = NULL
    setlist.n = 0

cdef inline int add(setlist_t *setlist, int n, int value) noexcept nogil:
    """
    Add a value to set `n`
    """

    cdef size_t i, sz
    cdef int *p

    if n < 0 or n >= setlist.n:
        return 1

    for i in range(setlist.sizes[n]):
        if setlist.sets[n][i] == value:
            return 0

    if setlist.sizes[n] >= setlist.alloc_sizes[n]:
        sz = 2*setlist.alloc_sizes[n] + 1
        p = <int*>libc.stdlib.realloc(<void*>setlist.sets[n], sz * sizeof(int))
        if p == NULL:
            return -1
        setlist.sets[n] = p
        setlist.alloc_sizes[n] = sz

    setlist.sets[n][setlist.sizes[n]] = value
    setlist.sizes[n] += 1

    return 0

cdef inline object tocsr(setlist_t *setlist):
    """
    Convert list of sets to CSR format

    Integers for set `i` reside in data[indptr[i]:indptr[i+1]]

    Returns
    -------
    indptr
        CSR indptr
    data
        CSR data

    """
    cdef size_t i, j, pos
    cdef size_t total_size
    cdef np.ndarray[np.npy_int, ndim=1] indptr, data

    total_size = 0
    for j in range(setlist.n):
        total_size += setlist.sizes[j]

    indptr = np.empty((setlist.n+1,), dtype=np.intc)
    data = np.empty((total_size,), dtype=np.intc)

    pos = 0
    for i in range(setlist.n):
        indptr[i] = pos
        for j in range(setlist.sizes[i]):
            data[pos] = setlist.sets[i][j]
            pos += 1
    indptr[setlist.n] = pos

    return indptr, data
