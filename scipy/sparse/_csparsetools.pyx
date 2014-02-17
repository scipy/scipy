"""
Fast snippets for sparse matrices
"""

cimport cython
cimport cpython.list
cimport cpython.int
cimport cpython
cimport numpy as cnp
import numpy as np


ctypedef fused idx_t:
    cnp.int32_t
    cnp.int64_t


ctypedef fused value_t:
    cnp.npy_bool
    cnp.npy_int8
    cnp.npy_uint8
    cnp.npy_int16
    cnp.npy_uint16
    cnp.npy_int32
    cnp.npy_uint32
    cnp.npy_int64
    cnp.npy_uint64
    cnp.npy_float32
    cnp.npy_float64
    float complex
    double complex


cpdef lil_get1(cnp.npy_intp M, cnp.npy_intp N, object[:] rows, object[:] datas, cnp.npy_intp i, cnp.npy_intp j):
    cdef list row, data

    if i < 0:
        i += M
    if i < 0 or i >= M:
        raise IndexError('row index out of bounds')

    if j < 0:
        j += N
    if j < 0 or j >= N:
        raise IndexError('column index out of bounds')

    row = rows[i]
    data = datas[i]
    pos = bisect_left(row, j)

    if pos != len(data) and row[pos] == j:
        return data[pos]
    else:
        return 0


cpdef lil_insert(cnp.npy_intp M, cnp.npy_intp N, object[:] rows, object[:] datas,
                 cnp.npy_intp i, cnp.npy_intp j, value_t x):
    cdef list row, data
    cdef int is_zero

    if i < 0:
        i += M
    if i < 0 or i >= M:
        raise IndexError('row index out of bounds')

    if j < 0:
        j += N
    if j < 0 or j >= N:
        raise IndexError('column index out of bounds')

    row = rows[i]
    data = datas[i]

    if x == 0:
        lil_deleteat_nocheck(rows[i], datas[i], j)
    else:
        lil_insertat_nocheck(rows[i], datas[i], j, x)


cdef lil_insertat_nocheck(list row, list data, cnp.npy_intp j, object x):
    cdef cnp.npy_intp pos

    pos = bisect_left(row, j)
    if pos == len(row):
        row.append(j)
        data.append(x)
    elif row[pos] != j:
        row.insert(pos, j)
        data.insert(pos, x)
    else:
        data[pos] = x


cdef lil_deleteat_nocheck(list row, list data, cnp.npy_intp j):
    cdef cnp.npy_intp pos
    pos = bisect_left(row, j)
    if pos < len(row) and row[pos] == j:
        del row[pos]
        del data[pos]


def lil_fancy_get(cnp.npy_intp M, cnp.npy_intp N,
                  object[:] rows,
                  object[:] data,
                  object[:] new_rows,
                  object[:] new_data,
                  cnp.ndarray[cnp.npy_intp, ndim=2] i_idx,
                  cnp.ndarray[cnp.npy_intp, ndim=2] j_idx):
    cdef cnp.npy_intp x, y, i, j
    cdef object value

    for x in range(i_idx.shape[0]):
        for y in range(i_idx.shape[1]):
            i = i_idx[x,y]
            j = j_idx[x,y]

            value = lil_get1(M, N, rows, data, i, j)

            if value is 0:
                # Object identity as shortcut
                continue

            lil_insertat_nocheck(new_rows[x], new_data[x],
                                 y, value)


def lil_fancy_set(cnp.npy_intp M, cnp.npy_intp N,
                  object[:] rows,
                  object[:] data,
                  cnp.ndarray[cnp.npy_intp, ndim=2] i_idx,
                  cnp.ndarray[cnp.npy_intp, ndim=2] j_idx,
                  value_t[:,:] values):

    cdef cnp.npy_intp x, y, i, j

    for x in range(i_idx.shape[0]):
        for y in range(i_idx.shape[1]):
            i = i_idx[x,y]
            j = j_idx[x,y]
            lil_insert[value_t](M, N, rows, data, i, j, values[x, y])


@cython.cdivision(True)
@cython.boundscheck(False)
@cython.wraparound(False)
cdef bisect_left(list a, cnp.npy_intp x):
    cdef cnp.npy_intp hi = len(a)
    cdef cnp.npy_intp lo = 0
    cdef cnp.npy_intp mid, v

    while lo < hi:
        mid = (lo + hi)//2
        v = a[mid]
        if v < x:
            lo = mid + 1
        else:
            hi = mid
    return lo
