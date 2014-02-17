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


cpdef lil_get1(cnp.intp_t M, cnp.intp_t N, object[:] rows, object[:] datas, cnp.intp_t i, cnp.intp_t j):
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


cdef lil_insert_nocheck(cnp.intp_t N, list row, list data, cnp.intp_t j, object x):
    cdef cnp.intp_t pos

    pos = bisect_left(row, j)
    if pos == len(row):
        row.append(j)
        data.append(x)
    elif row[pos] != j:
        row.insert(pos, j)
        row.insert(pos, x)
    else:
        data[pos] = x


def lil_fancy_get(cnp.intp_t M, cnp.intp_t N,
                  object[:] rows,
                  object[:] data,
                  object[:] new_rows,
                  object[:] new_data,
                  cnp.ndarray[cnp.intp_t, ndim=2] i_idx,
                  cnp.ndarray[cnp.intp_t, ndim=2] j_idx):
    cdef cnp.intp_t x, y, i, j
    cdef object value

    for x in range(i_idx.shape[0]):
        for y in range(i_idx.shape[1]):
            i = i_idx[x,y]
            j = j_idx[x,y]

            value = lil_get1(M, N, rows, data, i, j)

            if value is 0:
                # Object identity as shortcut
                continue

            lil_insert_nocheck(i_idx.shape[1],
                               new_rows[x], new_data[x],
                               y, value)


cdef bisect_left(list a, cnp.intp_t x):
    cdef cnp.intp_t hi = len(a)
    cdef cnp.intp_t lo = 0
    cdef cnp.intp_t mid_v
    cdef cpython.PyObject* v

    while lo < hi:
        mid = (lo + hi)//2
        v = cpython.list.PyList_GET_ITEM(a, mid)
        mid_v = cpython.int.PyInt_AsLong(<object>v)
        if mid_v < x:
            lo = mid + 1
        else:
            hi = mid
    return lo
