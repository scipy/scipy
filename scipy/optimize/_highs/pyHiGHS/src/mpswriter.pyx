# distutils: language=c++
# cython: language_level=3

from libc.stdio cimport FILE, fopen, fclose
from libcpp cimport bool
from libcpp.vector cimport vector
from libcpp.string cimport string

from HighsStatus cimport HighsStatus

from scipy.sparse import csc_matrix

cdef extern from "HMPSIO.h" nogil:
    HighsStatus writeMPS(FILE* logfile, const char* filename, const int& numRow, const int& numCol,
                         const int& numInt, const int& objSense, const double& objOffset,
                         const vector[int]& Astart, const vector[int]& Aindex,
                         const vector[double]& Avalue, const vector[double]& colCost,
                         const vector[double]& colLower, const vector[double]& colUpper,
                         const vector[double]& rowLower, const vector[double]& rowUpper,
                         const vector[int]& integerColumn, const vector[string]& col_names,
                         const vector[string]& row_names, const bool use_free_format)

def mpswriter(
        string filename,
        const double[::1] c,
        A,
        const double[::1] rhs,
        const double[::1] lhs,
        const double[::1] lb,
        const double[::1] ub,
        const int[::1] integer_valued,
        const bool use_free_format=True):

    cdef int ii = 0
    A = csc_matrix(A)

    cdef int numRow = A.shape[0]
    cdef int numCol = A.shape[1]
    cdef int numInt = A.nnz
    cdef int objSense = 1 # MIN for now
    cdef double objOffset = 0 # This is a RHS on cost row

    cdef vector[int] Astart = A.indptr
    cdef vector[int] Aindex = A.indices
    cdef vector[double] Avalue = A.data

    cdef vector[double] colCost
    for ii in range(numCol):
        colCost.push_back(c[ii])

    cdef vector[double] colLower
    cdef vector[double] colUpper
    cdef vector[double] rowLower
    cdef vector[double] rowUpper
    for ii in range(numCol):
        colLower.push_back(lb[ii])
        colUpper.push_back(ub[ii])
    for ii in range(numRow):
        rowLower.push_back(lhs[ii])
        rowUpper.push_back(rhs[ii])

    cdef vector[int] integerColumn
    integerColumn.push_back(0)
    integerColumn.clear()
    for ii in range(integer_valued.size):
        integerColumn.push_back(integer_valued[ii])

    cdef vector[string] row_names
    cdef vector[string] col_names
    cdef int jj = 0
    for ii in range(numRow):
        row_names.push_back(b'row%d' % ii)
    for jj in range(numCol):
        col_names.push_back(b'col%d' % jj)

    cdef FILE * logfile = fopen(filename.c_str(), 'w')

    writeMPS(
        logfile, filename.c_str(), numRow, numCol,
        numInt, objSense, objOffset,
        Astart, Aindex,
        Avalue, colCost,
        colLower, colUpper,
        rowLower, rowUpper,
        integerColumn, col_names,
        row_names, use_free_format)

    fclose(logfile)
