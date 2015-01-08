from scipy.linalg cimport lapack_pointers as lapack

def dlamch(cmach):
    cmach_bytes = bytes(cmach)
    cdef char* cmach_char = cmach_bytes
    return lapack.dlamch(cmach_char)

def slamch(cmach):
    cmach_bytes = bytes(cmach)
    cdef char* cmach_char = cmach_bytes
    return lapack.slamch(cmach_char)
