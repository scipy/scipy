# -*- python -*- like file
''' Utilities for generic processing of return arrays from read
'''

import numpy as np
cimport numpy as cnp
from mio_utils cimport FileReadOpts


cpdef size_t cproduct(tup):
    cdef size_t res = 1
    cdef int i
    for i in range(len(tup)):
        res *= tup[i]
    return res
           

cdef class FileReadOpts:
    def __new__(self,
                chars_as_strings,
                mat_dtype,
                squeeze_me):
        self.chars_as_strings = chars_as_strings
        self.mat_dtype = mat_dtype
        self.squeeze_me = squeeze_me
        

cpdef cnp.ndarray process_element(cnp.ndarray arr,
                    FileReadOpts file_read_opts,
                    object mat_dtype):
    cdef:
        int i
    if file_read_opts.chars_as_strings and arr.dtype.kind == 'U':
        arr = chars_to_strings(arr)
    if file_read_opts.mat_dtype:
        # Apply options to replicate matlab datatype on load
        if mat_dtype is not None:
            arr = arr.astype(mat_dtype)
    if file_read_opts.squeeze_me:
        if not arr.size:
            arr = np.array([])
        else:
            arr = np.squeeze(arr)
            if not arr.shape and arr.dtype.isbuiltin: # 0d coverted to scalar
                arr = arr.item()
    return arr


cpdef cnp.ndarray chars_to_strings(in_arr):
    ''' Convert final axis of char array to strings

    Parameters
    ----------
    arr : array
       dtype of 'U1'
       
    Returns
    -------
    str_arr : array
       dtype of 'UN' where N is the length of the last dimension of
       ``arr``
    '''
    # make numpy version of array.  Strangely, if we don't do this, the
    # array can change shape (1,1) to (1,) for example
    cdef cnp.ndarray arr = in_arr
    cdef int ndim = arr.ndim
    cdef cnp.npy_intp *dims = arr.shape
    cdef cnp.npy_intp last_dim = dims[ndim-1]
    cdef object new_dt_str
    if last_dim == 0: # deal with empty array case
        new_dt_str = arr.dtype.str
    else: # make new dtype string with N appended
        new_dt_str = arr.dtype.str[:-1] + str(last_dim)
    # Copy to deal with F ordered arrays
    arr = np.ascontiguousarray(arr)
    arr = arr.view(new_dt_str)
    return arr.reshape(in_arr.shape[:-1])
