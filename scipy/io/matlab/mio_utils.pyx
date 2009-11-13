# -*- python -*- like file
''' Utilities for generic processing of return arrays from read
'''

import numpy as np
cimport numpy as cnp
from mio_utils cimport FileReadOpts


def small_product(tup):
    return cproduct(tup)


cdef size_t cproduct(tup):
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
        

def process_element(cnp.ndarray arr,
                    FileReadOpts file_read_opts,
                    object mat_dtype):
    cdef:
        int i
    if file_read_opts.chars_as_strings and arr.dtype.kind == 'U':
        arr = cchars_to_strings(arr)
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


def chars_to_strings(arr):
    ''' Convert final axis of char array to strings

    Python version, for testing

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
    return cchars_to_strings(arr)


cdef cchars_to_strings(cnp.ndarray arr):
    # Convert char array to string or array of strings
    arr_np = arr
    dims = arr_np.shape
    last_dim = dims[-1]
    old_dt_str = arr.dtype.str[:-1]
    new_dt_str = old_dt_str + str(last_dim)
    arr = arr.view(new_dt_str).reshape(dims[:-1])
    return arr
