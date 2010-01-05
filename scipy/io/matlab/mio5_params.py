''' Constants and classes for matlab 5 read and write

See also mio5_utils.pyx where these same constants arise as c enums.

If you make changes in this file, don't forget to change mio5_utils.pyx
'''

import numpy as np


miINT8 = 1
miUINT8 = 2
miINT16 = 3
miUINT16 = 4
miINT32 = 5
miUINT32 = 6
miSINGLE = 7
miDOUBLE = 9
miINT64 = 12
miUINT64 = 13
miMATRIX = 14
miCOMPRESSED = 15
miUTF8 = 16
miUTF16 = 17
miUTF32 = 18

mxCELL_CLASS = 1
mxSTRUCT_CLASS = 2
# The March 2008 edition of "Matlab 7 MAT-File Format" says that
# mxOBJECT_CLASS = 3, whereas matrix.h says that mxLOGICAL = 3.
# Matlab 2008a appears to save logicals as type 9, so we assume that
# the document is correct.  See type 18, below.
mxOBJECT_CLASS = 3
mxCHAR_CLASS = 4
mxSPARSE_CLASS = 5
mxDOUBLE_CLASS = 6
mxSINGLE_CLASS = 7
mxINT8_CLASS = 8
mxUINT8_CLASS = 9
mxINT16_CLASS = 10
mxUINT16_CLASS = 11
mxINT32_CLASS = 12
mxUINT32_CLASS = 13
# The following are not in the March 2008 edition of "Matlab 7
# MAT-File Format," but were guessed from matrix.h.
mxINT64_CLASS = 14
mxUINT64_CLASS = 15
mxFUNCTION_CLASS = 16
# Not doing anything with these at the moment.
mxOPAQUE_CLASS = 17 # This appears to be a function workspace
# https://www-old.cae.wisc.edu/pipermail/octave-maintainers/2007-May/002824.html
mxOBJECT_CLASS_FROM_MATRIX_H = 18


class mat_struct(object):
    ''' Placeholder for holding read data from structs

    We deprecate this method of holding struct information, and will
    soon remove it, in favor of the recarray method (see loadmat
    docstring)
    '''
    pass


class MatlabObject(np.ndarray):
    ''' ndarray Subclass to contain matlab object '''
    def __new__(cls, input_array, classname=None):
        # Input array is an already formed ndarray instance
        # We first cast to be our class type
        obj = np.asarray(input_array).view(cls)
        # add the new attribute to the created instance
        obj.classname = classname
        # Finally, we must return the newly created object:
        return obj

    def __array_finalize__(self,obj):
        # reset the attribute from passed original object
        self.classname = getattr(obj, 'classname', None)
        # We do not need to return anything


class MatlabFunction(np.ndarray):
    ''' Subclass to signal this is a matlab function '''
    def __new__(cls, input_array):
        obj = np.asarray(input_array).view(cls)
        return obj


class MatlabOpaque(np.ndarray):
    ''' Subclass to signal this is a matlab opaque matrix '''
    def __new__(cls, input_array):
        obj = np.asarray(input_array).view(cls)
        return obj


OPAQUE_DTYPE = np.dtype(
    [('s0', 'O'), ('s1', 'O'), ('s2', 'O'), ('arr', 'O')])
