# Authors: Travis Oliphant, Matthew Brett

"""
Module for reading and writing matlab (TM) .mat files
"""

import os
import sys

from numpy import *

from mio4 import MatFile4Reader, MatFile4Writer
from mio5 import MatFile5Reader, MatFile5Writer


def mat_reader_factory(file_name, appendmat=True, **kwargs):
    """Create reader for matlab (TM) .mat format files
    
    See docstring for loadmat for input options
    """
    if appendmat and file_name[-4:] == ".mat":
        file_name = file_name[:-4]
    if os.sep in file_name:
        full_file_name = file_name
        if appendmat:
            full_name = file_name + ".mat"
    else:
        full_name = None
        junk,file_name = os.path.split(file_name)
        for path in sys.path:
            test_name = os.path.join(path,file_name)
            if appendmat:
                test_name += ".mat"
            try:
                fid = open(test_name,'rb')
                fid.close()
                full_name = test_name
                break
            except IOError:
                pass
        if full_name is None:
            raise IOError, "%s not found on the path." % file_name

    byte_stream = open(full_name, 'rb')
    MR = MatFile4Reader(byte_stream, **kwargs)
    if MR.format_looks_right():
        return MR
    return MatFile5Reader(byte_stream, **kwargs)

def loadmat(file_name,  mdict=None, appendmat=True, basename='raw', **kwargs):
    ''' Load Matlab(tm) file

    file_name          - Name of the mat file
                         (do not need .mat extension if appendmat==True)
                         If name not a full path name, search for the file on
                         the sys.path list and use the first one found (the
                         current directory is searched first).
    m_dict             - optional dictionary in which to insert matfile variables 
    appendmat          - True to append the .mat extension to the end of the
                         given filename.
    base_name          - base name for unnamed variables (unused in code)
    byte_order         - byte order ('native', 'little', 'BIG')
                          in ('native', '=')
                          or in ('little', '<')
                          or in ('BIG', '>')
    mat_dtype          - return arrays in same dtype as loaded into matlab
                          (instead of the dtype with which they are saved)
    squeeze_me         - whether to squeeze matrix dimensions or not
    chars_as_strings   - whether to convert char arrays to string arrays
    mat_dtype          - return matrices with datatype that matlab would load as
                          (rather than in the datatype matlab saves as)
    matlab_compatible   - returns matrices as would be loaded by matlab
                          (implies squeeze_me=False, chars_as_strings=False,
                          mat_dtype=True)

    v4 (Level 1.0), v6 and v7.1 matfiles are supported.  

    '''
    MR = mat_reader_factory(file_name, appendmat, **kwargs)
    matfile_dict = MR.get_variables()
    if mdict is not None:
        mdict.update(matfile_dict)
    else:
        mdict = matfile_dict
    return mdict

def savemat(file_name, mdict, appendmat=True):
    """Save a dictionary of names and arrays into the MATLAB-style .mat file.

    This saves the arrayobjects in the given dictionary to a matlab
    Version 4 style .mat file.
    
    @appendmat  - if true, appends '.mat' extension to filename, if not present
    """
    if appendmat and file_name[-4:] != ".mat":
        file_name = file_name + ".mat"
    file_stream = open(file_name, 'wb')
    MW = MatFile4Writer(file_stream)
    MW.put_variables(mdict)
    file_stream.close()
    
if __name__ == '__main__':
    D = savemat('test.mat', {'a': 1})
