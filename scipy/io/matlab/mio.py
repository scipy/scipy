# Authors: Travis Oliphant, Matthew Brett

"""
Module for reading and writing matlab (TM) .mat files
"""

import os
import sys
import warnings

from miobase import get_matfile_version
from mio4 import MatFile4Reader, MatFile4Writer
from mio5 import MatFile5Reader, MatFile5Writer

__all__ = ['find_mat_file', 'mat_reader_factory', 'loadmat', 'savemat']

def find_mat_file(file_name, appendmat=True):
    ''' Try to find .mat file on system path

    file_name     - file name string
    append_mat    - If True, and file_name does not end in '.mat', appends it
    '''
    if appendmat and file_name[-4:] == ".mat":
        file_name = file_name[:-4]
    if os.sep in file_name:
        full_name = file_name
        if appendmat:
            full_name = file_name + ".mat"
    else:
        full_name = None
        junk, file_name = os.path.split(file_name)
        for path in [os.curdir] + list(sys.path):
            test_name = os.path.join(path, file_name)
            if appendmat:
                test_name += ".mat"
            try:
                fid = open(test_name,'rb')
                fid.close()
                full_name = test_name
                break
            except IOError:
                pass
    return full_name

def mat_reader_factory(file_name, appendmat=True, **kwargs):
    """Create reader for matlab (TM) .mat format files

    See docstring for loadmat for input options
    """
    if isinstance(file_name, basestring):
        full_name = find_mat_file(file_name, appendmat)
        if full_name is None:
            raise IOError, "%s not found on the path." % file_name
        byte_stream = open(full_name, 'rb')
    else:
        try:
            file_name.read(0)
        except AttributeError:
            raise IOError, 'Reader needs file name or open file-like object'
        byte_stream = file_name

    mjv, mnv = get_matfile_version(byte_stream)
    if mjv == 0:
        return MatFile4Reader(byte_stream, **kwargs)
    elif mjv == 1:
        return MatFile5Reader(byte_stream, **kwargs)
    elif mjv == 2:
        raise NotImplementedError('Please use PyTables for matlab v7.3 (HDF) files')
    else:
        raise TypeError('Did not recognize version %s' % mv)

def loadmat(file_name,  mdict=None, appendmat=True, basename='raw', **kwargs):
    ''' Load Matlab(tm) file

    file_name          - Name of the mat file
                         (do not need .mat extension if appendmat==True)
                         If name not a full path name, search for the file on
                         the sys.path list and use the first one found (the
                         current directory is searched first).
                         Can also pass open file-like object
    m_dict             - optional dictionary in which to insert matfile variables
    appendmat          - True to append the .mat extension to the end of the
                         given filename, if not already present
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
    struct_as_record    - whether to load matlab structs as numpy record arrays, or
                          as old-style numpy arrays with dtype=object.
                          (warns if not set, and defaults to False.  non-recarrays
                          cannot be exported via savemat.)

    v4 (Level 1.0), v6 and v7.1 matfiles are supported.

    '''
    if not kwargs.get('struct_as_record', False):
        warnings.warn("loading matlab structures as arrays of dtype=object is deprecated",
                      DeprecationWarning, stacklevel=2)
    MR = mat_reader_factory(file_name, appendmat, **kwargs)
    matfile_dict = MR.get_variables()
    if mdict is not None:
        mdict.update(matfile_dict)
    else:
        mdict = matfile_dict
    return mdict

def savemat(file_name, mdict, appendmat=True, format='4'):
    """Save a dictionary of names and arrays into the MATLAB-style .mat file.

    This saves the arrayobjects in the given dictionary to a matlab
    style .mat file.

    appendmat  - if true, appends '.mat' extension to filename, if not present
    format     - '4' for matlab 4 mat files, '5' for matlab 5 onwards
    """
    file_is_string = isinstance(file_name, basestring)
    if file_is_string:
        if appendmat and file_name[-4:] != ".mat":
            file_name = file_name + ".mat"
        file_stream = open(file_name, 'wb')
    else:
        try:
            file_name.write('')
        except AttributeError:
            raise IOError, 'Writer needs file name or writeable '\
                           'file-like object'
        file_stream = file_name

    if format == '4':
        MW = MatFile4Writer(file_stream)
    elif format == '5':
        MW = MatFile5Writer(file_stream, unicode_strings=True)
    else:
        raise ValueError, 'Format should be 4 or 5'
    MW.put_variables(mdict)
    if file_is_string:
        file_stream.close()
