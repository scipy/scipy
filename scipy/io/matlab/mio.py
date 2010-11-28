"""Module for reading and writing MATLAB .mat files"""
# Authors: Travis Oliphant, Matthew Brett

"""
Module for reading and writing matlab (TM) .mat files
"""

import os
import sys
import warnings

from numpy.compat import asbytes

from miobase import get_matfile_version, docfiller
from mio4 import MatFile4Reader, MatFile4Writer
from mio5 import MatFile5Reader, MatFile5Writer

__all__ = ['find_mat_file', 'mat_reader_factory', 'loadmat', 'savemat']

@docfiller
def find_mat_file(file_name, appendmat=True):
    ''' Try to find .mat file on system path

    Parameters
    ----------
    file_name : str
       file name for mat file
    %(append_arg)s

    Returns
    -------
    full_name : string
       possibly modified name after path search
    '''
    warnings.warn('Searching for mat files on python system path will be ' +
                  'removed in next version of scipy',
                   DeprecationWarning, stacklevel=2)
    if appendmat and file_name.endswith(".mat"):
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


def _open_file(file_like, appendmat):
    ''' Open `file_like` and return as file-like object '''
    if isinstance(file_like, basestring):
        try:
            return open(file_like, 'rb')
        except IOError:
            pass
        if appendmat and not file_like.endswith('.mat'):
            try:
                return open(file_like + '.mat', 'rb')
            except IOError:
                pass
        # search the python path - we'll remove this soon
        full_name = find_mat_file(file_like, appendmat)
        if full_name is None:
            raise IOError("%s not found on the path."
                          % file_like)
        return open(full_name, 'rb')
    # not a string - maybe file-like object
    try:
        file_like.read(0)
    except AttributeError:
        raise IOError('Reader needs file name or open file-like object')
    return file_like


@docfiller
def mat_reader_factory(file_name, appendmat=True, **kwargs):
    """Create reader for matlab .mat format files

    Parameters
    ----------
    %(file_arg)s
    %(append_arg)s
    %(load_args)s
    %(struct_arg)s

    Returns
    -------
    matreader : MatFileReader object
       Initialized instance of MatFileReader class matching the mat file
       type detected in `filename`.
    """
    byte_stream = _open_file(file_name, appendmat)
    mjv, mnv = get_matfile_version(byte_stream)
    if mjv == 0:
        return MatFile4Reader(byte_stream, **kwargs)
    elif mjv == 1:
        return MatFile5Reader(byte_stream, **kwargs)
    elif mjv == 2:
        raise NotImplementedError('Please use HDF reader for matlab v7.3 files')
    else:
        raise TypeError('Did not recognize version %s' % mjv)

@docfiller
def loadmat(file_name,  mdict=None, appendmat=True, **kwargs):
    """
    Load MATLAB file

    Parameters
    ----------
    %(file_arg)s
    m_dict : dict, optional
        Dictionary in which to insert matfile variables.
    %(append_arg)s
    %(load_args)s
    %(struct_arg)s
    variable_names : None or sequence
        If None (the default) - read all variables in file. Otherwise
        `variable_names` should be a sequence of strings, giving names of the
        matlab variables to read from the file.  The reader will skip any
        variable with a name not in this sequence, possibly saving some read
        processing.

    Returns
    -------
    mat_dict : dict
       dictionary with variable names as keys, and loaded matrices as
       values

    Notes
    -----
    v4 (Level 1.0), v6 and v7 to 7.2 matfiles are supported.

    You will need an HDF5 python library to read matlab 7.3 format mat
    files.  Because scipy does not supply one, we do not implement the
    HDF5 / 7.3 interface here.

    """
    variable_names = kwargs.pop('variable_names', None)
    MR = mat_reader_factory(file_name, appendmat, **kwargs)
    matfile_dict = MR.get_variables(variable_names)
    if mdict is not None:
        mdict.update(matfile_dict)
    else:
        mdict = matfile_dict
    if isinstance(file_name, basestring):
        MR.mat_stream.close()
    return mdict

@docfiller
def savemat(file_name, mdict,
            appendmat=True,
            format='5',
            long_field_names=False,
            do_compression=False,
            oned_as=None):
    """
    Save a dictionary of names and arrays into a MATLAB-style .mat file.

    This saves the array objects in the given dictionary to a MATLAB-
    style .mat file.

    Parameters
    ----------
    file_name : str or file-like object
        Name of the .mat file (.mat extension not needed if ``appendmat ==
        True``).
        Can also pass open file_like object.
    m_dict : dict
        Dictionary from which to save matfile variables.
    %(append_arg)s
    format : {'5', '4'}, string, optional
        '5' (the default) for MATLAB 5 and up (to 7.2),
        '4' for MATLAB 4 .mat files
    %(long_fields)s
    %(do_compression)s
    %(oned_as)s

    See also
    --------
    mio4.MatFile4Writer
    mio5.MatFile5Writer

    Notes
    -----
    If ``format == '4'``, `mio4.MatFile4Writer` is called, which sets
    `oned_as` to 'row' if it had been None.  If ``format == '5'``,
    `mio5.MatFile5Writer` is called, which sets `oned_as` to 'column' if
    it had been None, but first it executes:

    ``warnings.warn("Using oned_as default value ('column')" +``
                  ``" This will change to 'row' in future versions",``
                  ``FutureWarning, stacklevel=2)``

    without being more specific as to precisely when the change will take
    place.

    """
    file_is_string = isinstance(file_name, basestring)
    if file_is_string:
        if appendmat and file_name[-4:] != ".mat":
            file_name = file_name + ".mat"
        file_stream = open(file_name, 'wb')
    else:
        try:
            file_name.write(asbytes(''))
        except AttributeError:
            raise IOError('Writer needs file name or writeable '
                           'file-like object')
        file_stream = file_name

    if format == '4':
        if long_field_names:
            raise ValueError("Long field names are not available for version 4 files")
        MW = MatFile4Writer(file_stream, oned_as)
    elif format == '5':
        MW = MatFile5Writer(file_stream,
                            do_compression=do_compression,
                            unicode_strings=True,
                            long_field_names=long_field_names,
                            oned_as=oned_as)
    else:
        raise ValueError("Format should be '4' or '5'")
    MW.put_variables(mdict)
    if file_is_string:
        file_stream.close()
