"""
Module to read / write Fortran unformatted files.

This module allows one to read sequential unformatted Fortran files, provided the structure of the file is known. Sequential files typically have a record size written, the contents of the record, and the size a second time. Files from compilers without both sizes are not supported yet.

This is in the spirit of code written by Neil Martinsen-Burrell and Joe Zuntz.

As Joe Zuntz states, "Like all Fortran-later language interoperability code documentation should, this docstring ends with an exhortation to please stop using fortran."

"""

__all__ = ['fortran_file']

import numpy as np
import warnings

class fortran_file(object):
    """
    A file object for unformatted Fortran files.

    Parameters
    ----------
    filename: file or str
        Open file object or filename.
    mode : {'r', 'w'}, optional
        Read-write mode, default is 'r'.
    header_dtype : data-type
        Data type of the header. Size and endiness must match the input/output file.

    Examples
    --------
    To create an unformatted Fortran file:

        >>> from scipy.io import fortran
        >>> f = fortran.fortran_file('test.unf', 'w')
        >>> f.writeRecord(np.array([1,2,3,4,5],dtype=np.int32))
        >>> f.writeRecord(np.linspace(0,1,20).reshape((5,-1)))
        >>> f.close()

    To read this file:

        >>> from scipy.io import fortran
        >>> f = fortran.fortran_file('test.unf', 'r')
        >>> print(f.readInts(dtype=np.int32))
        [1 2 3 4 5]
        >>> print(f.readReals(dtype=np.float).reshape((5,-1)))
        [[ 0.          0.26315789  0.52631579  0.78947368]
         [ 0.05263158  0.31578947  0.57894737  0.84210526]
         [ 0.10526316  0.36842105  0.63157895  0.89473684]
         [ 0.15789474  0.42105263  0.68421053  0.94736842]
         [ 0.21052632  0.47368421  0.73684211  1.        ]]
        >>> f.close()

    """
    def __init__(self, filename, mode='r', header_dtype=np.int64):
        header_dtype = np.dtype(header_dtype)
        if mode not in 'rw' or len(mode) != 1:
            raise ValueError('mode must be either r or w')

        if hasattr(filename, 'seek'):
            self._fp = filename
        else:
            self._fp = open(filename, '%sb' % mode)

        self._header_dtype=header_dtype

    def _readSize(self):
        return np.fromfile(self._fp, dtype=self._header_dtype, count=1)

    def writeRecord(self, s):
        """
        Write a record (including header/footer) to the file.

        Parameters
        ----------
        s : array_like
           The data to write.
        """
        s = np.array(s, order='F')
        np.array([s.nbytes],dtype=self._header_dtype).tofile(self._fp)
        s.tofile(self._fp)
        np.array([s.nbytes],dtype=self._header_dtype).tofile(self._fp)

    def readRecord(self, dtype=None):
        """
        Reads a record of a given type from the file.

        Parameters
        ----------
        dtype : data-type
            Data type specifying the size and endiness of the data.

        Returns
        -------
        data : ndarray
            A one-diminsional array object.

        Notes
        -----
        If the record contains a multi-dimensional array, calling reshape or resize will restructure the array to the correct size.
        Since Fortran multidiminsional arrays are stored in column-major format, this may have some non-intuitive consequences. If the variable was declared as 'INTEGER var(5,4)', for example, var could be read with 'readRecord(dtype=np.integer).reshape( (5,4) )' since Python uses row-major ordering of indices. For two-diminsional arrays, one can transpose to obtain the indices in the same order as in Fortran.

        See Also
        --------
        readReals
        readInts
        """
        if dtype is None:
            raise ValueError('Must specify dtype')
        dtype = np.dtype(dtype)

        firstSize = self._readSize()
        if firstSize % dtype.itemsize != 0:
            raise ValueError('Size obtained is not a multiple of the dtype given.')
        data = np.fromfile(self._fp, dtype=dtype, count=int(firstSize / dtype.itemsize))
        secondSize = self._readSize()
        if firstSize != secondSize:
            raise IOError('Sizes do not agree in the header and footer for this record - check header dtype')
        return data

    def readInts(self, dtype='i4'):
        """
        Reads a record of a given type from the file, defaulting to an integer type (INTEGER*4 in Fortran)

        Parameters
        ----------
        dtype : data-type
            Data type specifying the size and endiness of the data.

        Returns
        -------
        data : ndarray
            A one-diminsional array object.

        See Also
        --------
        readReals
        readRecord
        """
        return self.readRecord(dtype=dtype)

    def readReals(self, dtype='f8'):
        """
        Reads a record of a given type from the file, defaulting to a floating point number (real*8 in Fortran)

        Parameters
        ----------
        dtype : data-type
            Data type specifying the size and endiness of the data.

        Returns
        -------
        data : ndarray
            A one-diminsional array object.

        See Also
        --------
        readInts
        readRecord
        """
        return self.readRecord(dtype=dtype)

    def close(self):
        """
        Closes the file. It is unsupported to call any other methods off this object after closing it. Note that this class supports the 'with' statement in modern versions of Python, to call this automatically
        """
        self._fp.close()

    def __enter__(self):
        return self

    def __exit__(self, type, value, tb):
        self.close()
