# Copyright (C) 2022-2023 Adam Lugowski. All rights reserved.
# Use of this source code is governed by the BSD 2-clause license found in the LICENSE.txt file.
# SPDX-License-Identifier: BSD-2-Clause
"""
Matrix Market I/O with a C++ backend.
See http://math.nist.gov/MatrixMarket/formats.html
for information about the Matrix Market format.

.. versionadded:: 1.12.0
"""
import io
import os

import numpy as np
import scipy.sparse

from . import _core  # type: ignore

__version__ = _core.__version__
header = _core.header

PARALLELISM = 0
"""
Default value for the parallelism argument to mmread() and mmwrite().
0 means number of CPUs in the system.
"""

ALWAYS_FIND_SYMMETRY = False
"""
Whether mmwrite() with symmetry='AUTO' will always search for symmetry inside the matrix.
This is scipy.io._mmio.mmwrite()'s default behavior, but has a significant performance cost on large matrices.
"""

_field_to_dtype = {
    "integer": "int64",
    "unsigned-integer": "uint64",
    "real": "float64",
    "complex": "complex",
    "pattern": "float64",

    "long-integer": "int64",
    "long-unsigned-integer": "uint64",
    "long-real": "longdouble",
    "long-complex": "longcomplex",
    "long-pattern": "longdouble",
}


class _TextToBytesWrapper(io.BufferedReader):
    """
    Convert a TextIOBase string stream to a byte stream.
    """

    def __init__(self, text_io_buffer, encoding=None, errors=None, **kwargs):
        super(_TextToBytesWrapper, self).__init__(text_io_buffer, **kwargs)
        self.encoding = encoding or text_io_buffer.encoding or 'utf-8'
        self.errors = errors or text_io_buffer.errors or 'strict'

    def __del__(self):
        # do not close the wrapped stream
        self.detach()

    def _encoding_call(self, method_name, *args, **kwargs):
        raw_method = getattr(self.raw, method_name)
        val = raw_method(*args, **kwargs)
        return val.encode(self.encoding, errors=self.errors)

    def read(self, size=-1):
        return self._encoding_call('read', size)

    def read1(self, size=-1):
        return self._encoding_call('read1', size)

    def peek(self, size=-1):
        return self._encoding_call('peek', size)


def _read_body_array(cursor, long_type):
    """
    Read MatrixMarket array body
    """
    vals = np.zeros(cursor.header.shape, dtype=_field_to_dtype.get(("long-" if long_type else "")+cursor.header.field))
    _core.read_body_array(cursor, vals)
    return vals


def _read_body_coo(cursor, long_type, generalize_symmetry=True):
    """
    Read MatrixMarket coordinate body
    """
    index_dtype = "int32"
    if cursor.header.nrows >= 2**31 or cursor.header.ncols >= 2**31:
        # Dimensions are too large to fit in int32
        index_dtype = "int64"

    i = np.zeros(cursor.header.nnz, dtype=index_dtype)
    j = np.zeros(cursor.header.nnz, dtype=index_dtype)
    data = np.zeros(cursor.header.nnz, dtype=_field_to_dtype.get(("long-" if long_type else "")+cursor.header.field))

    _core.read_body_coo(cursor, i, j, data)

    if generalize_symmetry and cursor.header.symmetry != "general":
        off_diagonal_mask = (i != j)
        off_diagonal_rows = i[off_diagonal_mask]
        off_diagonal_cols = j[off_diagonal_mask]
        off_diagonal_data = data[off_diagonal_mask]

        if cursor.header.symmetry == "skew-symmetric":
            off_diagonal_data *= -1
        elif cursor.header.symmetry == "hermitian":
            off_diagonal_data = off_diagonal_data.conjugate()

        i = np.concatenate((i, off_diagonal_cols))
        j = np.concatenate((j, off_diagonal_rows))
        data = np.concatenate((data, off_diagonal_data))

    return (data, (i, j)), cursor.header.shape


def _get_read_cursor(source, parallelism=None):
    """
    Open file for reading.
    """
    if parallelism is None:
        parallelism = PARALLELISM

    try:
        source = os.fspath(source)
        # It's a file path
        is_path = True
    except TypeError:
        is_path = False

    if is_path:
        path = str(source)
        if path.endswith('.gz'):
            import gzip
            source = gzip.GzipFile(path, 'r')
        elif path.endswith('.bz2'):
            import bz2
            source = bz2.BZ2File(path, 'rb')
        else:
            return _core.open_read_file(path, parallelism)

    # Stream object.
    if hasattr(source, "read"):
        if isinstance(source, io.TextIOBase):
            source = _TextToBytesWrapper(source)
        return _core.open_read_stream(source, parallelism)
    else:
        raise TypeError("Unknown source type")


def _get_write_cursor(target, h=None, comment=None, parallelism=None, symmetry="general", precision=None):
    """
    Open file for writing.
    """
    if parallelism is None:
        parallelism = PARALLELISM
    if comment is None:
        comment = ''
    if symmetry is None:
        symmetry = "general"
    if precision is None:
        precision = -1

    if not h:
        h = header(comment=comment, symmetry=symmetry)

    try:
        target = os.fspath(target)
        # It's a file path
        return _core.open_write_file(str(target), h, parallelism, precision)
    except TypeError:
        pass

    if hasattr(target, "write"):
        # Stream object.
        if isinstance(target, io.TextIOBase):
            raise TypeError("target stream must be open in binary mode.")
        return _core.open_write_stream(target, h, parallelism, precision)
    else:
        raise TypeError("Unknown source object")


def _apply_field(data, field, no_pattern=False):
    """
    Ensure that a NumPy array has a dtype compatible with the specified MatrixMarket field type.

    Parameters
    ----------
    data : ndarray
        Input array.

    field : str
        Matrix Market field, such as 'real', 'complex', 'integer', 'pattern'.

    no_pattern : bool, optional
        Whether an empty array may be returned for a 'pattern' field.

    Returns
    -------
    data : ndarray
        Input data if no conversion necessary, or a converted version
    """

    if field is None:
        return data
    if field == "pattern":
        if no_pattern:
            return data
        else:
            return np.zeros(0)

    dtype = _field_to_dtype.get(field, None)
    if dtype is None:
        raise ValueError("Invalid field.")

    return np.asarray(data, dtype=dtype)


def _validate_symmetry(symmetry):
    """
    Check that the symmetry parameter is one that MatrixMarket allows..
    """
    if symmetry is None:
        return "general"

    symmetry = str(symmetry).lower()
    symmetries = ["general", "symmetric", "skew-symmetric", "hermitian"]
    if symmetry not in symmetries:
        raise ValueError("Invalid symmetry. Must be one of: " + ", ".join(symmetries))

    return symmetry


def mmread(source, parallelism=None, long_type=False):
    """
    Reads the contents of a Matrix Market file-like 'source' into a matrix (C++ implementation).

    Parameters
    ----------
    source : str or file-like
        Matrix Market filename (extensions .mtx, .mtz.gz)
        or open file-like object.
    parallelism : int, optional
        Number of threads to use. Default is the number of CPUs in the system.

        .. versionadded:: 1.12.0

    long_type : bool, optional
        Whether to use 'longdouble' and 'longcomplex' extended-precision floating-point number types.
        Default is False.

        .. versionadded:: 1.12.0

    Returns
    -------
    a : ndarray or coo_matrix
        Dense or sparse matrix depending on the matrix format in the
        Matrix Market file.

    Notes
    -----
    .. versionadded:: 1.12.0

    Examples
    --------
    >>> from io import StringIO
    >>> from scipy.io._fast_matrix_market import mmread

    >>> text = '''%%MatrixMarket matrix coordinate real general
    ...  5 5 7
    ...  2 3 1.0
    ...  3 4 2.0
    ...  3 5 3.0
    ...  4 1 4.0
    ...  4 2 5.0
    ...  4 3 6.0
    ...  4 4 7.0
    ... '''

    ``mmread(source)`` returns the data as sparse matrix in COO format.

    >>> m = mmread(StringIO(text))
    >>> m
    <5x5 sparse matrix of type '<class 'numpy.float64'>'
    with 7 stored elements in COOrdinate format>
    >>> m.A
    array([[0., 0., 0., 0., 0.],
           [0., 0., 1., 0., 0.],
           [0., 0., 0., 2., 3.],
           [4., 5., 6., 7., 0.],
           [0., 0., 0., 0., 0.]])
    """
    cursor = _get_read_cursor(source, parallelism)

    if cursor.header.format == "array":
        return _read_body_array(cursor, long_type=long_type)
    else:
        from scipy.sparse import coo_matrix
        triplet, shape = _read_body_coo(cursor, long_type=long_type, generalize_symmetry=True)
        return coo_matrix(triplet, shape=shape)


def mmwrite(target, a, comment=None, field=None, precision=None, symmetry="AUTO", parallelism=None):
    r"""
    Writes the sparse or dense array `a` to Matrix Market file-like `target` (C++ implementation).

    Parameters
    ----------
    target : str or file-like
        Matrix Market filename (extension .mtx) or open file-like object.
    a : array like
        Sparse or dense 2-D array.
    comment : str, optional
        Comments to be prepended to the Matrix Market file.
    field : None or str, optional
        Either 'real', 'complex', 'pattern', or 'integer'.
    precision : None or int, optional
        Number of digits to display for real or complex values.
    symmetry : None or str, optional
        Either 'AUTO', 'general', 'symmetric', 'skew-symmetric', or 'hermitian'.
        If symmetry is None the symmetry type of 'a' is determined by its
        values. If symmetry is 'AUTO' the symmetry type of 'a' is either
        determined or set to 'general', at mmwrite's discretion.
    parallelism : int, optional
        Number of threads to use. Default is the number of CPUs in the system.

        .. versionadded:: 1.12.0

    Returns
    -------
    None

    Notes
    -----
    .. versionadded:: 1.12.0

    Examples
    --------
    >>> from io import BytesIO
    >>> import numpy as np
    >>> from scipy.sparse import coo_matrix
    >>> from scipy.io._fast_matrix_market import mmwrite

    Write a small NumPy array to a matrix market file.  The file will be
    written in the ``'array'`` format.

    >>> a = np.array([[1.0, 0, 0, 0], [0, 2.5, 0, 6.25]])
    >>> target = BytesIO()
    >>> mmwrite(target, a)
    >>> print(target.getvalue().decode('latin1'))
    %%MatrixMarket matrix array real general
    %
    2 4
    1.0000000000000000e+00
    0.0000000000000000e+00
    0.0000000000000000e+00
    2.5000000000000000e+00
    0.0000000000000000e+00
    0.0000000000000000e+00
    0.0000000000000000e+00
    6.2500000000000000e+00

    Add a comment to the output file, and set the precision to 3.

    >>> target = BytesIO()
    >>> mmwrite(target, a, comment='\n Some test data.\n', precision=3)
    >>> print(target.getvalue().decode('latin1'))
    %%MatrixMarket matrix array real general
    %
    % Some test data.
    %
    2 4
    1.000e+00
    0.000e+00
    0.000e+00
    2.500e+00
    0.000e+00
    0.000e+00
    0.000e+00
    6.250e+00

    Convert to a sparse matrix before calling ``mmwrite``.  This will
    result in the output format being ``'coordinate'`` rather than
    ``'array'``.

    >>> target = BytesIO()
    >>> mmwrite(target, coo_matrix(a), precision=3)
    >>> print(target.getvalue().decode('latin1'))
    %%MatrixMarket matrix coordinate real general
    %
    2 4 3
    1 1 1.00e+00
    2 2 2.50e+00
    2 4 6.25e+00

    Write a complex Hermitian array to a matrix market file.  Note that
    only six values are actually written to the file; the other values
    are implied by the symmetry.

    >>> z = np.array([[3, 1+2j, 4-3j], [1-2j, 1, -5j], [4+3j, 5j, 2.5]])
    >>> z
    array([[ 3. +0.j,  1. +2.j,  4. -3.j],
           [ 1. -2.j,  1. +0.j, -0. -5.j],
           [ 4. +3.j,  0. +5.j,  2.5+0.j]])

    >>> target = BytesIO()
    >>> mmwrite(target, z, precision=2)
    >>> print(target.getvalue().decode('latin1'))
    %%MatrixMarket matrix array complex hermitian
    %
    3 3
    3.00e+00 0.00e+00
    1.00e+00 -2.00e+00
    4.00e+00 3.00e+00
    1.00e+00 0.00e+00
    0.00e+00 5.00e+00
    2.50e+00 0.00e+00

    """

    if isinstance(a, list) or isinstance(a, tuple) or hasattr(a, "__array__"):
        a = np.asarray(a)

    if symmetry == "AUTO":
        if ALWAYS_FIND_SYMMETRY or (hasattr(a, "shape") and max(a.shape) < 100):
            symmetry = None
        else:
            symmetry = "general"

    if symmetry is None:
        symmetry = scipy.io._mmio.MMFile()._get_symmetry(a)

    symmetry = _validate_symmetry(symmetry)
    cursor = _get_write_cursor(target, comment=comment, parallelism=parallelism, precision=precision, symmetry=symmetry)

    if isinstance(a, np.ndarray):
        # Write dense numpy arrays
        a = _apply_field(a, field, no_pattern=True)
        _core.write_array(cursor, a)
        return

    if scipy.sparse.isspmatrix(a):
        # Write sparse scipy matrices
        if symmetry is not None and symmetry != "general":
            # A symmetric matrix only specifies the elements below the diagonal.
            # Ensure that the matrix satisfies this requirement.
            from scipy.sparse import coo_matrix
            a = a.tocoo()
            lower_triangle_mask = a.row >= a.col
            a = coo_matrix((a.data[lower_triangle_mask],
                            (a.row[lower_triangle_mask],
                             a.col[lower_triangle_mask])), shape=a.shape)

        # CSC and CSR have specialized writers.
        is_compressed = (isinstance(a, scipy.sparse.csc_matrix) or isinstance(a, scipy.sparse.csr_matrix))

        if not is_compressed:
            # convert everything except CSC/CSR to coo
            a = a.tocoo()

        data = _apply_field(a.data, field)

        if is_compressed:
            # CSC and CSR can be written directly
            is_csr = isinstance(a, scipy.sparse.csr_matrix)
            _core.write_csc(cursor, a.shape, a.indptr, a.indices, data, is_csr)
        else:
            _core.write_coo(cursor, a.shape, a.row, a.col, data)
        return

    raise ValueError("unknown matrix type: %s" % type(a))


def mminfo(source):
    """
    Return size and storage parameters from Matrix Market file-like 'source'.

    A C++ implementation of :ref:`scipy.io.mminfo <scipy.io.mminfo>`.

    Parameters
    ----------
    source : str or file-like
        Matrix Market filename (extension .mtx) or open file-like object

    Returns
    -------
    rows : int
        Number of matrix rows.
    cols : int
        Number of matrix columns.
    entries : int
        Number of non-zero entries of a sparse matrix
        or rows*cols for a dense matrix.
    format : str
        Either 'coordinate' or 'array'.
    field : str
        Either 'real', 'complex', 'pattern', or 'integer'.
    symmetry : str
        Either 'general', 'symmetric', 'skew-symmetric', or 'hermitian'.

    Examples
    --------
    >>> from io import StringIO
    >>> from scipy.io._fast_matrix_market import mminfo

    >>> text = '''%%MatrixMarket matrix coordinate real general
    ...  5 5 7
    ...  2 3 1.0
    ...  3 4 2.0
    ...  3 5 3.0
    ...  4 1 4.0
    ...  4 2 5.0
    ...  4 3 6.0
    ...  4 4 7.0
    ... '''


    ``mminfo(source)`` returns the number of rows, number of columns,
    format, field type and symmetry attribute of the source file.

    >>> mminfo(StringIO(text))
    (5, 5, 7, 'coordinate', 'real', 'general')
    """
    h = _get_read_cursor(source, 1).header
    return h.nrows, h.ncols, h.nnz, h.format, h.field, h.symmetry
