import numpy as np
import scipy as sp

__all__ = ['save_npz', 'load_npz', 'load_npz_array']


# Make loading safe vs. malicious input
PICKLE_KWARGS = dict(allow_pickle=False)


def save_npz(file, matrix, compressed=True):
    """ Save a sparse matrix or array to a file using ``.npz`` format.

    Parameters
    ----------
    file : str or file-like object
        Either the file name (string) or an open file (file-like object)
        where the data will be saved. If file is a string, the ``.npz``
        extension will be appended to the file name if it is not already
        there.
    matrix: spmatrix or sparray
        The sparse matrix or array to save.
        Supported formats: ``csc``, ``csr``, ``bsr``, ``dia`` or coo``.
    compressed : bool, optional
        Allow compressing the file. Default: True

    See Also
    --------
    scipy.sparse.load_npz: Load a sparse matrix from a file using ``.npz`` format.
    numpy.savez: Save several arrays into a ``.npz`` archive.
    numpy.savez_compressed : Save several arrays into a compressed ``.npz`` archive.

    Examples
    --------
    Store sparse matrix to disk, and load it again:

    >>> import numpy as np
    >>> import scipy as sp
    >>> sparse_matrix = sp.sparse.csc_matrix(np.array([[0, 0, 3], [4, 0, 0]]))
    >>> sparse_matrix
    <2x3 sparse matrix of type '<class 'numpy.int64'>'
       with 2 stored elements in Compressed Sparse Column format>
    >>> sparse_matrix.toarray()
    array([[0, 0, 3],
           [4, 0, 0]], dtype=int64)

    >>> sp.sparse.save_npz('/tmp/sparse_matrix.npz', sparse_matrix)
    >>> sparse_matrix = sp.sparse.load_npz('/tmp/sparse_matrix.npz')

    >>> sparse_matrix
    <2x3 sparse matrix of type '<class 'numpy.int64'>'
       with 2 stored elements in Compressed Sparse Column format>
    >>> sparse_matrix.toarray()
    array([[0, 0, 3],
           [4, 0, 0]], dtype=int64)
    """
    arrays_dict = {}
    if matrix.format in ('csc', 'csr', 'bsr'):
        arrays_dict.update(indices=matrix.indices, indptr=matrix.indptr)
    elif matrix.format == 'dia':
        arrays_dict.update(offsets=matrix.offsets)
    elif matrix.format == 'coo':
        arrays_dict.update(row=matrix.row, col=matrix.col)
    else:
        raise NotImplementedError(f'Save is not implemented for sparse matrix of format {matrix.format}.')
    arrays_dict.update(
        format=matrix.format.encode('ascii'),
        shape=matrix.shape,
        data=matrix.data
    )
    if compressed:
        np.savez_compressed(file, **arrays_dict)
    else:
        np.savez(file, **arrays_dict)


def load_npz_array(file):
    """ Load a sparse array from a file using ``.npz`` format.

    Parameters
    ----------
    file : str or file-like object
        Either the file name (string) or an open file (file-like object)
        where the data will be loaded.

    Returns
    -------
    result : csc_array, csr_array, bsr_array, dia_array or coo_array
        A sparse array containing the loaded data.

    Raises
    ------
    OSError
        If the input file does not exist or cannot be read.

    See Also
    --------
    scipy.sparse.save_npz: Save a sparse array to a file using ``.npz`` format.
    numpy.load: Load several arrays from a ``.npz`` archive.

    Examples
    --------
    Store sparse array to disk, and load it again:

    >>> import numpy as np
    >>> import scipy as sp
    >>> sparse_array = sp.sparse.csc_array(np.array([[0, 0, 3], [4, 0, 0]]))
    >>> sparse_array
    <2x3 sparse array of type '<class 'numpy.int64'>'
       with 2 stored elements in Compressed Sparse Column format>
    >>> sparse_array.toarray()
    array([[0, 0, 3],
           [4, 0, 0]], dtype=int64)

    >>> sp.sparse.save_npz('/tmp/sparse_array.npz', sparse_array)
    >>> sparse_array = sp.sparse.load_npz_array('/tmp/sparse_array.npz')

    >>> sparse_array
    <2x3 sparse array of type '<class 'numpy.int64'>'
        with 2 stored elements in Compressed Sparse Column format>
    >>> sparse_array.toarray()
    array([[0, 0, 3],
           [4, 0, 0]], dtype=int64)
    """
    return _load_npz(file)


def _load_npz(file, as_sparray=True):
    with np.load(file, **PICKLE_KWARGS) as loaded:
        try:
            sparse_format = loaded['format']
        except KeyError as e:
            raise ValueError(f'The file {file} does not contain a sparse matrix.') from e

        sparse_format = sparse_format.item()
        arr_mat = 'array' if as_sparray else 'matrix'

        if not isinstance(sparse_format, str):
            # Play safe with Python 2 vs 3 backward compatibility;
            # files saved with SciPy < 1.0.0 may contain unicode or bytes.
            sparse_format = sparse_format.decode('ascii')

        try:
            cls = getattr(sp.sparse, f'{sparse_format}_{arr_mat}')
        except AttributeError as e:
            raise ValueError(f'Unknown matrix format "{sparse_format}"') from e

        if sparse_format in ('csc', 'csr', 'bsr'):
            return cls((loaded['data'], loaded['indices'], loaded['indptr']), shape=loaded['shape'])
        elif sparse_format == 'dia':
            return cls((loaded['data'], loaded['offsets']), shape=loaded['shape'])
        elif sparse_format == 'coo':
            return cls((loaded['data'], (loaded['row'], loaded['col'])), shape=loaded['shape'])
        else:
            raise NotImplementedError('Load is not implemented for '
                                      'sparse matrix of format {}.'.format(sparse_format))


def load_npz(file):
    """ Load a sparse matrix from a file using ``.npz`` format.

    Parameters
    ----------
    file : str or file-like object
        Either the file name (string) or an open file (file-like object)
        where the data will be loaded.

    Returns
    -------
    result : csc_matrix, csr_matrix, bsr_matrix, dia_matrix or coo_matrix
        A sparse matrix containing the loaded data.

    Raises
    ------
    OSError
        If the input file does not exist or cannot be read.

    See Also
    --------
    scipy.sparse.save_npz: Save a sparse matrix to a file using ``.npz`` format.
    numpy.load: Load several arrays from a ``.npz`` archive.

    Examples
    --------
    Store sparse matrix to disk, and load it again:

    >>> import numpy as np
    >>> import scipy.sparse
    >>> sparse_matrix = scipy.sparse.csc_matrix(np.array([[0, 0, 3], [4, 0, 0]]))
    >>> sparse_matrix
    <2x3 sparse matrix of type '<class 'numpy.int64'>'
       with 2 stored elements in Compressed Sparse Column format>
    >>> sparse_matrix.toarray()
    array([[0, 0, 3],
           [4, 0, 0]], dtype=int64)

    >>> scipy.sparse.save_npz('/tmp/sparse_matrix.npz', sparse_matrix)
    >>> sparse_matrix = scipy.sparse.load_npz('/tmp/sparse_matrix.npz')

    >>> sparse_matrix
    <2x3 sparse matrix of type '<class 'numpy.int64'>'
        with 2 stored elements in Compressed Sparse Column format>
    >>> sparse_matrix.toarray()
    array([[0, 0, 3],
           [4, 0, 0]], dtype=int64)
    """
    return _load_npz(file, as_sparray=False)
