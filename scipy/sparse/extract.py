"""Functions to extract parts of sparse matrices
"""

from __future__ import division, print_function, absolute_import

__docformat__ = "restructuredtext en"

__all__ = ['find', 'tril', 'triu', 'unique']

from itertools import izip

import numpy as np

from .coo import coo_matrix


def find(A):
    """Return the indices and values of the nonzero elements of a matrix

    Parameters
    ----------
    A : dense or sparse matrix
        Matrix whose nonzero elements are desired.

    Returns
    -------
    (I,J,V) : tuple of arrays
        I,J, and V contain the row indices, column indices, and values
        of the nonzero matrix entries.


    Examples
    --------
    >>> from scipy.sparse import csr_matrix, find
    >>> A = csr_matrix([[7.0, 8.0, 0],[0, 0, 9.0]])
    >>> find(A)
    (array([0, 0, 1], dtype=int32), array([0, 1, 2], dtype=int32), array([ 7.,  8.,  9.]))

    """

    A = coo_matrix(A, copy=True)
    A.sum_duplicates()
    # remove explicit zeros
    nz_mask = A.data != 0
    return A.row[nz_mask], A.col[nz_mask], A.data[nz_mask]


def tril(A, k=0, format=None):
    """Return the lower triangular portion of a matrix in sparse format

    Returns the elements on or below the k-th diagonal of the matrix A.
        - k = 0 corresponds to the main diagonal
        - k > 0 is above the main diagonal
        - k < 0 is below the main diagonal

    Parameters
    ----------
    A : dense or sparse matrix
        Matrix whose lower trianglar portion is desired.
    k : integer : optional
        The top-most diagonal of the lower triangle.
    format : string
        Sparse format of the result, e.g. format="csr", etc.

    Returns
    -------
    L : sparse matrix
        Lower triangular portion of A in sparse format.

    See Also
    --------
    triu : upper triangle in sparse format

    Examples
    --------
    >>> from scipy.sparse import csr_matrix, tril
    >>> A = csr_matrix([[1, 2, 0, 0, 3], [4, 5, 0, 6, 7], [0, 0, 8, 9, 0]],
    ...                dtype='int32')
    >>> A.toarray()
    array([[1, 2, 0, 0, 3],
           [4, 5, 0, 6, 7],
           [0, 0, 8, 9, 0]])
    >>> tril(A).toarray()
    array([[1, 0, 0, 0, 0],
           [4, 5, 0, 0, 0],
           [0, 0, 8, 0, 0]])
    >>> tril(A).nnz
    4
    >>> tril(A, k=1).toarray()
    array([[1, 2, 0, 0, 0],
           [4, 5, 0, 0, 0],
           [0, 0, 8, 9, 0]])
    >>> tril(A, k=-1).toarray()
    array([[0, 0, 0, 0, 0],
           [4, 0, 0, 0, 0],
           [0, 0, 0, 0, 0]])
    >>> tril(A, format='csc')
    <3x5 sparse matrix of type '<class 'numpy.int32'>'
            with 4 stored elements in Compressed Sparse Column format>

    """

    # convert to COOrdinate format where things are easy
    A = coo_matrix(A, copy=False)
    mask = A.row + k >= A.col
    return _masked_coo(A, mask).asformat(format)


def triu(A, k=0, format=None):
    """Return the upper triangular portion of a matrix in sparse format

    Returns the elements on or above the k-th diagonal of the matrix A.
        - k = 0 corresponds to the main diagonal
        - k > 0 is above the main diagonal
        - k < 0 is below the main diagonal

    Parameters
    ----------
    A : dense or sparse matrix
        Matrix whose upper trianglar portion is desired.
    k : integer : optional
        The bottom-most diagonal of the upper triangle.
    format : string
        Sparse format of the result, e.g. format="csr", etc.

    Returns
    -------
    L : sparse matrix
        Upper triangular portion of A in sparse format.

    See Also
    --------
    tril : lower triangle in sparse format

    Examples
    --------
    >>> from scipy.sparse import csr_matrix, triu
    >>> A = csr_matrix([[1, 2, 0, 0, 3], [4, 5, 0, 6, 7], [0, 0, 8, 9, 0]],
    ...                dtype='int32')
    >>> A.toarray()
    array([[1, 2, 0, 0, 3],
           [4, 5, 0, 6, 7],
           [0, 0, 8, 9, 0]])
    >>> triu(A).toarray()
    array([[1, 2, 0, 0, 3],
           [0, 5, 0, 6, 7],
           [0, 0, 8, 9, 0]])
    >>> triu(A).nnz
    8
    >>> triu(A, k=1).toarray()
    array([[0, 2, 0, 0, 3],
           [0, 0, 0, 6, 7],
           [0, 0, 0, 9, 0]])
    >>> triu(A, k=-1).toarray()
    array([[1, 2, 0, 0, 3],
           [4, 5, 0, 6, 7],
           [0, 0, 8, 9, 0]])
    >>> triu(A, format='csc')
    <3x5 sparse matrix of type '<class 'numpy.int32'>'
            with 8 stored elements in Compressed Sparse Column format>

    """

    # convert to COOrdinate format where things are easy
    A = coo_matrix(A, copy=False)
    mask = A.row + k <= A.col
    return _masked_coo(A, mask).asformat(format)


def _masked_coo(A, mask):
    row = A.row[mask]
    col = A.col[mask]
    data = A.data[mask]
    return coo_matrix((data, (row, col)), shape=A.shape, dtype=A.dtype)


def unique(mat, return_indices=False, return_inverse=False, return_counts=False):
    """
    Find the unique elements of a sparse matrix

    Returns the sorted unique elements of a sparse matrix. There are two
    optional outputs in addition to the unique elements: the indices of the
    input matrix that give the unique values, and the indices of the unique
    matrix that reconstruct the input matrix.

    Parameters
    ----------
    mat : sparse matrix
        Input matrix. This will be converted to the CSR representation
        internally.
    return_indices : bool, optional
        If True, also return the indices of `mat` that result in the unique
        array.
    return_inverse : bool, optional
        If True, also return the indices of the unique array that can be used
        to reconstruct `mat`.
    return_counts : bool, optional
        If True, also return the number of times each unique value comes up
        in `mat`.

    Returns
    -------
    unique : ndarray
        The sorted unique values.
    unique_indices : ndarray, optional
        The indices of the first occurrences of the unique values in the
        (flattened) original array. Only provided if `return_indices` is True.
    unique_inverse : ndarray, optional
        The indices to reconstruct the (flattened) original array from the
        unique array. Only provided if `return_inverse` is True.

        Note that, because the matrix is sparse, the full array of indices is
        not returned. Instead, a 2 x n array i is returned such that the
        following code will reproduce the original matrix:
            m = np.zeros(np.prod(mat.shape), dtype=mat.dtype) # or sparse matrix
            m[i[0]] = unique[i[1]]
            np.reshape(m, mat.shape)
    unique_counts : ndarray, optional
        The number of times each of the unique values comes up in the
        original array. Only provided if `return_counts` is True.

    See Also
    --------
    numpy.lib.arraysetops.unique : Basis for this function, but only works for
                                   dense matrices/arrays.
    """
    # Convert to CSR because this format has a .data matrix, which is the main
    # thing we need, and which is stored in sorted order of rows -> columns.
    # This means that np.unique returns the indices and inverse in terms of a
    # sensibly linearized matrix. (The nonzero indices are also returned in
    # row -> column order, which is useful for the return_inverse and
    # return_indices options.) Also, CSR is fairly memory-efficient and quick to
    # convert to from other formats.
    mat = mat.tocsr()
    size = mat.shape[0] * mat.shape[1]  # mat.size just gives nnz

    unique_data = np.unique(mat.data, return_indices, return_inverse,
                            return_counts)

    # If there are no zeros, we can just pretend we're operating on a normal
    # dense array. All we have to do then is adapt the inverse return value, if
    # there is one, to our special sparse inverse format.
    if mat.nnz == size:
        if return_inverse:
            inv_index = (2 if return_indices else 1)
            inverse = np.vstack((range(size), unique_data[inv_index]))
            unique_data = list(unique_data)
            unique_data[inv_index] = inverse
            unique_data = tuple(unique_data)
        return unique_data

    # OK, there are some zeros.
    # Our lives are simplest if the only thing requested was the unique values.
    if not isinstance(unique_data, tuple):
        # We got here because there are zeros, so we know 0 should be in the
        # list of unique values.
        return np.insert(unique_data, 0, 0.0)

    # If more values were requested, process other return values in the tuple
    # as necessary.

    unique_data = list(reversed(unique_data))
    unique_values = unique_data.pop()
    unique_values = np.insert(unique_values, 0, 0.0)
    ret = (unique_values,)

    # Offset returned indices to account for missing zero entries.
    if return_indices or return_inverse:
        if return_indices:
            indices = unique_data.pop()
        if return_inverse:
            inverse = unique_data.pop()

            # We're going to use inverse[0] as the array indices at which
            # values in the original matrix reside, and inverse[1] as the
            # indices in the unique array from which to draw those values.
            # We must add 1 to inverse[1] to account for the 0 in the initial
            # position.

            # The indices for the inverse matrix aren't accounting for the
            # presence of a zero value at the start of the uniques list.
            inverse_unique_indices = inverse + 1
            # Initialize positions in original matrix to values' current
            # positions in the inverse array. As we detect 0 values in the
            # original matrix, we'll increase these indices accordingly.
            inverse_orig_pos_indices = np.array(range(len(inverse)))

        first_zero = None
        offset = 0
        mat.sort_indices()
        nonzero = mat.nonzero()

        for i, (row, col) in enumerate(izip(nonzero[0], nonzero[1])):
            offset_i = i + offset
            flattened_index = row * mat.shape[1] + col
            difference = flattened_index - offset_i
            if difference > 0:  # We've found one or more zero entries!
                if return_indices:
                    indices[np.where(indices >= offset_i)] += difference
                    if first_zero is None:
                        first_zero = i
                        indices = np.insert(indices, 0, first_zero)
                if return_inverse:
                    inverse_orig_pos_indices[
                        np.where(inverse_orig_pos_indices >= offset_i)
                        ] += difference
                offset += difference

        if return_indices:
            if first_zero is None:  # Didn't find a zero amidst nonzeros;
                try:                # all zeros are trailing
                    first_zero = (nonzero[0][-1] * mat.shape[0] +
                                  nonzero[1][-1] + 1)
                except IndexError:  # nonzero is empty; mat is all zeros!
                    indices.dtype = int  # Fix dtype from np.unique call on []
                    first_zero = 0
                indices = np.insert(indices, 0, first_zero)
            ret += (indices,)

        if return_inverse:
            inverse = np.vstack((inverse_orig_pos_indices,
                                 inverse_unique_indices))
            # Fix dtype, just in case mat is all zeros and unique returned a
            # weird dtype for inverse.
            inverse.dtype = int
            ret += (inverse,)

    # Add counts for 0 value.
    if return_counts:
        counts = unique_data.pop()
        counts = np.insert(counts, 0, size - mat.nnz)
        ret += (counts,)

    return ret
