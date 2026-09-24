import os
import numpy as np
import tempfile

import pytest
from pytest import raises as assert_raises
from numpy.testing import assert_equal, assert_

from scipy.sparse import (save_npz, load_npz, csc_array, csr_array, bsr_array,
                          dia_array, coo_array, dok_array, lil_array)


DATA_DIR = os.path.join(os.path.dirname(__file__), 'data')


def _save_and_load(matrix):
    fd, tmpfile = tempfile.mkstemp(suffix='.npz')
    os.close(fd)
    try:
        save_npz(tmpfile, matrix)
        loaded_matrix = load_npz(tmpfile)
    finally:
        os.remove(tmpfile)
    return loaded_matrix

def _check_save_and_load(dense_matrix):
    for matrix_class in [csc_array, csr_array, bsr_array, dia_array, coo_array]:
        matrix = matrix_class(dense_matrix)
        loaded_matrix = _save_and_load(matrix)
        assert_(type(loaded_matrix) is matrix_class)
        assert_(loaded_matrix.shape == dense_matrix.shape)
        assert_(loaded_matrix.dtype == dense_matrix.dtype)
        assert_equal(loaded_matrix.toarray(), dense_matrix)

def test_save_and_load_random():
    N = 10
    np.random.seed(0)
    dense_matrix = np.random.random((N, N))
    dense_matrix[dense_matrix > 0.7] = 0
    _check_save_and_load(dense_matrix)

def test_save_and_load_empty():
    dense_matrix = np.zeros((4,6))
    _check_save_and_load(dense_matrix)

def test_save_and_load_one_entry():
    dense_matrix = np.zeros((4,6))
    dense_matrix[1,2] = 1
    _check_save_and_load(dense_matrix)

@pytest.mark.parametrize("value", [0, 1.2])
@pytest.mark.parametrize("ndim", [1, 2, 3])
def test_nd_coo_format(ndim, value):
    A = coo_array([value]).reshape((1,) * ndim)

    #save/load array
    fd, tmpfile = tempfile.mkstemp(suffix='.npz')
    os.close(fd)
    try:
        save_npz(tmpfile, A)
        loaded_A = load_npz(tmpfile)
    finally:
        os.remove(tmpfile)

    assert isinstance(loaded_A, coo_array)
    assert_(loaded_A.shape == A.shape)
    assert_equal(A.toarray(), loaded_A.toarray())

def test_malicious_load():
    class Executor:
        def __reduce__(self):
            return (assert_, (False, 'unexpected code execution'))

    fd, tmpfile = tempfile.mkstemp(suffix='.npz')
    os.close(fd)
    try:
        np.savez(tmpfile, format=Executor())

        # Should raise a ValueError, not execute code
        assert_raises(ValueError, load_npz, tmpfile)
    finally:
        os.remove(tmpfile)

@pytest.mark.parametrize("container", [dok_array, lil_array])
def test_implemented_error(container):
    # Attempts to save an unsupported type and checks that an
    # NotImplementedError is raised.

    x = container((2,3))
    x[0,1] = 1

    with pytest.raises(NotImplementedError, match="convert.*before saving"):
        save_npz("x.npz", x)
