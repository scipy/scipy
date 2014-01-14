from __future__ import division, print_function, absolute_import

from numpy.testing import assert_raises

from scipy.sparse import sparsetools

def test_exception():
    assert_raises(MemoryError, sparsetools.csr.test_throw_error)
