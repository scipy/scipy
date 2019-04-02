import numpy
from pytest import raises as assert_raises

import scipy.ndimage as sndi


def test_binary_erosion_noninteger_iterations():
    # regression test for gh-9905, gh-9909: ValueError for
    # non integer iterations
    data = numpy.ones([1])
    assert_raises(TypeError, sndi.binary_erosion, data, iterations=0.5)


def test_binary_dilation_noninteger_iterations():
    # regression test for gh-9905, gh-9909: ValueError for
    # non integer iterations
    data = numpy.ones([1])
    assert_raises(TypeError, sndi.binary_dilation, data, iterations=0.5)


def test_binary_opening_noninteger_iterations():
    # regression test for gh-9905, gh-9909: ValueError for
    # non integer iterations
    data = numpy.ones([1])
    assert_raises(TypeError, sndi.binary_opening, data, iterations=0.5)


def test_binary_closing_noninteger_iterations():
    # regression test for gh-9905, gh-9909: ValueError for
    # non integer iterations
    data = numpy.ones([1])
    assert_raises(TypeError, sndi.binary_closing, data, iterations=0.5)
