''' Some tests for the documenting decorator '''

import numpy as np

from numpy.testing import assert_equal, assert_raises

from nose.tools import assert_true

from scipy.ndimage.doccer import docformat, filldoc

def test_docformat():
    docstring = \
    """Docstring
    %(strtest)s
    """
    param_doc = \
"""Another test
   with some indent"""
    formatted = docformat(docstring, {'strtest':param_doc})
    expected = \
"""Docstring
    Another test
       with some indent
    """
    assert_equal(formatted, expected)
