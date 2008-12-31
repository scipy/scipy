''' Some tests for the documenting decorator and support functions '''

import numpy as np

from numpy.testing import assert_equal, assert_raises

from nose.tools import assert_true

import scipy.ndimage.doccer as sndd

docstring = \
"""Docstring
    %(strtest1)s
        %(strtest2)s
     %(strtest3)s
"""
param_doc1 = \
"""Another test
   with some indent"""

param_doc2 = \
"""Another test, one line"""

param_doc3 = \
"""    Another test
       with some indent"""

doc_dict = {'strtest1':param_doc1,
            'strtest2':param_doc2,
            'strtest3':param_doc3}

filled_docstring = \
"""Docstring
    Another test
       with some indent
        Another test, one line
     Another test
       with some indent
"""


def test_unindent():
    yield assert_equal, sndd.unindent_string(param_doc1), param_doc1
    yield assert_equal, sndd.unindent_string(param_doc2), param_doc2
    yield assert_equal, sndd.unindent_string(param_doc3), param_doc1


def test_unindent_dict():
    d2 = sndd.unindent_dict(doc_dict)
    yield assert_equal, d2['strtest1'], doc_dict['strtest1']
    yield assert_equal, d2['strtest2'], doc_dict['strtest2']
    yield assert_equal, d2['strtest3'], doc_dict['strtest1']


def test_docformat():
    udd = sndd.unindent_dict(doc_dict)
    formatted = sndd.docformat(docstring, udd)
    yield assert_equal, formatted, filled_docstring
    single_doc = 'Single line doc %(strtest1)s'
    formatted = sndd.docformat(single_doc, doc_dict)
    # Note - initial indent of format string does not
    # affect subsequent indent of inserted parameter
    yield assert_equal, formatted, """Single line doc Another test
   with some indent"""


def test_decorator():
    # with unindentation of parameters
    decorator = sndd.filldoc(doc_dict, True)
    @decorator
    def func():
        """ Docstring
        %(strtest3)s
        """
    yield assert_equal, func.__doc__, """ Docstring
        Another test
           with some indent
        """
    # without unindentation of parameters
    decorator = sndd.filldoc(doc_dict, False)
    @decorator
    def func():
        """ Docstring
        %(strtest3)s
        """
    yield assert_equal, func.__doc__, """ Docstring
            Another test
               with some indent
        """
