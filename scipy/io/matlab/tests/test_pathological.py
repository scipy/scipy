""" Test reading of files not conforming to matlab specification

We try and read any file that matlab reads, these files included
"""
from os.path import dirname, join as pjoin
import sys

if sys.version_info[0] >= 3:
    from io import BytesIO
    cStringIO = BytesIO
else:
    from cStringIO import StringIO as cStringIO
    from StringIO import StringIO as BytesIO

import numpy as np

from nose.tools import assert_true, assert_false, \
     assert_equal, assert_raises

from numpy.testing import assert_array_equal, assert_array_almost_equal, \
     run_module_suite

from scipy.io.matlab.mio import loadmat

TEST_DATA_PATH = pjoin(dirname(__file__), 'data')

def test_multiple_fieldnames():
    # Example provided by Dharhas Pothina
    # Extracted using mio5.varmats_from_mat
    multi_fname = pjoin(TEST_DATA_PATH, 'nasty_duplicate_fieldnames.mat')
    vars = loadmat(multi_fname)
    funny_names = vars['Summary'].dtype.names
    assert_true(set(['_1_Station_Q', '_2_Station_Q',
                     '_3_Station_Q']).issubset(funny_names))
