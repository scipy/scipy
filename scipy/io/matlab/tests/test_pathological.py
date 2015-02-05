""" Test reading of files not conforming to matlab specification

We try and read any file that matlab reads, these files included
"""
from __future__ import division, print_function, absolute_import

from os.path import dirname, join as pjoin
import sys

from nose.tools import assert_true

from numpy.testing import run_module_suite

from scipy.io.matlab.mio import loadmat

TEST_DATA_PATH = pjoin(dirname(__file__), 'data')


def test_multiple_fieldnames():
    # Example provided by Dharhas Pothina
    # Extracted using mio5.varmats_from_mat
    multi_fname = pjoin(TEST_DATA_PATH, 'nasty_duplicate_fieldnames.mat')
    vars = loadmat(multi_fname)
    funny_names = set(vars['Summary'].dtype.names)
    assert_true(set(['_1_Station_Q', '_2_Station_Q', '_3_Station_Q']) < funny_names)


if __name__ == "__main__":
    run_module_suite()
