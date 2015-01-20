#!/usr/bin/env python
''' Jottings to work out format for __function_workspace__ matrix at end
of mat file.

'''
from __future__ import division, print_function, absolute_import

from os.path import join as pjoin, dirname
import sys

from io import BytesIO

from numpy.testing import \
     assert_array_equal, \
     assert_array_almost_equal, \
     assert_equal, \
     assert_raises, run_module_suite

from nose.tools import assert_true

import numpy as np
from numpy.compat import asstr

from scipy.io.matlab.mio5 import MatlabObject, MatFile5Writer, \
      MatFile5Reader, MatlabFunction

test_data_path = pjoin(dirname(__file__), 'data')


def read_workspace_vars(fname):
    fp = open(fname, 'rb')
    rdr = MatFile5Reader(fp, struct_as_record=True)
    return rdr.minimat_reader().get_variables()


def test_jottings():
    # example
    fname = pjoin(test_data_path, 'parabola.mat')
    ws_vars = read_workspace_vars(fname)


if __name__ == "__main__":
    run_module_suite(argv=sys.argv)
