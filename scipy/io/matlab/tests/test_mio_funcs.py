#!/usr/bin/env python
''' Jottings to work out format for __function_workspace__ matrix at end
of mat file.

'''
from os.path import join as pjoin, dirname
import sys

if sys.version_info[0] >= 3:
    from io import BytesIO
else:
    from cStringIO import StringIO as BytesIO

from numpy.testing import \
     assert_array_equal, \
     assert_array_almost_equal, \
     assert_equal, \
     assert_raises, run_module_suite

from nose.tools import assert_true

import numpy as np
from numpy.compat import asbytes, asstr

from scipy.io.matlab.mio5 import MatlabObject, MatFile5Writer, \
      MatFile5Reader, MatlabFunction

test_data_path = pjoin(dirname(__file__), 'data')

def read_minimat_vars(rdr):
    rdr.initialize_read()
    mdict = {'__globals__': []}
    i = 0
    while not rdr.end_of_stream():
        hdr, next_position = rdr.read_var_header()
        name = asstr(hdr.name)
        if name == '':
            name = 'var_%d' % i
            i += 1
        res = rdr.read_var_array(hdr, process=False)
        rdr.mat_stream.seek(next_position)
        mdict[name] = res
        if hdr.is_global:
            mdict['__globals__'].append(name)
    return mdict

def read_workspace_vars(fname):
    fp = open(fname, 'rb')
    rdr = MatFile5Reader(fp, struct_as_record=True)
    vars = rdr.get_variables()
    fws = vars['__function_workspace__']
    ws_bs = BytesIO(fws.tostring())
    ws_bs.seek(2)
    rdr.mat_stream = ws_bs
    # Guess byte order.
    mi = rdr.mat_stream.read(2)
    rdr.byte_order = mi == asbytes('IM') and '<' or '>'
    rdr.mat_stream.read(4) # presumably byte padding
    mdict = read_minimat_vars(rdr)
    fp.close()
    return mdict


def test_jottings():
    # example
    fname = pjoin(test_data_path, 'parabola.mat')
    ws_vars = read_workspace_vars(fname)

if __name__ == "__main__":
    run_module_suite()
