""" Testing 

"""

import os

import StringIO
import cStringIO
from tempfile import mkstemp

import numpy as np

from nose.tools import assert_true, assert_false, \
     assert_equal, assert_raises

from numpy.testing import assert_array_equal, assert_array_almost_equal

from scipy.io.matlab.streams import make_stream, \
    GenericStream, cStringStream, FileStream, \
    _read_into, _read_string


def setup():
    val = 'a\x00string'
    global fs, gs, cs, fname
    fd, fname = mkstemp()
    fs = os.fdopen(fd, 'wb')
    fs.write(val)
    fs.close()
    fs = open(fname)
    gs = StringIO.StringIO(val)
    cs = cStringIO.StringIO(val)


def teardown():
    global fname, fs
    del fs
    os.unlink(fname)


def test_make_stream():
    global fs, gs, cs
    # test stream initialization
    yield assert_true, isinstance(make_stream(gs), GenericStream)
    yield assert_true, isinstance(make_stream(cs), cStringStream)
    yield assert_true, isinstance(make_stream(fs), FileStream)


def test_tell_seek():
    global fs, gs, cs
    for s in (fs, gs, cs):
        st = make_stream(s)
        res = st.seek(0)
        yield assert_equal, res, 0
        yield assert_equal, st.tell(), 0
        res = st.seek(5)
        yield assert_equal, res, 0
        yield assert_equal, st.tell(), 5
        res = st.seek(2, 1)
        yield assert_equal, res, 0
        yield assert_equal, st.tell(), 7
        res = st.seek(-2, 2)
        yield assert_equal, res, 0
        yield assert_equal, st.tell(), 6


def test_read():
    global fs, gs, cs
    for s in (fs, gs, cs):
        st = make_stream(s)
        st.seek(0)
        res = st.read(-1)
        yield assert_equal, res, 'a\x00string'
        st.seek(0)
        res = st.read(4)
        yield assert_equal, res, 'a\x00st'
        # read into
        st.seek(0)
        res = _read_into(st, 4)
        yield assert_equal, res, 'a\x00st'
        res = _read_into(st, 4)
        yield assert_equal, res, 'ring'
        yield assert_raises, IOError, _read_into, st, 2
        # read alloc
        st.seek(0)
        res = _read_string(st, 4)
        yield assert_equal, res, 'a\x00st'
        res = _read_string(st, 4)
        yield assert_equal, res, 'ring'
        yield assert_raises, IOError, _read_string, st, 2
