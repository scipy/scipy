""" Testing

"""

from __future__ import division, print_function, absolute_import

import os
import sys
import zlib

from io import BytesIO

if sys.version_info[0] >= 3:
    cStringIO = BytesIO
else:
    from cStringIO import StringIO as cStringIO

from tempfile import mkstemp

import numpy as np

from nose.tools import assert_true, assert_false, \
     assert_equal, assert_raises

from numpy.testing import assert_array_equal, assert_array_almost_equal, \
     run_module_suite

from scipy.io.matlab.streams import make_stream, \
    GenericStream, cStringStream, FileStream, ZlibInputStream, \
    _read_into, _read_string


fs = None
gs = None
cs = None
fname = None


def setup():
    val = b'a\x00string'
    global fs, gs, cs, fname
    fd, fname = mkstemp()
    fs = os.fdopen(fd, 'wb')
    fs.write(val)
    fs.close()
    fs = open(fname, 'rb')
    gs = BytesIO(val)
    cs = cStringIO(val)


def teardown():
    global fname, fs
    fs.close()
    del fs
    os.unlink(fname)


def test_make_stream():
    global fs, gs, cs
    # test stream initialization
    assert_true(isinstance(make_stream(gs), GenericStream))
    if sys.version_info[0] < 3:
        assert_true(isinstance(make_stream(cs), cStringStream))
        assert_true(isinstance(make_stream(fs), FileStream))


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
        yield assert_equal, res, b'a\x00string'
        st.seek(0)
        res = st.read(4)
        yield assert_equal, res, b'a\x00st'
        # read into
        st.seek(0)
        res = _read_into(st, 4)
        yield assert_equal, res, b'a\x00st'
        res = _read_into(st, 4)
        yield assert_equal, res, b'ring'
        yield assert_raises, IOError, _read_into, st, 2
        # read alloc
        st.seek(0)
        res = _read_string(st, 4)
        yield assert_equal, res, b'a\x00st'
        res = _read_string(st, 4)
        yield assert_equal, res, b'ring'
        yield assert_raises, IOError, _read_string, st, 2


class TestZlibInputStream(object):
    def _get_data(self, size):
        data = np.random.randint(0, 256, size).astype(np.uint8).tostring()
        compressed_data = zlib.compress(data)
        stream = BytesIO(compressed_data)
        return stream, len(compressed_data), data

    def test_read(self):
        block_size = 131072

        SIZES = [0, 1, 10, block_size//2, block_size-1,
                 block_size, block_size+1, 2*block_size-1]

        READ_SIZES = [block_size//2, block_size-1,
                      block_size, block_size+1]

        def check(size, read_size):
            compressed_stream, compressed_data_len, data = self._get_data(size)
            stream = ZlibInputStream(compressed_stream, compressed_data_len)
            data2 = b''
            so_far = 0
            while True:
                block = stream.read(min(read_size,
                                        size - so_far))
                if not block:
                    break
                so_far += len(block)
                data2 += block
            assert_equal(data, data2)

        for size in SIZES:
            for read_size in READ_SIZES:
                yield check, size, read_size

    def test_read_max_length(self):
        size = 1234
        data = np.random.randint(0, 256, size).astype(np.uint8).tostring()
        compressed_data = zlib.compress(data)
        compressed_stream = BytesIO(compressed_data + b"abbacaca")
        stream = ZlibInputStream(compressed_stream, len(compressed_data))

        stream.read(len(data))
        assert_equal(compressed_stream.tell(), len(compressed_data))

        assert_equal(stream.read(1), b"")

    def test_seek(self):
        compressed_stream, compressed_data_len, data = self._get_data(1024)

        stream = ZlibInputStream(compressed_stream, compressed_data_len)

        stream.seek(123)
        p = 123
        assert_equal(stream.tell(), p)
        d1 = stream.read(11)
        assert_equal(d1, data[p:p+11])

        stream.seek(321, 1)
        p = 123+11+321
        assert_equal(stream.tell(), p)
        d2 = stream.read(21)
        assert_equal(d2, data[p:p+21])

        stream.seek(641, 0)
        p = 641
        assert_equal(stream.tell(), p)
        d3 = stream.read(11)
        assert_equal(d3, data[p:p+11])

        assert_raises(IOError, stream.seek, 10, 2)
        assert_raises(IOError, stream.seek, -1, 1)
        assert_raises(ValueError, stream.seek, 1, 123)

        stream.seek(10000, 1)
        assert_equal(stream.read(12), b"")

    def test_all_data_read(self):
        compressed_stream, compressed_data_len, data = self._get_data(1024)
        stream = ZlibInputStream(compressed_stream, compressed_data_len)
        assert_false(stream.all_data_read())
        stream.seek(512)
        assert_false(stream.all_data_read())
        stream.seek(1024)
        assert_true(stream.all_data_read())


if __name__ == "__main__":
    run_module_suite()
