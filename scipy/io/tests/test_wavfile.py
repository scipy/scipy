from __future__ import division, print_function, absolute_import

import os
import tempfile
import warnings

import numpy as np
from numpy.testing import assert_equal, assert_, assert_raises, assert_array_equal
from numpy.testing.utils import WarningManager
from scipy.io import wavfile


def datafile(fn):
    return os.path.join(os.path.dirname(__file__), 'data', fn)

def test_read_1():
    for mmap in [False, True]:
        warn_ctx = WarningManager()
        warn_ctx.__enter__()
        try:
            warnings.simplefilter('ignore', wavfile.WavFileWarning)
            rate, data = wavfile.read(datafile('test-44100-le-1ch-4bytes.wav'),
                                      mmap=mmap)
        finally:
            warn_ctx.__exit__()

        assert_equal(rate, 44100)
        assert_(np.issubdtype(data.dtype, np.int32))
        assert_equal(data.shape, (4410,))

        del data

def test_read_2():
    for mmap in [False, True]:
        rate, data = wavfile.read(datafile('test-8000-le-2ch-1byteu.wav'),
                                  mmap=mmap)
        assert_equal(rate, 8000)
        assert_(np.issubdtype(data.dtype, np.uint8))
        assert_equal(data.shape, (800, 2))

        del data

def test_read_fail():
    for mmap in [False, True]:
        fp = open(datafile('example_1.nc'))
        assert_raises(ValueError, wavfile.read, fp, mmap=mmap)
        fp.close()

def _check_roundtrip(rate, dtype, channels):
    fd, tmpfile = tempfile.mkstemp(suffix='.wav')
    try:
        os.close(fd)

        data = np.random.rand(100, channels)
        if channels == 1:
            data = data[:,0]
        data = (data*128).astype(dtype)

        wavfile.write(tmpfile, rate, data)

        for mmap in [False, True]:
            rate2, data2 = wavfile.read(tmpfile, mmap=mmap)

            assert_equal(rate, rate2)
            assert_(data2.dtype.byteorder in ('<', '=', '|'), msg=data2.dtype)
            assert_array_equal(data, data2)

            del data2
    finally:
        os.unlink(tmpfile)

def test_write_roundtrip():
    for signed in ('i', 'u'):
        for size in (1, 2, 4, 8):
            if size == 1 and signed == 'i':
                # signed 8-bit integer PCM is not allowed
                continue
            for endianness in ('>', '<'):
                if size == 1 and endianness == '<':
                    continue
                for rate in (8000, 32000):
                    for channels in (1, 2, 5):
                        dt = np.dtype('%s%s%d' % (endianness, signed, size))
                        yield _check_roundtrip, rate, dt, channels
