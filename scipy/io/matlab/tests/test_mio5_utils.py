""" Testing

"""
from __future__ import division, print_function, absolute_import

import sys

from io import BytesIO
cStringIO = BytesIO

import numpy as np

from nose.tools import assert_true, assert_false, \
     assert_equal, assert_raises

from numpy.testing import assert_array_equal, assert_array_almost_equal, \
     run_module_suite

from scipy.lib.six import u

import scipy.io.matlab.byteordercodes as boc
import scipy.io.matlab.streams as streams
from scipy.io.matlab.mio5 import MatFile5Reader
import scipy.io.matlab.mio5_params as mio5p
import scipy.io.matlab.mio5_utils as m5u


def test_byteswap():
    for val in (
        1,
        0x100,
        0x10000):
        a = np.array(val, dtype=np.uint32)
        b = a.byteswap()
        c = m5u.byteswap_u4(a)
        yield assert_equal, b.item(), c
        d = m5u.byteswap_u4(c)
        yield assert_equal, a.item(), d


def _make_tag(base_dt, val, mdtype, sde=False):
    ''' Makes a simple matlab tag, full or sde '''
    base_dt = np.dtype(base_dt)
    bo = boc.to_numpy_code(base_dt.byteorder)
    byte_count = base_dt.itemsize
    if not sde:
        udt = bo + 'u4'
        padding = (-byte_count) % 8
        all_dt = [('mdtype', udt),
                  ('byte_count', udt),
                  ('val', base_dt)]
        if padding:
            all_dt.append(('padding', 'u1', padding))
    else:  # is sde
        udt = bo + 'u2'
        padding = 4 - byte_count
        if bo == '<':  # little endian
            all_dt = [('mdtype', udt),
                      ('byte_count', udt),
                      ('val', base_dt)]
        else:  # big endian
            all_dt = [('byte_count', udt),
                      ('mdtype', udt),
                      ('val', base_dt)]
        if padding:
            all_dt.append(('padding', 'u1', padding))
    tag = np.zeros((1,), dtype=all_dt)
    tag['mdtype'] = mdtype
    tag['byte_count'] = byte_count
    tag['val'] = val
    return tag


def _write_stream(stream, *strings):
    stream.truncate(0)
    stream.seek(0)
    for s in strings:
        stream.write(s)
    stream.seek(0)


def _make_readerlike(stream, byte_order=boc.native_code):
    class R(object):
        pass
    r = object.__new__(MatFile5Reader)
    r._stream = stream
    r._endian = byte_order
    r._struct_as_record = True
    r._chars_as_strings = False
    r._mat_dtype = False
    r._squeeze_me = False
    return r


def test_read_tag():
    # mainly to test errors
    # make reader-like thing
    str_io = BytesIO()
    r = _make_readerlike(str_io)
    # This works for StringIO but _not_ cStringIO
    yield assert_raises, StopIteration, next, r._read_iter()
    # bad SDE
    tag = _make_tag('i4', 1, mio5p.miINT32, sde=True)
    tag['byte_count'] = 5
    _write_stream(str_io, tag.tostring())
    yield assert_raises, ValueError, next, r._read_iter()


def test_read_stream():
    tag = _make_tag('i4', 1, mio5p.miINT32, sde=True)
    tag_str = tag.tostring()
    str_io = cStringIO(tag_str)
    st = streams.make_stream(str_io)
    s = streams._read_into(st, tag.itemsize)
    yield assert_equal, s, tag.tostring()


def test_read_numeric():
    # make reader-like thing
    str_io = cStringIO()
    r = _make_readerlike(str_io)
    # check simplest of tags
    for base_dt, val, mdtype in (('u2', 30, mio5p.miUINT16),
                                 ('i4', 1, mio5p.miINT32),
                                 ('i2', -1, mio5p.miINT16)):
        for byte_code in ('<', '>'):
            r._endian = byte_code
            for sde_f in (False, True):
                dt = np.dtype(base_dt).newbyteorder(byte_code)
                a = _make_tag(dt, val, mdtype, sde_f)
                a_str = a.tostring()
                _write_stream(str_io, a_str)
                el = next(r._read_iter())
                yield assert_equal, el, val
                # two sequential reads
                _write_stream(str_io, a_str, a_str)
                el = next(r._read_iter())
                yield assert_equal, el, val
                el = next(r._read_iter())
                yield assert_equal, el, val


def test_read_numeric_writeable():
    # make reader-like thing
    str_io = cStringIO()
    r = _make_readerlike(str_io, '<')
    dt = np.dtype('<u2')
    a = _make_tag(dt, 30, mio5p.miUINT16, sde=False)
    a_str = a.tostring()
    _write_stream(str_io, a_str)
    el = next(r._read_iter())
    yield assert_true, el.flags.writeable


# FIXME read_char is not available independently anymore.
#def test_zero_byte_string():
    ## Tests hack to allow chars of non-zero length, but 0 bytes
    ## make reader-like thing
    #str_io = cStringIO()
    #r = _make_readerlike(str_io, boc.native_code)
    #tag_dt = np.dtype([('mdtype', 'u4'), ('byte_count', 'u4')])
    #tag = np.zeros((1,), dtype=tag_dt)
    #tag['mdtype'] = mio5p.miINT8
    #tag['byte_count'] = 1
    #hdr = m5u.VarHeader5()
    ## Try when string is 1 length
    #hdr.set_dims([1,])
    #_write_stream(str_io, tag.tostring() + b'        ')
    #str_io.seek(0)
    #_, val = r._read1()
    #assert_equal(val, u(' '))
    ## Now when string has 0 bytes 1 length
    #tag['byte_count'] = 0
    #_write_stream(str_io, tag.tostring())
    #str_io.seek(0)
    #_, val = r._read1()
    #assert_equal(val, u(' '))
    ## Now when string has 0 bytes 4 length
    #str_io.seek(0)
    #hdr.set_dims([4,])
    #_, val = r._read1()
    #assert_array_equal(val, [u(' ')] * 4)


if __name__ == "__main__":
    run_module_suite(argv=sys.argv)
