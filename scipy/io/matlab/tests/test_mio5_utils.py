""" Testing 

"""
import cStringIO
import StringIO

import numpy as np

from nose.tools import assert_true, assert_false, \
     assert_equal, assert_raises

from numpy.testing import assert_array_equal, assert_array_almost_equal

import scipy.io.matlab.byteordercodes as boc
import scipy.io.matlab.streams as streams
import scipy.io.matlab.miobase as miob
import scipy.io.matlab.mio5 as mio5
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
        padding = 8 - (byte_count % 8)
        all_dt = [('mdtype', udt),
                  ('byte_count', udt),
                  ('val', base_dt)]
        if padding:
            all_dt.append(('padding', 'u1', padding))
    else: # is sde
        udt = bo + 'u2'
        padding = 4-byte_count
        if bo == '<': # little endian
            all_dt = [('mdtype', udt),
                      ('byte_count', udt),
                      ('val', base_dt)]
        else: # big endian
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
    for s in strings:
        stream.write(s)
    stream.seek(0)


def _make_readerlike():
    class R(object):
        pass
    r = R()
    r.byte_order = boc.native_code
    r.dtypes = {}
    r.class_dtypes = {}
    r.codecs = {}
    r.struct_as_record = True
    r.uint16_codec = None
    r.chars_as_strings = False
    r.mat_dtype = False
    r.squeeze_me = False
    return r


def test_read_tag():
    # mainly to test errors
    # make reader-like thing
    str_io = StringIO.StringIO()
    r = _make_readerlike()
    r.mat_stream = str_io
    c_reader = m5u.VarReader5(r)
    # This works for StringIO but _not_ cStringIO
    yield assert_raises, IOError, c_reader.read_tag
    # bad SDE
    tag = _make_tag('i4', 1, mio5.miINT32, sde=True)
    tag['byte_count'] = 5
    _write_stream(str_io, tag.tostring())
    yield assert_raises, ValueError, c_reader.read_tag


def test_read_stream():
    tag = _make_tag('i4', 1, mio5.miINT32, sde=True)
    tag_str = tag.tostring()
    str_io = cStringIO.StringIO(tag_str)
    st = streams.make_stream(str_io)
    s = streams._read_into(st, tag.itemsize)
    yield assert_equal, s, tag.tostring()


def test_read_numeric():
    # make reader-like thing
    str_io = cStringIO.StringIO()
    r = _make_readerlike()
    r.mat_stream = str_io
    # check simplest of tags
    for base_dt, val, mdtype in (
        ('u2', 30, mio5.miUINT16),
        ('i4', 1, mio5.miINT32),
        ('i2', -1, mio5.miINT16)):
        for byte_code in ('<', '>'):
            r.byte_order = byte_code
            r.dtypes = miob.convert_dtypes(mio5.mdtypes_template, byte_code)
            c_reader = m5u.VarReader5(r)
            yield assert_equal, c_reader.little_endian, byte_code == '<'
            yield assert_equal, c_reader.is_swapped, byte_code != boc.native_code
            for sde_f in (False, True):
                dt = np.dtype(base_dt).newbyteorder(byte_code)
                a = _make_tag(dt, val, mdtype, sde_f)
                a_str = a.tostring()
                _write_stream(str_io, a_str)
                el = c_reader.read_numeric()
                yield assert_equal, el, val
                # two sequential reads
                _write_stream(str_io, a_str, a_str)
                el = c_reader.read_numeric()
                yield assert_equal, el, val
                el = c_reader.read_numeric()
                yield assert_equal, el, val
    

def test_read_numeric_writeable():
    # make reader-like thing
    str_io = cStringIO.StringIO()
    r = _make_readerlike()
    r.mat_stream = str_io
    r.byte_order = '<'
    r.dtypes = miob.convert_dtypes(mio5.mdtypes_template, '<')
    c_reader = m5u.VarReader5(r)
    dt = np.dtype('<u2')
    a = _make_tag(dt, 30, mio5.miUINT16, 0)
    a_str = a.tostring()
    _write_stream(str_io, a_str)
    el = c_reader.read_numeric()
    yield assert_true, el.flags.writeable
    
