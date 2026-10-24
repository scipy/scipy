""" Test reading of files not conforming to matlab specification

We try and read any file that matlab reads, these files included
"""
from io import BytesIO
from os.path import dirname, join as pjoin

import pytest
from numpy.testing import assert_
from pytest import raises as assert_raises

from scipy.io.matlab._mio import loadmat, savemat

TEST_DATA_PATH = pjoin(dirname(__file__), 'data')


def test_multiple_fieldnames():
    # Example provided by Dharhas Pothina
    # Extracted using mio5.varmats_from_mat
    multi_fname = pjoin(TEST_DATA_PATH, 'nasty_duplicate_fieldnames.mat')
    vars = loadmat(multi_fname)
    funny_names = vars['Summary'].dtype.names
    assert_({'_1_Station_Q', '_2_Station_Q',
                     '_3_Station_Q'}.issubset(funny_names))


# savemat is not thread-safe (class-level VarWriter5.mat_tag), as in test_mio.py
@pytest.mark.thread_unsafe
def test_duplicate_fieldnames_numbering():
    # Repeats of a field name are numbered per name, in order of appearance.
    # savemat cannot write duplicates, so patch unique names in the file.
    stream = BytesIO()
    savemat(stream, {'s': {'a0': 0, 'b0': 1, 'a1': 2, 'a2': 3, 'b1': 4, 'c0': 5}})
    contents = stream.getvalue()
    for name in (b'a0', b'b0', b'a1', b'a2', b'b1', b'c0'):
        assert contents.count(name + b'\x00') == 1
        contents = contents.replace(name + b'\x00', name[:1] + b'\x00\x00')
    s = loadmat(BytesIO(contents), squeeze_me=True)['s']
    assert s.dtype.names == ('a', 'b', '_1_a', '_2_a', '_1_b', 'c')
    assert [s[name] for name in s.dtype.names] == [0, 1, 2, 3, 4, 5]


def test_malformed1():
    # Example from gh-6072
    # Contains malformed header data, which previously resulted into a
    # buffer overflow.
    #
    # Should raise an exception, not segfault
    fname = pjoin(TEST_DATA_PATH, 'malformed1.mat')
    with open(fname, 'rb') as f:
        assert_raises(ValueError, loadmat, f)
