""" Test reading of files not conforming to matlab specification

We try and read any file that matlab reads, these files included
"""
import struct
from io import BytesIO
from os.path import dirname, join as pjoin

import numpy as np
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


def test_malformed1():
    # Example from gh-6072
    # Contains malformed header data, which previously resulted into a
    # buffer overflow.
    #
    # Should raise an exception, not segfault
    fname = pjoin(TEST_DATA_PATH, 'malformed1.mat')
    with open(fname, 'rb') as f:
        assert_raises(ValueError, loadmat, f)


def _patch_mdtype(payload, marker, mdtype):
    # Overwrite the 4-byte mdtype of the full element whose data is `marker`.
    data = bytearray(payload)
    idx = data.find(marker)
    assert idx != -1
    struct.pack_into('<I', data, idx - 8, mdtype)
    return bytes(data)


def test_bad_data_type_code():
    # A data element tag carries the mdtype straight from the file. An unknown
    # or out-of-range code must not be used to index the dtypes array.
    #
    # Should raise an exception, not segfault
    num = BytesIO()
    savemat(num, {'a': np.array([[1234.5]])}, do_compression=False)
    txt = BytesIO()
    savemat(txt, {'s': 'WXYZ0123'}, do_compression=False)
    for payload, marker in ((num.getvalue(), struct.pack('<d', 1234.5)),
                            (txt.getvalue(), b'WXYZ0123')):
        for code in (8, 9999):  # unused in-range slot, then past the array
            bad = _patch_mdtype(payload, marker, code)
            assert_raises(ValueError, loadmat, BytesIO(bad))
