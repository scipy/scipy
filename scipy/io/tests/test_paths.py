"""
Ensure that we can use pathlib.Path objects in all relevant IO functions.
"""
import os
import sys

try:
    from pathlib import Path
except ImportError:
    # Not available. No fallback import, since we'll skip the entire
    # test suite for Python < 3.6.
    pass

import numpy as np
from numpy.testing import assert_, assert_raises
import pytest

import scipy.io
from scipy._lib._tmpdirs import tempdir
import scipy.sparse

# Bit of a hack to keep the test runner from exploding in Python 2.7.
# FileNotFoundError was added in Python 3.3.
if sys.version_info < (3, 3):
    FileNotFoundError = IOError


@pytest.mark.skipif(sys.version_info < (3, 6),
                    reason='Passing path-like objects to IO functions requires Python >= 3.6')
class TestPaths(object):
    data = np.arange(5)

    def test_savemat(self):
        with tempdir() as temp_dir:
            path = Path(temp_dir) / 'data.mat'
            scipy.io.savemat(path, {'data': self.data})
            assert_(path.is_file())

    def test_loadmat(self):
        # Save data with string path, load with pathlib.Path
        with tempdir() as temp_dir:
            path = Path(temp_dir) / 'data.mat'
            scipy.io.savemat(str(path), {'data': self.data})

            mat_contents = scipy.io.loadmat(path)
            assert_((mat_contents['data'] == self.data).all())

    def test_whosmat(self):
        # Save data with string path, load with pathlib.Path
        with tempdir() as temp_dir:
            path = Path(temp_dir) / 'data.mat'
            scipy.io.savemat(str(path), {'data': self.data})

            contents = scipy.io.whosmat(path)
            assert_(contents[0] == ('data', (1, 5), 'int64'))

    def test_readsav(self):
        filename = os.path.join(os.path.dirname(__file__), 'data', 'scalar_string.sav')
        path = Path(filename)
        scipy.io.readsav(path)

    def test_hb_read(self):
        # Save data with string path, load with pathlib.Path
        with tempdir() as temp_dir:
            data = scipy.sparse.csr_matrix(scipy.sparse.eye(3))
            path = Path(temp_dir) / 'data.hb'
            scipy.io.harwell_boeing.hb_write(str(path), data)

            data_new = scipy.io.harwell_boeing.hb_read(path)
            assert_((data_new != data).nnz == 0)

    def test_hb_write(self):
        with tempdir() as temp_dir:
            data = scipy.sparse.csr_matrix(scipy.sparse.eye(3))
            path = Path(temp_dir) / 'data.hb'
            scipy.io.harwell_boeing.hb_write(path, data)
            assert_(path.is_file())
