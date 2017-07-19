"""
Ensure that we can use pathlib.Path objects in all relevant IO functions.
"""
import sys
import unittest

try:
    from pathlib import Path
except ImportError:
    # Not available. No fallback import, since we'll skip the entire
    # test suite for Python < 3.6.
    pass

import numpy as np
from numpy.testing import assert_, assert_raises

import scipy.io
from scipy._lib._tmpdirs import tempdir

@unittest.skipIf(sys.version_info < (3, 6),
                 'Passing path-like objects to IO functions requires Python >= 3.6')
class TestPaths(unittest.TestCase):
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

            data = scipy.io.loadmat(path)
            assert_((data == self.data).all())

    def test_whosmat(self):
        # Save data with string path, load with pathlib.Path
        with tempdir() as temp_dir:
            path = Path(temp_dir) / 'data.mat'
            scipy.io.savemat(str(path), {'data': self.data})

            contents = scipy.io.whosmat(path)
            assert_(contents[0] == ('data', (1, 5), 'int64'))

    def test_readsav(self):
        # I have no idea how to create a valid .sav file from IDL, so let's
        # just try to open something that doesn't exist, and ensure that we get
        # a FileNotFoundError (as opposed to a TypeError, if `open` can't
        # understand the argument, or another type of exception if `readsav` is
        # unnecessarily strict in which type it accepts.

        # FileNotFoundError was added in Python 3.3, and this would throw an
        # IOError in older Python, but this entire test class is skipped in
        # Python < 3.6 so it's safe to use the newer exception class here.
        with self.assertRaises(FileNotFoundError):
            path = Path('nonexistent.sav')
            scipy.io.readsav(path)
