from __future__ import division, print_function, absolute_import

import os.path
import tempfile
import shutil
import numpy as np
import warnings

from numpy.testing import (assert_, assert_equal, dec, decorate_methods,
                           TestCase, run_module_suite, assert_allclose)

from scipy import misc

try:
    import PIL.Image
except ImportError:
    _have_PIL = False
else:
    _have_PIL = True


# Function / method decorator for skipping PIL tests on import failure
_pilskip = dec.skipif(not _have_PIL, 'Need to import PIL for this test')

datapath = os.path.dirname(__file__)


class TestPILUtil(TestCase):
    def test_imresize(self):
        im = np.random.random((10,20))
        for T in np.sctypes['float'] + [float]:
            # 1.1 rounds to below 1.1 for float16, 1.101 works
            im1 = misc.imresize(im,T(1.101))
            assert_equal(im1.shape,(11,22))

    def test_imresize2(self):
        im = np.random.random((20,30))
        im2 = misc.imresize(im, (30,40), interp='bicubic')
        assert_equal(im2.shape, (30,40))

    def test_imresize3(self):
        im = np.random.random((15,30))
        im2 = misc.imresize(im, (30,60), interp='nearest')
        assert_equal(im2.shape, (30,60))

    def test_imresize4(self):
        im = np.array([[1, 2],
                       [3, 4]])
        # Check that resizing by target size, float and int are the same
        im2 = misc.imresize(im, (4,4), mode='F')  # output size
        im3 = misc.imresize(im, 2., mode='F')  # fraction
        im4 = misc.imresize(im, 200, mode='F')  # percentage
        assert_equal(im2, im3)
        assert_equal(im2, im4)

    def test_bytescale(self):
        x = np.array([0,1,2], np.uint8)
        y = np.array([0,1,2])
        assert_equal(misc.bytescale(x), x)
        assert_equal(misc.bytescale(y), [0,127,255])

    def test_bytescale_keywords(self):
        x = np.array([40, 60, 120, 200, 300, 500])
        res_lowhigh = misc.bytescale(x, low=10, high=143)
        assert_equal(res_lowhigh, [10, 16, 33, 56, 85, 143])
        res_cmincmax = misc.bytescale(x, cmin=60, cmax=300)
        assert_equal(res_cmincmax, [0, 0, 64, 149, 255, 255])
        assert_equal(misc.bytescale(np.array([3, 3, 3]), low=4), [4, 4, 4])

    def test_imsave(self):
        with warnings.catch_warnings(record=True):  # PIL ResourceWarning
            img = misc.imread(os.path.join(datapath, 'data', 'icon.png'))
        tmpdir = tempfile.mkdtemp()
        try:
            fn1 = os.path.join(tmpdir, 'test.png')
            fn2 = os.path.join(tmpdir, 'testimg')
            with warnings.catch_warnings(record=True):  # PIL ResourceWarning
                misc.imsave(fn1, img)
                misc.imsave(fn2, img, 'PNG')

            with warnings.catch_warnings(record=True):  # PIL ResourceWarning
                data1 = misc.imread(fn1)
                data2 = misc.imread(fn2)

            assert_allclose(data1, img)
            assert_allclose(data2, img)
        finally:
            shutil.rmtree(tmpdir)

    def test_imwrite(self):
        with warnings.catch_warnings(record=True):  # PIL ResourceWarning
            img = misc.imread(os.path.join(datapath, 'data', 'icon.png'))
        tmpdir = tempfile.mkdtemp()
        try:
            fn1 = os.path.join(tmpdir, 'test.png')
            fn2 = os.path.join(tmpdir, 'testimg')
            with warnings.catch_warnings(record=True):  # PIL ResourceWarning
                misc.imwrite(fn1, img)
                misc.imwrite(fn2, img, 'PNG')

            with warnings.catch_warnings(record=True):  # PIL ResourceWarning
                data1 = misc.imread(fn1)
                data2 = misc.imread(fn2)

            assert_allclose(data1, img)
            assert_allclose(data2, img)
        finally:
            shutil.rmtree(tmpdir)


def tst_fromimage(filename, irange):
    fp = open(filename, "rb")
    img = misc.fromimage(PIL.Image.open(fp))
    fp.close()
    imin,imax = irange
    assert_(img.min() >= imin)
    assert_(img.max() <= imax)


@_pilskip
def test_fromimage():
    # Test generator for parametric tests
    data = {'icon.png':(0,255),
            'icon_mono.png':(0,2),
            'icon_mono_flat.png':(0,1)}
    for fn, irange in data.items():
        yield tst_fromimage, os.path.join(datapath,'data',fn), irange

decorate_methods(TestPILUtil, _pilskip)


if __name__ == "__main__":
    run_module_suite()
