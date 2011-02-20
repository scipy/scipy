import os.path
import numpy as np

from numpy.testing import assert_, assert_equal, \
        dec, decorate_methods, TestCase, run_module_suite

try:
    import PIL.Image
except ImportError:
    _have_PIL = False
else:
    _have_PIL = True
    import scipy.misc.pilutil as pilutil

# Function / method decorator for skipping PIL tests on import failure
_pilskip = dec.skipif(not _have_PIL, 'Need to import PIL for this test')

datapath = os.path.dirname(__file__)

class TestPILUtil(TestCase):
    def test_imresize(self):
        im = np.random.random((10,20))
        for T in np.sctypes['float'] + [float]:
            # 1.1 rounds to below 1.1 for float16, 1.101 works
            im1 = pilutil.imresize(im,T(1.101))
            assert_equal(im1.shape,(11,22))

    def test_imresize2(self):
        im = np.random.random((20,30))
        im2 = pilutil.imresize(im, (30,40), interp='bicubic')
        assert_equal(im2.shape, (30,40))

    def test_imresize3(self):
        im = np.random.random((15,30))
        im2 = pilutil.imresize(im, (30,60), interp='nearest')
        assert_equal(im2.shape, (30,60))

    def test_bytescale(self):
        x = np.array([0,1,2],np.uint8)
        y = np.array([0,1,2])
        assert_equal(pilutil.bytescale(x),x)
        assert_equal(pilutil.bytescale(y),[0,127,255])

def tst_fromimage(filename, irange):
    fp = open(filename, "rb")
    img = pilutil.fromimage(PIL.Image.open(fp))
    fp.close()
    imin,imax = irange
    assert_(img.min() >= imin)
    assert_(img.max() <= imax)

@_pilskip
def test_fromimage():
    ''' Test generator for parametric tests '''
    data = {'icon.png':(0,255),
            'icon_mono.png':(0,2),
            'icon_mono_flat.png':(0,1)}
    for fn, irange in data.iteritems():
        yield tst_fromimage, os.path.join(datapath,'data',fn), irange

decorate_methods(TestPILUtil, _pilskip)

if __name__ == "__main__":
    run_module_suite()
