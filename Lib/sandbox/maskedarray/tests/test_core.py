# pylint: disable-msg=W0611, W0612, W0511,R0201
"""Tests suite for MaskedArray.
Adapted from the original test_ma by Pierre Gerard-Marchant

:author: Pierre Gerard-Marchant
:contact: pierregm_at_uga_dot_edu
:version: $Id$
"""
__author__ = "Pierre GF Gerard-Marchant ($Author$)"
__version__ = '1.0'
__revision__ = "$Revision$"
__date__     = '$Date$'

import types

import numpy as N
import numpy.core.fromnumeric  as fromnumeric
from numpy.testing import NumpyTest, NumpyTestCase
from numpy.testing.utils import build_err_msg

import maskedarray.testutils
reload(maskedarray.testutils)
from maskedarray.testutils import *

import maskedarray.core
reload(maskedarray.core)
from maskedarray.core import *

pi = N.pi

#..............................................................................
class test_ma(NumpyTestCase):
    "Base test class for MaskedArrays."
    def __init__(self, *args, **kwds):
        NumpyTestCase.__init__(self, *args, **kwds)
        self.setUp()

    def setUp (self):
        "Base data definition."
        x = N.array([1.,1.,1.,-2., pi/2.0, 4., 5., -10., 10., 1., 2., 3.])
        y = N.array([5.,0.,3., 2., -1., -4., 0., -10., 10., 1., 0., 3.])
        a10 = 10.
        m1 = [1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0]
        m2 = [0, 0, 1, 0, 0, 1, 1, 0, 0, 0 ,0, 1]
        xm = masked_array(x, mask=m1)
        ym = masked_array(y, mask=m2)
        z = N.array([-.5, 0., .5, .8])
        zm = masked_array(z, mask=[0,1,0,0])
        xf = N.where(m1, 1.e+20, x)
        xm.set_fill_value(1.e+20)  
        self.d = (x, y, a10, m1, m2, xm, ym, z, zm, xf)
    #........................
    def check_testBasic1d(self):
        "Test of basic array creation and properties in 1 dimension."
        (x, y, a10, m1, m2, xm, ym, z, zm, xf) = self.d
        assert(not isMaskedArray(x))
        assert(isMaskedArray(xm))
        assert((xm-ym).filled(0).any())
        fail_if_equal(xm.mask.astype(int_), ym.mask.astype(int_))
        s = x.shape
        assert_equal(N.shape(xm), s)
        assert_equal(xm.shape, s)
        assert_equal(xm.dtype, x.dtype)
        assert_equal(zm.dtype, z.dtype)
        assert_equal(xm.size , reduce(lambda x,y:x*y, s))
        assert_equal(count(xm) , len(m1) - reduce(lambda x,y:x+y, m1))
        assert_array_equal(xm, xf)
        assert_array_equal(filled(xm, 1.e20), xf)
        assert_array_equal(x, xm)   
    #........................
    def check_testBasic2d(self):
        "Test of basic array creation and properties in 2 dimensions."
        (x, y, a10, m1, m2, xm, ym, z, zm, xf) = self.d
        for s in [(4,3), (6,2)]:
            x.shape = s
            y.shape = s
            xm.shape = s
            ym.shape = s
            xf.shape = s

            assert(not isMaskedArray(x))
            assert(isMaskedArray(xm))
            assert_equal(shape(xm), s)
            assert_equal(xm.shape, s)
            assert_equal( xm.size , reduce(lambda x,y:x*y, s))
            assert_equal( count(xm) , len(m1) - reduce(lambda x,y:x+y, m1))
            assert_equal(xm, xf)
            assert_equal(filled(xm, 1.e20), xf)
            assert_equal(x, xm)
    #........................
    def check_testArithmetic (self):
        "Test of basic arithmetic."
        (x, y, a10, m1, m2, xm, ym, z, zm, xf) = self.d
        a2d = array([[1,2],[0,4]])
        a2dm = masked_array(a2d, [[0,0],[1,0]])
        assert_equal(a2d * a2d, a2d * a2dm)
        assert_equal(a2d + a2d, a2d + a2dm)
        assert_equal(a2d - a2d, a2d - a2dm)
        for s in [(12,), (4,3), (2,6)]:
            x = x.reshape(s)
            y = y.reshape(s)
            xm = xm.reshape(s)
            ym = ym.reshape(s)
            xf = xf.reshape(s)
            assert_equal(-x, -xm)
            assert_equal(x + y, xm + ym)
            assert_equal(x - y, xm - ym)
            assert_equal(x * y, xm * ym)
            assert_equal(x / y, xm / ym)
            assert_equal(a10 + y, a10 + ym)
            assert_equal(a10 - y, a10 - ym)
            assert_equal(a10 * y, a10 * ym)
            assert_equal(a10 / y, a10 / ym)
            assert_equal(x + a10, xm + a10)
            assert_equal(x - a10, xm - a10)
            assert_equal(x * a10, xm * a10)
            assert_equal(x / a10, xm / a10)
            assert_equal(x**2, xm**2)
            assert_equal(abs(x)**2.5, abs(xm) **2.5)
            assert_equal(x**y, xm**ym)
            assert_equal(N.add(x,y), add(xm, ym))
            assert_equal(N.subtract(x,y), subtract(xm, ym))
            assert_equal(N.multiply(x,y), multiply(xm, ym))
            assert_equal(N.divide(x,y), divide(xm, ym))
    #........................
    def check_testMixedArithmetic(self):
        "Tests mixed arithmetics."
        na = N.array([1])
        ma = array([1])
        self.failUnless(isinstance(na + ma, MaskedArray))
        self.failUnless(isinstance(ma + na, MaskedArray))    
    #.........................
    def check_testUfuncs1 (self):
        "Test various functions such as sin, cos."
        (x, y, a10, m1, m2, xm, ym, z, zm, xf) = self.d
        assert_equal(N.cos(x), cos(xm))
        assert_equal(N.cosh(x), cosh(xm))
        assert_equal(N.sin(x), sin(xm))
        assert_equal(N.sinh(x), sinh(xm))
        assert_equal(N.tan(x), tan(xm))
        assert_equal(N.tanh(x), tanh(xm))
        assert_equal(N.sqrt(abs(x)), sqrt(xm))
        assert_equal(N.log(abs(x)), log(xm))
        assert_equal(N.log10(abs(x)), log10(xm))
        assert_equal(N.exp(x), exp(xm))
        assert_equal(N.arcsin(z), arcsin(zm))
        assert_equal(N.arccos(z), arccos(zm))
        assert_equal(N.arctan(z), arctan(zm))
        assert_equal(N.arctan2(x, y), arctan2(xm, ym))
        assert_equal(N.absolute(x), absolute(xm))
        assert_equal(N.equal(x,y), equal(xm, ym))
        assert_equal(N.not_equal(x,y), not_equal(xm, ym))
        assert_equal(N.less(x,y), less(xm, ym))
        assert_equal(N.greater(x,y), greater(xm, ym))
        assert_equal(N.less_equal(x,y), less_equal(xm, ym))
        assert_equal(N.greater_equal(x,y), greater_equal(xm, ym))
        assert_equal(N.conjugate(x), conjugate(xm))
        assert_equal(N.concatenate((x,y)), concatenate((xm,ym)))
        assert_equal(N.concatenate((x,y)), concatenate((x,y)))
        assert_equal(N.concatenate((x,y)), concatenate((xm,y)))
        assert_equal(N.concatenate((x,y,x)), concatenate((x,ym,x)))
    #........................
    def check_xtestCount (self):
        "Tests count"
        ott = array([0.,1.,2.,3.], mask=[1,0,0,0])
        assert( isinstance(count(ott), types.IntType))
        assert_equal(3, count(ott))
        assert_equal(1, count(1))
        assert_equal(0, array(1,mask=[1]))
        ott = ott.reshape((2,2))
        assert isMaskedArray(count(ott,0))
        assert isinstance(count(ott), types.IntType)
        assert_equal(3, count(ott))
        assert getmask(count(ott,0)) is nomask
        assert_equal([1,2],count(ott,0))
    #........................
    def check_testMinMax (self):
        "Tests minimum and maximum."
        (x, y, a10, m1, m2, xm, ym, z, zm, xf) = self.d
        xr = N.ravel(x) #max doesn't work if shaped
        xmr = ravel(xm)
        assert_equal(max(xr), maximum(xmr)) #true because of careful selection of data
        assert_equal(min(xr), minimum(xmr)) #true because of careful selection of data
        #
        assert_equal(minimum([1,2,3],[4,0,9]), [1,0,3])
        assert_equal(maximum([1,2,3],[4,0,9]), [4,2,9])
        x = arange(5)
        y = arange(5) - 2
        x[3] = masked
        y[0] = masked
        assert_equal(minimum(x,y), where(less(x,y), x, y))
        assert_equal(maximum(x,y), where(greater(x,y), x, y))
        assert minimum(x) == 0
        assert maximum(x) == 4
    #........................
    def check_testAddSumProd (self):
        "Tests add, sum, product."
        (x, y, a10, m1, m2, xm, ym, z, zm, xf) = self.d
        assert_equal(N.add.reduce(x), add.reduce(x))
        assert_equal(N.add.accumulate(x), add.accumulate(x))
        assert_equal(4, sum(array(4),axis=0))
        assert_equal(4, sum(array(4), axis=0))
        assert_equal(N.sum(x,axis=0), sum(x,axis=0))
        assert_equal(N.sum(filled(xm,0),axis=0), sum(xm,axis=0))
        assert_equal(N.sum(x,0), sum(x,0))
        assert_equal(N.product(x,axis=0), product(x,axis=0))
        assert_equal(N.product(x,0), product(x,0))
        assert_equal(N.product(filled(xm,1),axis=0), product(xm,axis=0))
        s = (3,4)
        x.shape = y.shape = xm.shape = ym.shape = s
        if len(s) > 1:
            assert_equal(N.concatenate((x,y),1), concatenate((xm,ym),1))
            assert_equal(N.add.reduce(x,1), add.reduce(x,1))
            assert_equal(N.sum(x,1), sum(x,1))
            assert_equal(N.product(x,1), product(x,1))
    #.........................
    def check_concat(self):
        "Tests concatenations."
        (x, y, a10, m1, m2, xm, ym, z, zm, xf) = self.d
        s = (3,4)
        x.shape = y.shape = xm.shape = ym.shape = s
        assert_equal(xm.mask, N.reshape(m1, s))
        assert_equal(ym.mask, N.reshape(m2, s))
        xmym = concatenate((xm,ym),1)
        assert_equal(N.concatenate((x,y),1), xmym)
        assert_equal(N.concatenate((xm.mask,ym.mask),1), xmym._mask)
    #........................
    def check_testCI(self):
        "Tests conversions and indexing"
        x1 = N.array([1,2,4,3])
        x2 = array(x1, mask=[1,0,0,0])
        x3 = array(x1, mask=[0,1,0,1])
        x4 = array(x1)
    # test conversion to strings
        junk, garbage = str(x2), repr(x2)
        assert_equal(N.sort(x1),sort(x2, fill_value=0))
    # tests of indexing
        assert type(x2[1]) is type(x1[1])
        assert x1[1] == x2[1]
        assert x2[0] is masked
        assert_equal(x1[2],x2[2])
        assert_equal(x1[2:5],x2[2:5])
        assert_equal(x1[:],x2[:])
        assert_equal(x1[1:], x3[1:])
        x1[2] = 9
        x2[2] = 9
        assert_equal(x1,x2)
        x1[1:3] = 99
        x2[1:3] = 99
        assert_equal(x1,x2)
        x2[1] = masked
        assert_equal(x1,x2)
        x2[1:3] = masked
        assert_equal(x1,x2)
        x2[:] = x1
        x2[1] = masked
        assert allequal(getmask(x2),array([0,1,0,0]))
        x3[:] = masked_array([1,2,3,4],[0,1,1,0])
        assert allequal(getmask(x3), array([0,1,1,0]))
        x4[:] = masked_array([1,2,3,4],[0,1,1,0])
        assert allequal(getmask(x4), array([0,1,1,0]))
        assert allequal(x4, array([1,2,3,4]))
        x1 = N.arange(5)*1.0
        x2 = masked_values(x1, 3.0)
        assert_equal(x1,x2)
        assert allequal(array([0,0,0,1,0],MaskType), x2.mask)
#FIXME: Well, eh, fill_value is now a property        assert_equal(3.0, x2.fill_value())
        assert_equal(3.0, x2.fill_value)
        x1 = array([1,'hello',2,3],object)
        x2 = N.array([1,'hello',2,3],object)
        s1 = x1[1]
        s2 = x2[1]
        assert_equal(type(s2), str)
        assert_equal(type(s1), str)
        assert_equal(s1, s2)
        assert x1[1:1].shape == (0,)
    #........................
    def check_testCopySize(self):
        "Tests of some subtle points of copying and sizing."
        n = [0,0,1,0,0]
        m = make_mask(n)
        m2 = make_mask(m)
        assert(m is m2)
        m3 = make_mask(m, copy=1)
        assert(m is not m3)

        x1 = N.arange(5)
        y1 = array(x1, mask=m)
        assert( y1._data is x1)
        assert( allequal(x1,y1.raw_data()))
        assert( y1.mask is m)

        y1a = array(y1)
        assert( y1a.raw_data() is y1.raw_data())
        assert( y1a.mask is y1.mask)

        y2 = array(x1, mask=m)
        assert( y2.raw_data() is x1)
        assert( y2.mask is m)
        assert( y2[2] is masked)
        y2[2] = 9
        assert( y2[2] is not masked)
        assert( y2.mask is not m)
        assert( allequal(y2.mask, 0))

        y3 = array(x1*1.0, mask=m)
        assert(filled(y3).dtype is (x1*1.0).dtype)

        x4 = arange(4)
        x4[2] = masked
        y4 = resize(x4, (8,))
        assert_equal(concatenate([x4,x4]), y4)
        assert_equal(getmask(y4),[0,0,1,0,0,0,1,0])
        y5 = repeat(x4, (2,2,2,2), axis=0)
        assert_equal(y5, [0,0,1,1,2,2,3,3])
        y6 = repeat(x4, 2, axis=0)
        assert_equal(y5, y6)
        y7 = x4.repeat((2,2,2,2), axis=0)
        assert_equal(y5,y7)
        y8 = x4.repeat(2,0)
        assert_equal(y5,y8)
        
        y9 = x4.copy()
        assert_equal(y9._data, x4._data)
        assert_equal(y9._mask, x4._mask)
        
    #........................
    def check_testOddFeatures_1(self):
        "Test of other odd features"
        x = arange(20)
        x = x.reshape(4,5)
        x.flat[5] = 12
        assert x[1,0] == 12
        z = x + 10j * x
        assert_equal(z.real, x)
        assert_equal(z.imag, 10*x)
        assert_equal((z*conjugate(z)).real, 101*x*x)
        z.imag[...] = 0.0

        x = arange(10)
        x[3] = masked
        assert str(x[3]) == str(masked)
        c = x >= 8
        assert count(where(c,masked,masked)) == 0
        assert shape(where(c,masked,masked)) == c.shape
        #
        z = where(c , x, masked)
        assert z.dtype is x.dtype
        assert z[3] is masked
        assert z[4] is masked
        assert z[7] is masked
        assert z[8] is not masked
        assert z[9] is not masked
        assert_equal(x,z)
        #
        z = where(c , masked, x)
        assert z.dtype is x.dtype
        assert z[3] is masked
        assert z[4] is not masked
        assert z[7] is not masked
        assert z[8] is masked
        assert z[9] is masked
        #
        z = masked_where(c, x)
        assert z.dtype is x.dtype
        assert z[3] is masked
        assert z[4] is not masked
        assert z[7] is not masked
        assert z[8] is masked
        assert z[9] is masked
        assert_equal(x,z)
        #
    #........................
    def check_testOddFeatures_2(self):
        "Tests some more features."
        x = array([1.,2.,3.,4.,5.])
        c = array([1,1,1,0,0])
        x[2] = masked
        z = where(c, x, -x)
        assert_equal(z, [1.,2.,0., -4., -5])
        c[0] = masked
        z = where(c, x, -x)
        assert_equal(z, [1.,2.,0., -4., -5])
        assert z[0] is masked
        assert z[1] is not masked
        assert z[2] is masked
        #
        x = arange(6)
        x[5] = masked
        y = arange(6)*10
        y[2] = masked
        c = array([1,1,1,0,0,0], mask=[1,0,0,0,0,0])
        cm = c.filled(1)
        z = where(c,x,y)
        zm = where(cm,x,y)
        assert_equal(z, zm)
        assert getmask(zm) is nomask
        assert_equal(zm, [0,1,2,30,40,50])
        z = where(c, masked, 1)
        assert_equal(z, [99,99,99,1,1,1])
        z = where(c, 1, masked)
        assert_equal(z, [99, 1, 1, 99, 99, 99])
    #........................
    def check_testOddFeatures_3(self):
        """Tests some generic features."""
        atest = ones((10,10,10), dtype=float_)
        btest = zeros(atest.shape, MaskType)
        ctest = masked_where(btest,atest)
        assert_equal(atest,ctest)
    #........................
    def check_maskingfunctions(self):       
        "Tests masking functions."
        x = array([1.,2.,3.,4.,5.])
        x[2] = masked
        assert_equal(masked_where(greater(x, 2), x), masked_greater(x,2))
        assert_equal(masked_where(greater_equal(x, 2), x), masked_greater_equal(x,2))
        assert_equal(masked_where(less(x, 2), x), masked_less(x,2))
        assert_equal(masked_where(less_equal(x, 2), x), masked_less_equal(x,2))
        assert_equal(masked_where(not_equal(x, 2), x), masked_not_equal(x,2))
        assert_equal(masked_where(equal(x, 2), x), masked_equal(x,2))
        assert_equal(masked_where(not_equal(x,2), x), masked_not_equal(x,2))
        assert_equal(masked_inside(range(5), 1, 3), [0, 199, 199, 199, 4])
        assert_equal(masked_outside(range(5), 1, 3),[199,1,2,3,199])
        assert_equal(masked_inside(array(range(5), mask=[1,0,0,0,0]), 1, 3).mask, [1,1,1,1,0])
        assert_equal(masked_outside(array(range(5), mask=[0,1,0,0,0]), 1, 3).mask, [1,1,0,0,1])
        assert_equal(masked_equal(array(range(5), mask=[1,0,0,0,0]), 2).mask, [1,0,1,0,0])
        assert_equal(masked_not_equal(array([2,2,1,2,1], mask=[1,0,0,0,0]), 2).mask, [1,0,1,0,1])
        assert_equal(masked_where([1,1,0,0,0], [1,2,3,4,5]), [99,99,3,4,5])
    #........................
    def check_testTakeTransposeInnerOuter(self):
        "Test of take, transpose, inner, outer products"
        x = arange(24)
        y = N.arange(24)
        x[5:6] = masked
        x = x.reshape(2,3,4)
        y = y.reshape(2,3,4)
        assert_equal(N.transpose(y,(2,0,1)), transpose(x,(2,0,1)))
        assert_equal(N.take(y, (2,0,1), 1), take(x, (2,0,1), 1))
        assert_equal(N.inner(filled(x,0),filled(y,0)),
                            inner(x, y))
        assert_equal(N.outer(filled(x,0),filled(y,0)),
                            outer(x, y))
        y = array(['abc', 1, 'def', 2, 3], object)
        y[2] = masked
        t = take(y,[0,3,4])
        assert t[0] == 'abc'
        assert t[1] == 2
        assert t[2] == 3
    #........................
    def check_testInplace(self):
        """Test of inplace operations and rich comparisons"""
        y = arange(10)

        x = arange(10)
        xm = arange(10)
        xm[2] = masked
        x += 1
        assert_equal(x, y+1)
        xm += 1
        assert_equal(xm, y+1)

        x = arange(10)
        xm = arange(10)
        xm[2] = masked
        x -= 1
        assert_equal(x, y-1)
        xm -= 1
        assert_equal(xm, y-1)

        x = arange(10)*1.0
        xm = arange(10)*1.0
        xm[2] = masked
        x *= 2.0
        assert_equal(x, y*2)
        xm *= 2.0
        assert_equal(xm, y*2)

        x = arange(10)*2
        xm = arange(10)*2
        xm[2] = masked
        x /= 2
        assert_equal(x, y)
        xm /= 2
        assert_equal(xm, y)

        x = arange(10)*1.0
        xm = arange(10)*1.0
        xm[2] = masked
        x /= 2.0
        assert_equal(x, y/2.0)
        xm /= arange(10)
        assert_equal(xm, ones((10,)))

        x = arange(10).astype(float_)
        xm = arange(10)
        xm[2] = masked
        id1 = id(x.raw_data())
        x += 1.
        assert id1 == id(x.raw_data())
        assert_equal(x, y+1.)
        
        x = arange(10, dtype=float_)
        xm = arange(10, dtype=float_)
        xm[2] = masked
        m = xm.mask
        a = arange(10, dtype=float_)
        a[-1] = masked
        x += a
        xm += a
        assert_equal(x,y+a)
        assert_equal(xm,y+a)
        assert_equal(xm.mask, mask_or(m,a.mask))
        
        x = arange(10, dtype=float_)
        xm = arange(10, dtype=float_)
        xm[2] = masked
        m = xm.mask
        a = arange(10, dtype=float_)
        a[-1] = masked
        x -= a
        xm -= a
        assert_equal(x,y-a)
        assert_equal(xm,y-a)
        assert_equal(xm.mask, mask_or(m,a.mask))        
        
        x = arange(10, dtype=float_)
        xm = arange(10, dtype=float_)
        xm[2] = masked
        m = xm.mask
        a = arange(10, dtype=float_)
        a[-1] = masked
        x *= a
        xm *= a
        assert_equal(x,y*a)
        assert_equal(xm,y*a)
        assert_equal(xm.mask, mask_or(m,a.mask))        
                
        x = arange(10, dtype=float_)
        xm = arange(10, dtype=float_)
        xm[2] = masked
        m = xm.mask
        a = arange(10, dtype=float_)
        a[-1] = masked
        x /= a
        xm /= a
        assert_equal(x,y/a)
        assert_equal(xm,y/a)
        assert_equal(xm.mask, mask_or(mask_or(m,a.mask), (a==0)))
    #........................
    def check_testPickle(self):
        "Test of pickling"
        import pickle
        x = arange(12)
        x[4:10:2] = masked
        x = x.reshape(4,3)
        s = pickle.dumps(x)
        y = pickle.loads(s)
        assert_equal(x,y)
    #.......................
    def check_testMasked(self):
        "Test of masked element"
        x = arange(6)
        x[1] = masked
        assert(str(masked) ==  '--')
        assert(x[1] is masked)
        assert_equal(filled(x[1], 0), 0)
        # don't know why these should raise an exception...
        #self.failUnlessRaises(Exception, lambda x,y: x+y, masked, masked)
        #self.failUnlessRaises(Exception, lambda x,y: x+y, masked, 2)
        #self.failUnlessRaises(Exception, lambda x,y: x+y, masked, xx)
        #self.failUnlessRaises(Exception, lambda x,y: x+y, xx, masked)
    #........................
    def check_testToPython(self):
        "Tests some communication issues with Python."
        assert_equal(1, int(array(1)))
        assert_equal(1.0, float(array(1)))
        assert_equal(1, int(array([[[1]]])))
        assert_equal(1.0, float(array([[1]])))
        self.failUnlessRaises(ValueError, float, array([1,1]))
        assert N.isnan(float(array([1],mask=[1])))
#TODO: Check how bool works...        
#TODO:        self.failUnless(bool(array([0,1])))
#TODO:        self.failUnless(bool(array([0,0],mask=[0,1])))
#TODO:        self.failIf(bool(array([0,0])))
#TODO:        self.failIf(bool(array([0,0],mask=[0,0])))
    #..........................
    def check_testScalarArithmetic(self):
        "Tests some scalar arithmetics on MaskedArrays."
        xm = array(0, mask=1)
        assert((1/array(0)).mask)
        assert((1 + xm).mask)
        assert((-xm).mask)
        assert((-xm).mask)
        assert(maximum(xm, xm).mask)
        assert(minimum(xm, xm).mask)
        assert(xm.filled().dtype is xm.data.dtype)
        x = array(0, mask=0)
        assert(x.filled() is x.data)
        assert_equal(str(xm), str(masked_print_option))
    #........................
    def check_testArrayMethods(self):
        "Tests some MaskedArray methods."
        a = array([1,3,2])
        b = array([1,3,2], mask=[1,0,1])
        assert_equal(a.any(), a.data.any())
        assert_equal(a.all(), a.data.all())
        assert_equal(a.argmax(), a.data.argmax())
        assert_equal(a.argmin(), a.data.argmin())
        assert_equal(a.choose(0,1,2,3,4), a.data.choose(0,1,2,3,4))
        assert_equal(a.compress([1,0,1]), a.data.compress([1,0,1]))
        assert_equal(a.conj(), a.data.conj())
        assert_equal(a.conjugate(), a.data.conjugate())
        #
        m = array([[1,2],[3,4]])
        assert_equal(m.diagonal(), m.data.diagonal())
        assert_equal(a.sum(), a.data.sum())
        assert_equal(a.take([1,2]), a.data.take([1,2]))
        assert_equal(m.transpose(), m.data.transpose())
    #........................
    def check_testArrayAttributes(self):
        "Tests some basic array attributes."
        a = array([1,3,2])
        b = array([1,3,2], mask=[1,0,1])
        assert_equal(a.ndim, 1)
        assert_equal(b.ndim, 1)
        assert_equal(a.size, 3)
        assert_equal(b.size, 3)
        assert_equal(a.shape, (3,))
        assert_equal(b.shape, (3,))
    #........................
    def check_testSingleElementSubscript(self):
        "Tests single element subscripts of Maskedarrays."
        a = array([1,3,2])
        b = array([1,3,2], mask=[1,0,1])
        assert_equal(a[0].shape, ())
        assert_equal(b[0].shape, ())
        assert_equal(b[1].shape, ())
    #........................
    def check_maskcreation(self):
        "Tests how masks are initialized at the creation of Maskedarrays."
        data = arange(24, dtype=float_)
        data[[3,6,15]] = masked
        dma_1 = MaskedArray(data)
        assert_equal(dma_1.mask, data.mask)
        dma_2 = MaskedArray(dma_1)
        assert_equal(dma_2.mask, dma_1.mask)
        dma_3 = MaskedArray(dma_1, mask=[1,0,0,0]*6)
        fail_if_equal(dma_3.mask, dma_1.mask)
        
    def check_backwards(self):
        "Tests backward compatibility with numpy.core.ma"
        import numpy.core.ma as nma
        x = nma.arange(5)
        x[2] = nma.masked
        X = masked_array(x, mask=x._mask)
        assert_equal(X._mask, x.mask)
        assert_equal(X._data, x._data)
        X = masked_array(x)
        assert_equal(X._data, x._data)
        assert_equal(X._mask, x.mask)
        assert_equal(getmask(x), [0,0,1,0,0])
        
#..............................................................................
        
class test_ufuncs(NumpyTestCase):
    "Test class for the application of ufuncs on MaskedArrays."
    def setUp(self):
        "Base data definition."
        self.d = (array([1.0, 0, -1, pi/2]*2, mask=[0,1]+[0]*6),
                  array([1.0, 0, -1, pi/2]*2, mask=[1,0]+[0]*6),)

    def check_testUfuncRegression(self):
        "Tests new ufuncs on MaskedArrays."
        for f in ['sqrt', 'log', 'log10', 'exp', 'conjugate',
                  'sin', 'cos', 'tan',
                  'arcsin', 'arccos', 'arctan',
                  'sinh', 'cosh', 'tanh',
                  'arcsinh',
                  'arccosh',
                  'arctanh',
                  'absolute', 'fabs', 'negative',
                  # 'nonzero', 'around',
                  'floor', 'ceil',
                  # 'sometrue', 'alltrue',
                  'logical_not',
                  'add', 'subtract', 'multiply',
                  'divide', 'true_divide', 'floor_divide',
                  'remainder', 'fmod', 'hypot', 'arctan2',
                  'equal', 'not_equal', 'less_equal', 'greater_equal',
                  'less', 'greater',
                  'logical_and', 'logical_or', 'logical_xor',
                  ]:
            #print f
            try:
                uf = getattr(umath, f)
            except AttributeError:
                uf = getattr(fromnumeric, f)
            mf = getattr(maskedarray, f)
            args = self.d[:uf.nin]
            ur = uf(*args)
            mr = mf(*args)
            assert_equal(ur.filled(0), mr.filled(0), f)
            assert_mask_equal(ur.mask, mr.mask)
    #........................
    def test_reduce(self):
        "Tests reduce on MaskedArrays."
        a = self.d[0]
        assert(not alltrue(a,axis=0))
        assert(sometrue(a,axis=0))
        assert_equal(sum(a[:3],axis=0), 0)
        assert_equal(product(a,axis=0), 0)
    #........................
    def test_minmax(self):
        "Tests extrema on MaskedArrays."
        a = arange(1,13).reshape(3,4)
        amask = masked_where(a < 5,a)
        assert_equal(amask.max(), a.max())
        assert_equal(amask.min(), 5)
        assert_equal(amask.max(0), a.max(0))
        assert_equal(amask.min(0), [5,6,7,8])
        assert(amask.max(1)[0].mask)
        assert(amask.min(1)[0].mask)

#..............................................................................

class test_array_methods(NumpyTestCase):
    "Test class for miscellaneous MaskedArrays methods."
    def setUp(self):
        "Base data definition."
        x = N.array([ 8.375,  7.545,  8.828,  8.5  ,  1.757,  5.928,  
                      8.43 ,  7.78 ,  9.865,  5.878,  8.979,  4.732,  
                      3.012,  6.022,  5.095,  3.116,  5.238,  3.957,  
                      6.04 ,  9.63 ,  7.712,  3.382,  4.489,  6.479,
                      7.189,  9.645,  5.395,  4.961,  9.894,  2.893,  
                      7.357,  9.828,  6.272,  3.758,  6.693,  0.993])
        X = x.reshape(6,6)
        XX = x.reshape(3,2,2,3)
    
        m = N.array([0, 1, 0, 1, 0, 0, 
                     1, 0, 1, 1, 0, 1, 
                     0, 0, 0, 1, 0, 1, 
                     0, 0, 0, 1, 1, 1, 
                     1, 0, 0, 1, 0, 0, 
                     0, 0, 1, 0, 1, 0])
        mx = array(data=x,mask=m)
        mX = array(data=X,mask=m.reshape(X.shape))
        mXX = array(data=XX,mask=m.reshape(XX.shape))
    
        m2 = N.array([1, 1, 0, 1, 0, 0, 
                      1, 1, 1, 1, 0, 1, 
                      0, 0, 1, 1, 0, 1, 
                      0, 0, 0, 1, 1, 1, 
                      1, 0, 0, 1, 1, 0, 
                      0, 0, 1, 0, 1, 1])
        m2x = array(data=x,mask=m2)
        m2X = array(data=X,mask=m2.reshape(X.shape))
        m2XX = array(data=XX,mask=m2.reshape(XX.shape))
        self.d =  (x,X,XX,m,mx,mX,mXX,m2x,m2X,m2XX)

    #------------------------------------------------------
    def test_trace(self):
        "Tests trace on MaskedArrays."
        (x,X,XX,m,mx,mX,mXX,m2x,m2X,m2XX) = self.d
        mXdiag = mX.diagonal()
        assert_equal(mX.trace(), mX.diagonal().compressed().sum())
        assert_almost_equal(mX.trace(),
                            X.trace() - sum(mXdiag.mask*X.diagonal(),axis=0))

    def test_clip(self):
        "Tests clip on MaskedArrays."
        (x,X,XX,m,mx,mX,mXX,m2x,m2X,m2XX) = self.d
        clipped = mx.clip(2,8)
        assert_equal(clipped.mask,mx.mask)
        assert_equal(clipped.data,x.clip(2,8))
        assert_equal(clipped.data,mx.data.clip(2,8))

    def test_ptp(self):
        "Tests ptp on MaskedArrays."
        (x,X,XX,m,mx,mX,mXX,m2x,m2X,m2XX) = self.d
        (n,m) = X.shape
        assert_equal(mx.ptp(),mx.compressed().ptp())
        rows = N.zeros(n,N.float_)
        cols = N.zeros(m,N.float_)
        for k in range(m):
            cols[k] = mX[:,k].compressed().ptp()
        for k in range(n):
            rows[k] = mX[k].compressed().ptp()
        assert_equal(mX.ptp(0),cols)
        assert_equal(mX.ptp(1),rows)

    def test_swapaxes(self):
        "Tests swapaxes on MaskedArrays."
        (x,X,XX,m,mx,mX,mXX,m2x,m2X,m2XX) = self.d
        mXswapped = mX.swapaxes(0,1)
        assert_equal(mXswapped[-1],mX[:,-1])
        mXXswapped = mXX.swapaxes(0,2)
        assert_equal(mXXswapped.shape,(2,2,3,3))

    def test_cumsumprod(self):
        "Tests cumsum & cumprod on MaskedArrays."
        (x,X,XX,m,mx,mX,mXX,m2x,m2X,m2XX) = self.d
        mXcp = mX.cumsum(0)
        assert_equal(mXcp.data,mX.filled(0).cumsum(0))
        mXcp = mX.cumsum(1)
        assert_equal(mXcp.data,mX.filled(0).cumsum(1))
        #
        mXcp = mX.cumprod(0)
        assert_equal(mXcp.data,mX.filled(1).cumprod(0))
        mXcp = mX.cumprod(1)
        assert_equal(mXcp.data,mX.filled(1).cumprod(1))
                                                      
    def test_varstd(self):
        "Tests var & std on MaskedArrays."
        (x,X,XX,m,mx,mX,mXX,m2x,m2X,m2XX) = self.d
        assert_almost_equal(mX.var(axis=None),mX.compressed().var())
        assert_almost_equal(mX.std(axis=None),mX.compressed().std())
        assert_equal(mXX.var(axis=3).shape,XX.var(axis=3).shape)
        assert_equal(mX.var().shape,X.var().shape)
        (mXvar0,mXvar1) = (mX.var(axis=0), mX.var(axis=1))
        for k in range(6):
            assert_almost_equal(mXvar1[k],mX[k].compressed().var())
            assert_almost_equal(mXvar0[k],mX[:,k].compressed().var())
            assert_almost_equal(N.sqrt(mXvar0[k]), mX[:,k].compressed().std())
    
    def test_argmin(self):
        "Tests argmin & argmax on MaskedArrays."
        (x,X,XX,m,mx,mX,mXX,m2x,m2X,m2XX) = self.d
        #
        assert_equal(mx.argmin(),35)
        assert_equal(mX.argmin(),35)
        assert_equal(m2x.argmin(),4)
        assert_equal(m2X.argmin(),4)
        assert_equal(mx.argmax(),28)
        assert_equal(mX.argmax(),28)
        assert_equal(m2x.argmax(),31) 
        assert_equal(m2X.argmax(),31)    
        #
        assert_equal(mX.argmin(0), [2,2,2,5,0,5])
        assert_equal(m2X.argmin(0), [2,2,4,5,0,4])
        assert_equal(mX.argmax(0), [0,5,0,5,4,0])
        assert_equal(m2X.argmax(0), [5,5,0,5,1,0])
        #
        assert_equal(mX.argmin(1), [4,1,0,0,5,5,])
        assert_equal(m2X.argmin(1), [4,4,0,0,5,3])
        assert_equal(mX.argmax(1), [2,4,1,1,4,1])
        assert_equal(m2X.argmax(1), [2,4,1,1,1,1])
        
    def check_put(self):
        "Tests put."
        d = arange(5)
        n = [0,0,0,1,1]
        m = make_mask(n)
        x = array(d, mask = m)
        assert( x[3] is masked)
        assert( x[4] is masked)
        x[[1,4]] = [10,40]
        assert( x.mask is not m)
        assert( x[3] is masked)
        assert( x[4] is not masked)
        assert_equal(x, [0,10,2,-1,40])        
        #
        x = masked_array(arange(10), mask=[1,0,0,0,0]*2)
        i = [0,2,4,6]
        x.put(i, [6,4,2,0])
        assert_equal(x, asarray([6,1,4,3,2,5,0,7,8,9,]))
        assert_equal(x.mask, [0,0,0,0,0,1,0,0,0,0])
        x.put(i, masked_array([0,2,4,6],[1,0,1,0]))
        assert_array_equal(x, [0,1,2,3,4,5,6,7,8,9,])
        assert_equal(x.mask, [1,0,0,0,1,1,0,0,0,0])
        #
        x = masked_array(arange(10), mask=[1,0,0,0,0]*2)
        put(x, i, [6,4,2,0])
        assert_equal(x, asarray([6,1,4,3,2,5,0,7,8,9,]))
        assert_equal(x.mask, [0,0,0,0,0,1,0,0,0,0])
        put(x, i, masked_array([0,2,4,6],[1,0,1,0]))
        assert_array_equal(x, [0,1,2,3,4,5,6,7,8,9,])
        assert_equal(x.mask, [1,0,0,0,1,1,0,0,0,0])  
    
    def check_take(self):
        "Tests take"
        x = masked_array([10,20,30,40],[0,1,0,1])
        assert_equal(x.take([0,0,3]), masked_array([10, 10, 40], [0,0,1]) )
        assert_equal(x.take([0,0,3]), x[[0,0,3]])
        assert_equal(x.take([[0,1],[0,1]]), 
                     masked_array([[10,20],[10,20]], [[0,1],[0,1]]) )
        #
        x = array([[10,20,30],[40,50,60]], mask=[[0,0,1],[1,0,0,]])
        assert_equal(x.take([0,2], axis=1), 
                     array([[10,30],[40,60]], mask=[[0,1],[1,0]]))
        assert_equal(take(x, [0,2], axis=1), 
                      array([[10,30],[40,60]], mask=[[0,1],[1,0]]))       
        #........................
    def check_anyall(self):
        """Checks the any/all methods/functions."""
        x = N.array([[ 0.13,  0.26,  0.90],
                     [ 0.28,  0.33,  0.63],
                     [ 0.31,  0.87,  0.70]])
        m = N.array([[ True, False, False],
                     [False, False, False],
                     [True,  True, False]], dtype=N.bool_)
        mx = masked_array(x, mask=m)
        xbig = N.array([[False, False,  True],
                        [False, False,  True],
                        [False,  True,  True]], dtype=N.bool_)
        mxbig = (mx > 0.5)
        mxsmall = (mx < 0.5)
        #
        assert (mxbig.all()==False)
        assert (mxbig.any()==True)
        assert_equal(mxbig.all(0),[False, False, True])
        assert_equal(mxbig.all(1), [False, False, True])
        assert_equal(mxbig.any(0),[False, False, True])
        assert_equal(mxbig.any(1), [True, True, True])
        #
        assert (mxsmall.all()==False)
        assert (mxsmall.any()==True)
        assert_equal(mxsmall.all(0), [True,   True, False])
        assert_equal(mxsmall.all(1), [False, False, False])
        assert_equal(mxsmall.any(0), [True,   True, False])
        assert_equal(mxsmall.any(1), [True,   True, False])
        #
        X = N.matrix(x)
        mX = masked_array(X, mask=m)
        mXbig = (mX > 0.5)
        mXsmall = (mX < 0.5)
        #
        assert (mXbig.all()==False)
        assert (mXbig.any()==True)
        assert_equal(mXbig.all(0), N.matrix([False, False, True]))
        assert_equal(mXbig.all(1), N.matrix([False, False, True]).T)
        assert_equal(mXbig.any(0), N.matrix([False, False, True]))
        assert_equal(mXbig.any(1), N.matrix([ True,  True, True]).T)
        #
        assert (mXsmall.all()==False)
        assert (mXsmall.any()==True)
        assert_equal(mXsmall.all(0), N.matrix([True,   True, False]))
        assert_equal(mXsmall.all(1), N.matrix([False, False, False]).T)
        assert_equal(mXsmall.any(0), N.matrix([True,   True, False]))
        assert_equal(mXsmall.any(1), N.matrix([True,   True, False]).T)
   
    def check_keepmask(self):
        "Tests the keep mask flag"        
        x = masked_array([1,2,3], mask=[1,0,0])
        mx = masked_array(x)
        assert_equal(mx.mask, x.mask)
        mx = masked_array(x, mask=[0,1,0], keep_mask=False)
        assert_equal(mx.mask, [0,1,0])
        mx = masked_array(x, mask=[0,1,0], keep_mask=True)
        assert_equal(mx.mask, [1,1,0])   
        #We default to true
        mx = masked_array(x, mask=[0,1,0])
        assert_equal(mx.mask, [1,1,0])
#..............................................................................

#..............................................................................

class test_subclassing(NumpyTestCase):
    """Test suite for masked subclasses of ndarray."""
    
    class SubArray(N.ndarray):
        """Defines a generic N.ndarray subclass, that stores some metadata
        in the  dictionary `info`."""
        def __new__(cls,arr,info={}):
            x = N.array(arr).view(cls)
            x.info = info
            return x
        def __array_finalize__(self, obj):
            if hasattr(obj,'info'):
                self.info = obj.info
            return
            
    def check_subclassing(self):
        "Tests whether the subclass is kept."
        x = N.arange(5)
        m = [0,0,1,0,0]
        xsub = test_subclassing.SubArray(x)
        xmsub = masked_array(xsub, mask=m)
        assert isinstance(xmsub, MaskedArray)
        assert_equal(xmsub._data, xsub)
        assert isinstance(xmsub._data, test_subclassing.SubArray)

###############################################################################
#------------------------------------------------------------------------------
if __name__ == "__main__":
    NumpyTest().run()