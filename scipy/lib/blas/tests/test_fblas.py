# Test interfaces to fortran blas.
#
# The tests are more of interface than they are of the underlying blas.
# Only very small matrices checked -- N=3 or so.
#
# !! Complex calculations really aren't checked that carefully.
# !! Only real valued complex numbers are used in tests.

from numpy import zeros, transpose, newaxis, shape, float32, float64, \
                  complex64, complex128, arange, array, common_type, conjugate
from numpy.testing import *
from scipy.lib.blas import fblas

#decimal accuracy to require between Python and LAPACK/BLAS calculations
accuracy = 5

# Since numpy.dot likely uses the same blas, use this routine
# to check.
def matrixmultiply(a, b):
    if len(b.shape) == 1:
        b_is_vector = True
        b = b[:,newaxis]
    else:
        b_is_vector = False
    assert a.shape[1] == b.shape[0]
    c = zeros((a.shape[0], b.shape[1]), common_type(a, b))
    for i in xrange(a.shape[0]):
        for j in xrange(b.shape[1]):
            s = 0
            for k in xrange(a.shape[1]):
                s += a[i,k] * b[k, j]
            c[i,j] = s
    if b_is_vector:
        c = c.reshape((a.shape[0],))
    return c

##################################################
### Test blas ?axpy

class BaseAxpy(object):
    # Mixin class to test dtypes
    def test_default_a(self):
        x = arange(3.,dtype=self.dtype)
        y = arange(3.,dtype=x.dtype)
        real_y = x*1.+y
        self.blas_func(x,y)
        assert_array_almost_equal(real_y,y)
    def test_simple(self):
        x = arange(3.,dtype=self.dtype)
        y = arange(3.,dtype=x.dtype)
        real_y = x*3.+y
        self.blas_func(x,y,a=3.)
        assert_array_almost_equal(real_y,y)
    def test_x_stride(self):
        x = arange(6.,dtype=self.dtype)
        y = zeros(3,x.dtype)
        y = arange(3.,dtype=x.dtype)
        real_y = x[::2]*3.+y
        self.blas_func(x,y,a=3.,n=3,incx=2)
        assert_array_almost_equal(real_y,y)
    def test_y_stride(self):
        x = arange(3.,dtype=self.dtype)
        y = zeros(6,x.dtype)
        real_y = x*3.+y[::2]
        self.blas_func(x,y,a=3.,n=3,incy=2)
        assert_array_almost_equal(real_y,y[::2])
    def test_x_and_y_stride(self):
        x = arange(12.,dtype=self.dtype)
        y = zeros(6,x.dtype)
        real_y = x[::4]*3.+y[::2]
        self.blas_func(x,y,a=3.,n=3,incx=4,incy=2)
        assert_array_almost_equal(real_y,y[::2])
    def test_x_bad_size(self):
        x = arange(12.,dtype=self.dtype)
        y = zeros(6,x.dtype)
        try:
            self.blas_func(x,y,n=4,incx=5)
        except: # what kind of error should be caught?
            return
        # should catch error and never get here
        assert(0)
    def test_y_bad_size(self):
        x = arange(12.,dtype=complex64)
        y = zeros(6,x.dtype)
        try:
            self.blas_func(x,y,n=3,incy=5)
        except: # what kind of error should be caught?
            return
        # should catch error and never get here
        assert(0)

try:
    class TestSaxpy(TestCase, BaseAxpy):
        blas_func = fblas.saxpy
        dtype = float32
except AttributeError:
    class TestSaxpy: pass
class TestDaxpy(TestCase, BaseAxpy):
    blas_func = fblas.daxpy
    dtype = float64
try:
    class TestCaxpy(TestCase, BaseAxpy):
        blas_func = fblas.caxpy
        dtype = complex64
except AttributeError:
    class TestCaxpy: pass
class TestZaxpy(TestCase, BaseAxpy):
    blas_func = fblas.zaxpy
    dtype = complex128


##################################################
### Test blas ?scal

class BaseScal(object):
    # Mixin class for testing particular dtypes
    def test_simple(self):
        x = arange(3.,dtype=self.dtype)
        real_x = x*3.
        self.blas_func(3.,x)
        assert_array_almost_equal(real_x,x)
    def test_x_stride(self):
        x = arange(6.,dtype=self.dtype)
        real_x = x.copy()
        real_x[::2] = x[::2]*array(3.,self.dtype)
        self.blas_func(3.,x,n=3,incx=2)
        assert_array_almost_equal(real_x,x)
    def test_x_bad_size(self):
        x = arange(12.,dtype=self.dtype)
        try:
            self.blas_func(2.,x,n=4,incx=5)
        except: # what kind of error should be caught?
            return
        # should catch error and never get here
        assert(0)
try:
    class TestSscal(TestCase, BaseScal):
        blas_func = fblas.sscal
        dtype = float32
except AttributeError:
    class TestSscal: pass
class TestDscal(TestCase, BaseScal):
    blas_func = fblas.dscal
    dtype = float64
try:
    class TestCscal(TestCase, BaseScal):
        blas_func = fblas.cscal
        dtype = complex64
except AttributeError:
    class TestCscal: pass
class TestZscal(TestCase, BaseScal):
    blas_func = fblas.zscal
    dtype = complex128




##################################################
### Test blas ?copy

class BaseCopy(object):
    # Mixin class for testing dtypes
    def test_simple(self):
        x = arange(3.,dtype=self.dtype)
        y = zeros(shape(x),x.dtype)
        self.blas_func(x,y)
        assert_array_almost_equal(x,y)
    def test_x_stride(self):
        x = arange(6.,dtype=self.dtype)
        y = zeros(3,x.dtype)
        self.blas_func(x,y,n=3,incx=2)
        assert_array_almost_equal(x[::2],y)
    def test_y_stride(self):
        x = arange(3.,dtype=self.dtype)
        y = zeros(6,x.dtype)
        self.blas_func(x,y,n=3,incy=2)
        assert_array_almost_equal(x,y[::2])
    def test_x_and_y_stride(self):
        x = arange(12.,dtype=self.dtype)
        y = zeros(6,x.dtype)
        self.blas_func(x,y,n=3,incx=4,incy=2)
        assert_array_almost_equal(x[::4],y[::2])
    def test_x_bad_size(self):
        x = arange(12.,dtype=self.dtype)
        y = zeros(6,x.dtype)
        try:
            self.blas_func(x,y,n=4,incx=5)
        except: # what kind of error should be caught?
            return
        # should catch error and never get here
        assert(0)
    def test_y_bad_size(self):
        x = arange(12.,dtype=complex64)
        y = zeros(6,x.dtype)
        try:
            self.blas_func(x,y,n=3,incy=5)
        except: # what kind of error should be caught?
            return
        # should catch error and never get here
        assert(0)
    #def test_y_bad_type(self):
    ##   Hmmm. Should this work?  What should be the output.
    #    x = arange(3.,dtype=self.dtype)
    #    y = zeros(shape(x))
    #    self.blas_func(x,y)
    #    assert_array_almost_equal(x,y)

try:
    class TestScopy(TestCase, BaseCopy):
        blas_func = fblas.scopy
        dtype = float32
except AttributeError:
    class TestScopy: pass
class TestDcopy(TestCase, BaseCopy):
    blas_func = fblas.dcopy
    dtype = float64
try:
    class TestCcopy(TestCase, BaseCopy):
        blas_func = fblas.ccopy
        dtype = complex64
except AttributeError:
    class TestCcopy: pass
class TestZcopy(TestCase, BaseCopy):
    blas_func = fblas.zcopy
    dtype = complex128


##################################################
### Test blas ?swap

class BaseSwap(object):
    # Mixin class to implement test objects
    def test_simple(self):
        x = arange(3.,dtype=self.dtype)
        y = zeros(shape(x),x.dtype)
        desired_x = y.copy()
        desired_y = x.copy()
        self.blas_func(x,y)
        assert_array_almost_equal(desired_x,x)
        assert_array_almost_equal(desired_y,y)
    def test_x_stride(self):
        x = arange(6.,dtype=self.dtype)
        y = zeros(3,x.dtype)
        desired_x = y.copy()
        desired_y = x.copy()[::2]
        self.blas_func(x,y,n=3,incx=2)
        assert_array_almost_equal(desired_x,x[::2])
        assert_array_almost_equal(desired_y,y)
    def test_y_stride(self):
        x = arange(3.,dtype=self.dtype)
        y = zeros(6,x.dtype)
        desired_x = y.copy()[::2]
        desired_y = x.copy()
        self.blas_func(x,y,n=3,incy=2)
        assert_array_almost_equal(desired_x,x)
        assert_array_almost_equal(desired_y,y[::2])

    def test_x_and_y_stride(self):
        x = arange(12.,dtype=self.dtype)
        y = zeros(6,x.dtype)
        desired_x = y.copy()[::2]
        desired_y = x.copy()[::4]
        self.blas_func(x,y,n=3,incx=4,incy=2)
        assert_array_almost_equal(desired_x,x[::4])
        assert_array_almost_equal(desired_y,y[::2])
    def test_x_bad_size(self):
        x = arange(12.,dtype=self.dtype)
        y = zeros(6,x.dtype)
        try:
            self.blas_func(x,y,n=4,incx=5)
        except: # what kind of error should be caught?
            return
        # should catch error and never get here
        assert(0)
    def test_y_bad_size(self):
        x = arange(12.,dtype=complex64)
        y = zeros(6,x.dtype)
        try:
            self.blas_func(x,y,n=3,incy=5)
        except: # what kind of error should be caught?
            return
        # should catch error and never get here
        assert(0)

try:
    class TestSswap(TestCase, BaseSwap):
        blas_func = fblas.sswap
        dtype = float32
except AttributeError:
    class TestSswap: pass
class TestDswap(TestCase, BaseSwap):
    blas_func = fblas.dswap
    dtype = float64
try:
    class TestCswap(TestCase, BaseSwap):
        blas_func = fblas.cswap
        dtype = complex64
except AttributeError:
    class TestCswap: pass
class TestZswap(TestCase, BaseSwap):
    blas_func = fblas.zswap
    dtype = complex128

##################################################
### Test blas ?gemv
### This will be a mess to test all cases.

class BaseGemv(object):
    # Mixin class to test dtypes
    def get_data(self,x_stride=1,y_stride=1):
        mult = array(1, dtype = self.dtype)
        if self.dtype in [complex64, complex128]:
            mult = array(1+1j, dtype = self.dtype)
        from numpy.random import normal
        alpha = array(1., dtype = self.dtype) * mult
        beta = array(1.,dtype = self.dtype) * mult
        a = normal(0.,1.,(3,3)).astype(self.dtype) * mult
        x = arange(shape(a)[0]*x_stride,dtype=self.dtype) * mult
        y = arange(shape(a)[1]*y_stride,dtype=self.dtype) * mult
        return alpha,beta,a,x,y
    def test_simple(self):
        alpha,beta,a,x,y = self.get_data()
        desired_y = alpha*matrixmultiply(a,x)+beta*y
        y = self.blas_func(alpha,a,x,beta,y)
        assert_array_almost_equal(desired_y,y)
    def test_default_beta_y(self):
        alpha,beta,a,x,y = self.get_data()
        desired_y = matrixmultiply(a,x)
        y = self.blas_func(1,a,x)
        assert_array_almost_equal(desired_y,y)
    def test_simple_transpose(self):
        alpha,beta,a,x,y = self.get_data()
        desired_y = alpha*matrixmultiply(transpose(a),x)+beta*y
        y = self.blas_func(alpha,a,x,beta,y,trans=1)
        assert_array_almost_equal(desired_y,y)
    def test_simple_transpose_conj(self):
        alpha,beta,a,x,y = self.get_data()
        desired_y = alpha*matrixmultiply(transpose(conjugate(a)),x)+beta*y
        y = self.blas_func(alpha,a,x,beta,y,trans=2)
        assert_array_almost_equal(desired_y,y)
    def test_x_stride(self):
        alpha,beta,a,x,y = self.get_data(x_stride=2)
        desired_y = alpha*matrixmultiply(a,x[::2])+beta*y
        y = self.blas_func(alpha,a,x,beta,y,incx=2)
        assert_array_almost_equal(desired_y,y)
    def test_x_stride_transpose(self):
        alpha,beta,a,x,y = self.get_data(x_stride=2)
        desired_y = alpha*matrixmultiply(transpose(a),x[::2])+beta*y
        y = self.blas_func(alpha,a,x,beta,y,trans=1,incx=2)
        assert_array_almost_equal(desired_y,y)
    def test_x_stride_assert(self):
        # What is the use of this test?
        alpha,beta,a,x,y = self.get_data(x_stride=2)
        try:
            y = self.blas_func(1,a,x,1,y,trans=0,incx=3)
            assert(0)
        except:
            pass
        try:
            y = self.blas_func(1,a,x,1,y,trans=1,incx=3)
            assert(0)
        except:
            pass
    def test_y_stride(self):
        alpha,beta,a,x,y = self.get_data(y_stride=2)
        desired_y = y.copy()
        desired_y[::2] = alpha*matrixmultiply(a,x)+beta*y[::2]
        y = self.blas_func(alpha,a,x,beta,y,incy=2)
        assert_array_almost_equal(desired_y,y)
    def test_y_stride_transpose(self):
        alpha,beta,a,x,y = self.get_data(y_stride=2)
        desired_y = y.copy()
        desired_y[::2] = alpha*matrixmultiply(transpose(a),x)+beta*y[::2]
        y = self.blas_func(alpha,a,x,beta,y,trans=1,incy=2)
        assert_array_almost_equal(desired_y,y)
    def test_y_stride_assert(self):
        # What is the use of this test?
        alpha,beta,a,x,y = self.get_data(y_stride=2)
        try:
            y = self.blas_func(1,a,x,1,y,trans=0,incy=3)
            assert(0)
        except:
            pass
        try:
            y = self.blas_func(1,a,x,1,y,trans=1,incy=3)
            assert(0)
        except:
            pass

try:
    class TestSgemv(TestCase, BaseGemv):
        blas_func = fblas.sgemv
        dtype = float32
except AttributeError:
    class TestSgemv: pass
class TestDgemv(TestCase, BaseGemv):
    blas_func = fblas.dgemv
    dtype = float64
try:
    class TestCgemv(TestCase, BaseGemv):
        blas_func = fblas.cgemv
        dtype = complex64
except AttributeError:
    class TestCgemv: pass
class TestZgemv(TestCase, BaseGemv):
    blas_func = fblas.zgemv
    dtype = complex128

"""
##################################################
### Test blas ?ger
### This will be a mess to test all cases.

class BaseGer(TestCase):
    def get_data(self,x_stride=1,y_stride=1):
        from numpy.random import normal
        alpha = array(1., dtype = self.dtype)
        a = normal(0.,1.,(3,3)).astype(self.dtype)
        x = arange(shape(a)[0]*x_stride,dtype=self.dtype)
        y = arange(shape(a)[1]*y_stride,dtype=self.dtype)
        return alpha,a,x,y
    def test_simple(self):
        alpha,a,x,y = self.get_data()
        # tranpose takes care of Fortran vs. C(and Python) memory layout
        desired_a = alpha*transpose(x[:,newaxis]*y) + a
        self.blas_func(x,y,a)
        assert_array_almost_equal(desired_a,a)
    def test_x_stride(self):
        alpha,a,x,y = self.get_data(x_stride=2)
        desired_a = alpha*transpose(x[::2,newaxis]*y) + a
        self.blas_func(x,y,a,incx=2)
        assert_array_almost_equal(desired_a,a)
    def test_x_stride_assert(self):
        alpha,a,x,y = self.get_data(x_stride=2)
        try:
            self.blas_func(x,y,a,incx=3)
            assert(0)
        except:
            pass
    def test_y_stride(self):
        alpha,a,x,y = self.get_data(y_stride=2)
        desired_a = alpha*transpose(x[:,newaxis]*y[::2]) + a
        self.blas_func(x,y,a,incy=2)
        assert_array_almost_equal(desired_a,a)

    def test_y_stride_assert(self):
        alpha,a,x,y = self.get_data(y_stride=2)
        try:
            self.blas_func(a,x,y,incy=3)
            assert(0)
        except:
            pass

class TestSger(BaseGer):
    blas_func = fblas.sger
    dtype = float32
class TestDger(BaseGer):
    blas_func = fblas.dger
    dtype = float64
"""
##################################################
### Test blas ?gerc
### This will be a mess to test all cases.

"""
class BaseGerComplex(BaseGer):
    def get_data(self,x_stride=1,y_stride=1):
        from numpy.random import normal
        alpha = array(1+1j, dtype = self.dtype)
        a = normal(0.,1.,(3,3)).astype(self.dtype)
        a = a + normal(0.,1.,(3,3)) * array(1j, dtype = self.dtype)
        x = normal(0.,1.,shape(a)[0]*x_stride).astype(self.dtype)
        x = x + x * array(1j, dtype = self.dtype)
        y = normal(0.,1.,shape(a)[1]*y_stride).astype(self.dtype)
        y = y + y * array(1j, dtype = self.dtype)
        return alpha,a,x,y
    def test_simple(self):
        alpha,a,x,y = self.get_data()
        # tranpose takes care of Fortran vs. C(and Python) memory layout
        a = a * array(0.,dtype = self.dtype)
        #desired_a = alpha*transpose(x[:,newaxis]*self.transform(y)) + a
        desired_a = alpha*transpose(x[:,newaxis]*y) + a
        #self.blas_func(x,y,a,alpha = alpha)
        fblas.cgeru(x,y,a,alpha = alpha)
        assert_array_almost_equal(desired_a,a)

    #def test_x_stride(self):
    #    alpha,a,x,y = self.get_data(x_stride=2)
    #    desired_a = alpha*transpose(x[::2,newaxis]*self.transform(y)) + a
    #    self.blas_func(x,y,a,incx=2)
    #    assert_array_almost_equal(desired_a,a)
    #def test_y_stride(self):
    #    alpha,a,x,y = self.get_data(y_stride=2)
    #    desired_a = alpha*transpose(x[:,newaxis]*self.transform(y[::2])) + a
    #    self.blas_func(x,y,a,incy=2)
    #    assert_array_almost_equal(desired_a,a)

class TestCgeru(BaseGerComplex):
    blas_func = fblas.cgeru
    dtype = complex64
    def transform(self,x):
        return x
class TestZgeru(BaseGerComplex):
    blas_func = fblas.zgeru
    dtype = complex128
    def transform(self,x):
        return x

class TestCgerc(BaseGerComplex):
    blas_func = fblas.cgerc
    dtype = complex64
    def transform(self,x):
        return conjugate(x)

class TestZgerc(BaseGerComplex):
    blas_func = fblas.zgerc
    dtype = complex128
    def transform(self,x):
        return conjugate(x)
"""

if __name__ == "__main__":
    run_module_suite()
