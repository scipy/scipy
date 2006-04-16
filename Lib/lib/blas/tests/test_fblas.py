# Test interfaces to fortran blas.
#
# The tests are more of interface than they are of the underlying blas.
# Only very small matrices checked -- N=3 or so.
#
# These test really need to be checked on 64 bit architectures.
# What does complex32 become on such machines? complex64 I'll bet.
# If so, I think we are OK.
# Check when we have a machine to check on.
#
# !! Complex calculations really aren't checked that carefully.
# !! Only real valued complex numbers are used in tests.
#
# !! Uses matrixmultiply to check against blas.  If matrixmultiply is
# !! ever !replaced! by a blas call, we'll need to fill in a simple
# !! matrix multiply here to ensure integrity of tests.

from numpy import *

import sys
from numpy.testing import *
set_package_path()
from blas import fblas
restore_path()

#decimal accuracy to require between Python and LAPACK/BLAS calculations
accuracy = 5

##################################################
### Test blas ?axpy

class base_axpy(ScipyTestCase):
    def check_default_a(self):
        x = arange(3.,dtype=self.dtype)
        y = arange(3.,dtype=x.dtype)
        real_y = x*1.+y
        self.blas_func(x,y)
        assert_array_almost_equal(real_y,y)
    def check_simple(self):
        x = arange(3.,dtype=self.dtype)
        y = arange(3.,dtype=x.dtype)
        real_y = x*3.+y
        self.blas_func(x,y,a=3.)
        assert_array_almost_equal(real_y,y)
    def check_x_stride(self):
        x = arange(6.,dtype=self.dtype)
        y = zeros(3,x.dtype)
        y = arange(3.,dtype=x.dtype)
        real_y = x[::2]*3.+y
        self.blas_func(x,y,a=3.,n=3,incx=2)
        assert_array_almost_equal(real_y,y)
    def check_y_stride(self):
        x = arange(3.,dtype=self.dtype)
        y = zeros(6,x.dtype)
        real_y = x*3.+y[::2]
        self.blas_func(x,y,a=3.,n=3,incy=2)
        assert_array_almost_equal(real_y,y[::2])
    def check_x_and_y_stride(self):
        x = arange(12.,dtype=self.dtype)
        y = zeros(6,x.dtype)
        real_y = x[::4]*3.+y[::2]
        self.blas_func(x,y,a=3.,n=3,incx=4,incy=2)
        assert_array_almost_equal(real_y,y[::2])
    def check_x_bad_size(self):
        x = arange(12.,dtype=self.dtype)
        y = zeros(6,x.dtype)
        try:
            self.blas_func(x,y,n=4,incx=5)
        except: # what kind of error should be caught?
            return
        # should catch error and never get here
        assert(0)
    def check_y_bad_size(self):
        x = arange(12.,dtype=complex64)
        y = zeros(6,x.dtype)
        try:
            self.blas_func(x,y,n=3,incy=5)
        except: # what kind of error should be caught?
            return
        # should catch error and never get here
        assert(0)

try:
    class test_saxpy(base_axpy):
        blas_func = fblas.saxpy
        dtype = float32
except AttributeError:
    class test_saxpy: pass
class test_daxpy(base_axpy):
    blas_func = fblas.daxpy
    dtype = float64
try:
    class test_caxpy(base_axpy):
        blas_func = fblas.caxpy
        dtype = complex64
except AttributeError:
    class test_caxpy: pass
class test_zaxpy(base_axpy):
    blas_func = fblas.zaxpy
    dtype = complex128


##################################################
### Test blas ?scal

class base_scal(ScipyTestCase):
    def check_simple(self):
        x = arange(3.,dtype=self.dtype)
        real_x = x*3.
        self.blas_func(3.,x)
        assert_array_almost_equal(real_x,x)
    def check_x_stride(self):
        x = arange(6.,dtype=self.dtype)
        real_x = x.copy()
        real_x[::2] = x[::2]*array(3.,self.dtype)
        self.blas_func(3.,x,n=3,incx=2)
        assert_array_almost_equal(real_x,x)
    def check_x_bad_size(self):
        x = arange(12.,dtype=self.dtype)
        try:
            self.blas_func(2.,x,n=4,incx=5)
        except: # what kind of error should be caught?
            return
        # should catch error and never get here
        assert(0)
try:
    class test_sscal(base_scal):
        blas_func = fblas.sscal
        dtype = float32
except AttributeError:
    class test_sscal: pass
class test_dscal(base_scal):
    blas_func = fblas.dscal
    dtype = float64
try:
    class test_cscal(base_scal):
        blas_func = fblas.cscal
        dtype = complex64
except AttributeError:
    class test_cscal: pass
class test_zscal(base_scal):
    blas_func = fblas.zscal
    dtype = complex128




##################################################
### Test blas ?copy

class base_copy(ScipyTestCase):
    def check_simple(self):
        x = arange(3.,dtype=self.dtype)
        y = zeros(shape(x),x.dtype)
        self.blas_func(x,y)
        assert_array_almost_equal(x,y)
    def check_x_stride(self):
        x = arange(6.,dtype=self.dtype)
        y = zeros(3,x.dtype)
        self.blas_func(x,y,n=3,incx=2)
        assert_array_almost_equal(x[::2],y)
    def check_y_stride(self):
        x = arange(3.,dtype=self.dtype)
        y = zeros(6,x.dtype)
        self.blas_func(x,y,n=3,incy=2)
        assert_array_almost_equal(x,y[::2])
    def check_x_and_y_stride(self):
        x = arange(12.,dtype=self.dtype)
        y = zeros(6,x.dtype)
        self.blas_func(x,y,n=3,incx=4,incy=2)
        assert_array_almost_equal(x[::4],y[::2])
    def check_x_bad_size(self):
        x = arange(12.,dtype=self.dtype)
        y = zeros(6,x.dtype)
        try:
            self.blas_func(x,y,n=4,incx=5)
        except: # what kind of error should be caught?
            return
        # should catch error and never get here
        assert(0)
    def check_y_bad_size(self):
        x = arange(12.,dtype=complex64)
        y = zeros(6,x.dtype)
        try:
            self.blas_func(x,y,n=3,incy=5)
        except: # what kind of error should be caught?
            return
        # should catch error and never get here
        assert(0)
    #def check_y_bad_type(self):
    ##   Hmmm. Should this work?  What should be the output.
    #    x = arange(3.,dtype=self.dtype)
    #    y = zeros(shape(x))
    #    self.blas_func(x,y)
    #    assert_array_almost_equal(x,y)

try:
    class test_scopy(base_copy):
        blas_func = fblas.scopy
        dtype = float32
except AttributeError:
    class test_scopy: pass
class test_dcopy(base_copy):
    blas_func = fblas.dcopy
    dtype = float64
try:
    class test_ccopy(base_copy):
        blas_func = fblas.ccopy
        dtype = complex64
except AttributeError:
    class test_ccopy: pass
class test_zcopy(base_copy):
    blas_func = fblas.zcopy
    dtype = complex128


##################################################
### Test blas ?swap

class base_swap(ScipyTestCase):
    def check_simple(self):
        x = arange(3.,dtype=self.dtype)
        y = zeros(shape(x),x.dtype)
        desired_x = y.copy()
        desired_y = x.copy()
        self.blas_func(x,y)
        assert_array_almost_equal(desired_x,x)
        assert_array_almost_equal(desired_y,y)
    def check_x_stride(self):
        x = arange(6.,dtype=self.dtype)
        y = zeros(3,x.dtype)
        desired_x = y.copy()
        desired_y = x.copy()[::2]
        self.blas_func(x,y,n=3,incx=2)
        assert_array_almost_equal(desired_x,x[::2])
        assert_array_almost_equal(desired_y,y)
    def check_y_stride(self):
        x = arange(3.,dtype=self.dtype)
        y = zeros(6,x.dtype)
        desired_x = y.copy()[::2]
        desired_y = x.copy()
        self.blas_func(x,y,n=3,incy=2)
        assert_array_almost_equal(desired_x,x)
        assert_array_almost_equal(desired_y,y[::2])

    def check_x_and_y_stride(self):
        x = arange(12.,dtype=self.dtype)
        y = zeros(6,x.dtype)
        desired_x = y.copy()[::2]
        desired_y = x.copy()[::4]
        self.blas_func(x,y,n=3,incx=4,incy=2)
        assert_array_almost_equal(desired_x,x[::4])
        assert_array_almost_equal(desired_y,y[::2])
    def check_x_bad_size(self):
        x = arange(12.,dtype=self.dtype)
        y = zeros(6,x.dtype)
        try:
            self.blas_func(x,y,n=4,incx=5)
        except: # what kind of error should be caught?
            return
        # should catch error and never get here
        assert(0)
    def check_y_bad_size(self):
        x = arange(12.,dtype=complex64)
        y = zeros(6,x.dtype)
        try:
            self.blas_func(x,y,n=3,incy=5)
        except: # what kind of error should be caught?
            return
        # should catch error and never get here
        assert(0)

try:
    class test_sswap(base_swap):
        blas_func = fblas.sswap
        dtype = float32
except AttributeError:
    class test_sswap: pass
class test_dswap(base_swap):
    blas_func = fblas.dswap
    dtype = float64
try:
    class test_cswap(base_swap):
        blas_func = fblas.cswap
        dtype = complex64
except AttributeError:
    class test_cswap: pass
class test_zswap(base_swap):
    blas_func = fblas.zswap
    dtype = complex128

##################################################
### Test blas ?gemv
### This will be a mess to test all cases.

class base_gemv(ScipyTestCase):
    def get_data(self,x_stride=1,y_stride=1):
        mult = array(1, dtype = self.dtype)
        if self.dtype in ['F', 'D']:
            mult = array(1+1j, dtype = self.dtype)
        from numpy.random import normal
        alpha = array(1., dtype = self.dtype) * mult
        beta = array(1.,dtype = self.dtype) * mult
        a = normal(0.,1.,(3,3)).astype(self.dtype) * mult
        x = arange(shape(a)[0]*x_stride,dtype=self.dtype) * mult
        y = arange(shape(a)[1]*y_stride,dtype=self.dtype) * mult
        return alpha,beta,a,x,y
    def check_simple(self):
        alpha,beta,a,x,y = self.get_data()
        desired_y = alpha*matrixmultiply(a,x)+beta*y
        y = self.blas_func(alpha,a,x,beta,y)
        assert_array_almost_equal(desired_y,y)
    def check_default_beta_y(self):
        alpha,beta,a,x,y = self.get_data()
        desired_y = matrixmultiply(a,x)
        y = self.blas_func(1,a,x)
        assert_array_almost_equal(desired_y,y)
    def check_simple_transpose(self):
        alpha,beta,a,x,y = self.get_data()
        desired_y = alpha*matrixmultiply(transpose(a),x)+beta*y
        y = self.blas_func(alpha,a,x,beta,y,trans=1)
        assert_array_almost_equal(desired_y,y)
    def check_simple_transpose_conj(self):
        alpha,beta,a,x,y = self.get_data()
        desired_y = alpha*matrixmultiply(transpose(conjugate(a)),x)+beta*y
        y = self.blas_func(alpha,a,x,beta,y,trans=2)
        assert_array_almost_equal(desired_y,y)
    def check_x_stride(self):
        alpha,beta,a,x,y = self.get_data(x_stride=2)
        desired_y = alpha*matrixmultiply(a,x[::2])+beta*y
        y = self.blas_func(alpha,a,x,beta,y,incx=2)
        assert_array_almost_equal(desired_y,y)
    def check_x_stride_transpose(self):
        alpha,beta,a,x,y = self.get_data(x_stride=2)
        desired_y = alpha*matrixmultiply(transpose(a),x[::2])+beta*y
        y = self.blas_func(alpha,a,x,beta,y,trans=1,incx=2)
        assert_array_almost_equal(desired_y,y)
    def check_x_stride_assert(self):
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
    def check_y_stride(self):
        alpha,beta,a,x,y = self.get_data(y_stride=2)
        desired_y = y.copy()
        desired_y[::2] = alpha*matrixmultiply(a,x)+beta*y[::2]
        y = self.blas_func(alpha,a,x,beta,y,incy=2)
        assert_array_almost_equal(desired_y,y)
    def check_y_stride_transpose(self):
        alpha,beta,a,x,y = self.get_data(y_stride=2)
        desired_y = y.copy()
        desired_y[::2] = alpha*matrixmultiply(transpose(a),x)+beta*y[::2]
        y = self.blas_func(alpha,a,x,beta,y,trans=1,incy=2)
        assert_array_almost_equal(desired_y,y)
    def check_y_stride_assert(self):
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
    class test_sgemv(base_gemv):
        blas_func = fblas.sgemv
        dtype = float32
except AttributeError:
    class test_sgemv: pass
class test_dgemv(base_gemv):
    blas_func = fblas.dgemv
    dtype = float64
try:
    class test_cgemv(base_gemv):
        blas_func = fblas.cgemv
        dtype = complex64
except AttributeError:
    class test_cgemv: pass
class test_zgemv(base_gemv):
    blas_func = fblas.zgemv
    dtype = complex128

"""
##################################################
### Test blas ?ger
### This will be a mess to test all cases.

class base_ger(ScipyTestCase):
    def get_data(self,x_stride=1,y_stride=1):
        from numpy.random import normal
        alpha = array(1., dtype = self.dtype)
        a = normal(0.,1.,(3,3)).astype(self.dtype)
        x = arange(shape(a)[0]*x_stride,dtype=self.dtype)
        y = arange(shape(a)[1]*y_stride,dtype=self.dtype)
        return alpha,a,x,y
    def check_simple(self):
        alpha,a,x,y = self.get_data()
        # tranpose takes care of Fortran vs. C(and Python) memory layout
        desired_a = alpha*transpose(x[:,newaxis]*y) + a
        self.blas_func(x,y,a)
        assert_array_almost_equal(desired_a,a)
    def check_x_stride(self):
        alpha,a,x,y = self.get_data(x_stride=2)
        desired_a = alpha*transpose(x[::2,newaxis]*y) + a
        self.blas_func(x,y,a,incx=2)
        assert_array_almost_equal(desired_a,a)
    def check_x_stride_assert(self):
        alpha,a,x,y = self.get_data(x_stride=2)
        try:
            self.blas_func(x,y,a,incx=3)
            assert(0)
        except:
            pass
    def check_y_stride(self):
        alpha,a,x,y = self.get_data(y_stride=2)
        desired_a = alpha*transpose(x[:,newaxis]*y[::2]) + a
        self.blas_func(x,y,a,incy=2)
        assert_array_almost_equal(desired_a,a)

    def check_y_stride_assert(self):
        alpha,a,x,y = self.get_data(y_stride=2)
        try:
            self.blas_func(a,x,y,incy=3)
            assert(0)
        except:
            pass

class test_sger(base_ger):
    blas_func = fblas.sger
    dtype = float32
class test_dger(base_ger):
    blas_func = fblas.dger
    dtype = float64
"""
##################################################
### Test blas ?gerc
### This will be a mess to test all cases.

"""
class base_ger_complex(base_ger):
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
    def check_simple(self):
        alpha,a,x,y = self.get_data()
        # tranpose takes care of Fortran vs. C(and Python) memory layout
        a = a * array(0.,dtype = self.dtype)
        #desired_a = alpha*transpose(x[:,newaxis]*self.transform(y)) + a
        desired_a = alpha*transpose(x[:,newaxis]*y) + a
        #self.blas_func(x,y,a,alpha = alpha)
        fblas.cgeru(x,y,a,alpha = alpha)
        print x, y
        print desired_a.dtype,desired_a
        print
        print a.dtype,a
        assert_array_almost_equal(desired_a,a)

    #def check_x_stride(self):
    #    alpha,a,x,y = self.get_data(x_stride=2)
    #    desired_a = alpha*transpose(x[::2,newaxis]*self.transform(y)) + a
    #    self.blas_func(x,y,a,incx=2)
    #    assert_array_almost_equal(desired_a,a)
    #def check_y_stride(self):
    #    alpha,a,x,y = self.get_data(y_stride=2)
    #    desired_a = alpha*transpose(x[:,newaxis]*self.transform(y[::2])) + a
    #    self.blas_func(x,y,a,incy=2)
    #    assert_array_almost_equal(desired_a,a)

class test_cgeru(base_ger_complex):
    blas_func = fblas.cgeru
    dtype = complex64
    def transform(self,x):
        return x
class test_zgeru(base_ger_complex):
    blas_func = fblas.zgeru
    dtype = complex128
    def transform(self,x):
        return x

class test_cgerc(base_ger_complex):
    blas_func = fblas.cgerc
    dtype = complex64
    def transform(self,x):
        return conjugate(x)

class test_zgerc(base_ger_complex):
    blas_func = fblas.zgerc
    dtype = complex128
    def transform(self,x):
        return conjugate(x)
"""

if __name__ == "__main__":
    ScipyTest().run()
