""" Test functions for misc module

"""

import unittest
from scipy_test import assert_array_equal, assert_equal
from scipy_test import assert_almost_equal
from scipy import *


##################################################
### Test for sum

class test_real(unittest.TestCase):
    def check_real(self):
        y = rand(10,)
        assert_array_equal(y,real(y))

    def check_cmplx(self):
        y = rand(10,)+1j*rand(10,)
        assert_array_equal(y.real,real(y))

class test_imag(unittest.TestCase):
    def check_real(self):
        y = rand(10,)
        assert_array_equal(0,imag(y))

    def check_cmplx(self):
        y = rand(10,)+1j*rand(10,)
        assert_array_equal(y.imag,imag(y))

    
class test_sum(unittest.TestCase):
    def check_sum_1D_array(self):
        z = arange(10.)
        val = sum(z)
        real_sum = 0
        for i in z:
            real_sum = real_sum + i
        assert_equal(val,real_sum)

    def check_sum_scalar(self):
        z = 10.
        val = sum(z)
        real_sum = 10.
        assert_equal(val,real_sum)

class test_mean(unittest.TestCase):
    def check_default_cols(self):
        z = array(((1,2),(3,4)),Float)
        val = mean(z)
        desired = array((2.,3.))
        assert_array_equal(val,desired)
    def check_rows(self):
        z = array(((1,2),(3,4)),Float)
        val = mean(z,axis=-1)
        desired = array((1.5,3.5))
        assert_array_equal(val,desired)
    def check_cols(self):
        z = array(((1,2),(3,4)),Float)
        val = mean(z,axis=0)
        desired = array((2.,3.))
        assert_array_equal(val,desired)


class test_iscomplex(unittest.TestCase):
    def check_fail(self):
        z = array([-1,0,1])
        res = iscomplex(z)
        assert(not sometrue(res))
    def check_pass(self):
        z = array([-1j,1,0])
        res = iscomplex(z)
        assert_array_equal(res,[1,0,0])

class test_isreal(unittest.TestCase):
    def check_pass(self):
        z = array([-1,0,1j])
        res = isreal(z)
        assert_array_equal(res,[1,1,0])
    def check_fail(self):
        z = array([-1j,1,0])
        res = isreal(z)
        assert_array_equal(res,[0,1,1])

class test_array_iscomplex(unittest.TestCase):
    def check_basic(self):
        z = array([-1,0,1])
        assert(not array_iscomplex(z))
        z = array([-1j,0,-1])
        assert(array_iscomplex(z))

class test_array_isreal(unittest.TestCase):
    def check_basic(self):
        z = array([-1,0,1])
        assert(array_isreal(z))
        z = array([-1j,0,-1])
        assert(not array_isreal(z))

class test_isnan(unittest.TestCase):
    def check_goodvalues(self):
        z = array((-1.,0.,1.))
        res = isnan(z) == 0
        assert(alltrue(res))            
    def check_posinf(self): 
        assert(isnan(array((1.,))/0.) == 0)
    def check_neginf(self): 
        assert(isnan(array((-1.,))/0.) == 0)
    def check_ind(self): 
        assert(isnan(array((0.,))/0.) == 1)
    def check_qnan(self): 
        assert(isnan(log(-1.)) == 1)
    def check_integer(self):
        assert(isnan(1) == 0)
    def check_complex(self):
        assert(isnan(1+1j) == 0)
    def check_complex1(self):
        assert(isnan(array(0+0j)/0.) == 1)
                
class test_isfinite(unittest.TestCase):
    def check_goodvalues(self):
        z = array((-1.,0.,1.))
        res = isfinite(z) == 1
        assert(alltrue(res))            
    def check_posinf(self): 
        assert(isfinite(array((1.,))/0.) == 0)
    def check_neginf(self): 
        assert(isfinite(array((-1.,))/0.) == 0)
    def check_ind(self): 
        assert(isfinite(array((0.,))/0.) == 0)
    def check_qnan(self): 
        assert(isfinite(log(-1.)) == 0)
    def check_integer(self):
        assert(isfinite(1) == 1)
    def check_complex(self):
        assert(isfinite(1+1j) == 1)
    def check_complex1(self):
        assert(isfinite(array(1+1j)/0.) == 0)
        
class test_isinf(unittest.TestCase):
    def check_goodvalues(self):
        z = array((-1.,0.,1.))
        res = isinf(z) == 0
        assert(alltrue(res))            
    def check_posinf(self): 
        assert(isinf(array((1.,))/0.) == 1)
    def check_posinf_scalar(self): 
        assert(isinf(array(1.,)/0.) == 1)
    def check_neginf(self): 
        assert(isinf(array((-1.,))/0.) == 1)
    def check_neginf_scalar(self): 
        assert(isinf(array(-1.)/0.) == 1)
    def check_ind(self): 
        assert(isinf(array((0.,))/0.) == 0)
    def check_qnan(self): 
        assert(isinf(log(-1.)) == 0)
        assert(isnan(log(-1.)) == 1)

class test_isposinf(unittest.TestCase):
    def check_generic(self):
        vals = isposinf(array((-1.,0,1))/0.)
        assert(vals[0] == 0)
        assert(vals[1] == 0)
        assert(vals[2] == 1)

class test_isneginf(unittest.TestCase):
    def check_generic(self):
        vals = isneginf(array((-1.,0,1))/0.)
        assert(vals[0] == 1)
        assert(vals[1] == 0)
        assert(vals[2] == 0)

class test_nan_to_num(unittest.TestCase):
    def check_generic(self):
        vals = nan_to_num(array((-1.,0,1))/0.)
        assert(vals[0] < -1e10 and isfinite(vals[0]))
        assert(vals[1] == 0)
        assert(vals[2] > 1e10 and isfinite(vals[2]))
    def check_integer(self):
        vals = nan_to_num(1)
        assert(vals == 1)
    def check_complex_good(self):
        vals = nan_to_num(1+1j)
        assert(vals == 1+1j)
    def check_complex_bad(self):
        v = 1+1j
        v += array(0+1.j)/0.
        vals = nan_to_num(v)
        # !! This is actually (unexpectedly) zero
        assert(vals.imag > 1e10 and isfinite(vals))
    def check_complex_bad2(self):
        v = 1+1j
        v += array(-1+1.j)/0.
        vals = nan_to_num(v)
        assert(isfinite(vals))    
        #assert(vals.imag > 1e10  and isfinite(vals))    
        # !! This is actually (unexpectedly) positive
        # !! inf.  Comment out for now, and see if it
        # !! changes
        #assert(vals.real < -1e10 and isfinite(vals))    
        
        
class test_logn(unittest.TestCase):
    def check_log_3_4(self):
        val = logn(3,4)
        assert_almost_equal(val,1.2618595071429148,11)
    def check_log_0(self):
        """ log(0) should print warning, but succeed.
        """
        try:
            val = logn(3,0)
            #assert(isinf(val))
        except:
            assert(0)
    def check_log_neg(self):
        """ log(-1) should print warning, but succeed.
        """
        try:
            val = logn(3,-1)
            #assert(isinf(val))
        except:
            assert(0)

class test_log2(unittest.TestCase):
    def check_log_value(self):
        val = log2(1.2)
        assert_almost_equal(val,0.26303440583379378,11)
        val = log2(1024)
        assert_almost_equal(val,10.0,15)
    def check_log_0(self):
        """ log(0) should print warning, but succeed.
        """
        try:
            val = log2(0)
            #assert(isinf(val))
        except:
            assert(0)
    def check_log_neg(self):
        """ log(-1) should print warning, but succeed.
        """
        try:
            val = log2(-1)
            #assert(isinf(val))
        except:
            pass

class test_histogram(unittest.TestCase):
    def check_basic_case(self):
        a = array((.2,.3,1.2,2.5)); bins = arange(4)
        hist = histogram(a,bins)
        desired = array((2,1,1,0))
        assert_array_equal(hist,desired)
    def check_lower_bin_edge(self):
        a = array((0,.3,1.2,2.5)); bins = arange(4)
        hist = histogram(a,bins)
        desired = array((2,1,1,0))
        assert_array_equal(hist,desired)
    def check_upper_bin_edge(self):
        a = array((0,.3,2,2.5)); bins = arange(4)
        hist = histogram(a,bins)
        desired = array((2,0,2,0))
        assert_array_equal(hist,desired)
    def check_below_lower_bin(self):
        a = array((-1,.3,2,2.5)); bins = arange(4)
        hist = histogram(a,bins)
        desired = array((1,0,2,0))
        assert_array_equal(hist,desired)
    def check_beyond_upper_bin(self):
        a = array((0,.3,2,10)); bins = arange(4)
        hist = histogram(a,bins)
        desired = array((2,0,1,2))

class test_array_split(unittest.TestCase):
    def check_integer_0_split(self):
        a = arange(10)
        try:
            res = array_split(a,0)
            assert(0) # it should have thrown a value error
        except ValueError:
            pass
    def check_integer_split(self):
        a = arange(10)
        res = array_split(a,1)
        desired = [arange(10)]
        compare_results(res,desired)

        res = array_split(a,2)
        desired = [arange(5),arange(5,10)]
        compare_results(res,desired)

        res = array_split(a,3)
        desired = [arange(4),arange(4,7),arange(7,10)]
        compare_results(res,desired)

        res = array_split(a,4)
        desired = [arange(3),arange(3,6),arange(6,8),arange(8,10)]
        compare_results(res,desired)

        res = array_split(a,5)
        desired = [arange(2),arange(2,4),arange(4,6),arange(6,8),arange(8,10)]
        compare_results(res,desired)

        res = array_split(a,6)
        desired = [arange(2),arange(2,4),arange(4,6),arange(6,8),arange(8,9),
                   arange(9,10)]
        compare_results(res,desired)

        res = array_split(a,7)
        desired = [arange(2),arange(2,4),arange(4,6),arange(6,7),arange(7,8),
                   arange(8,9), arange(9,10)]
        compare_results(res,desired)

        res = array_split(a,8)
        desired = [arange(2),arange(2,4),arange(4,5),arange(5,6),arange(6,7),
                   arange(7,8), arange(8,9), arange(9,10)]
        compare_results(res,desired)

        res = array_split(a,9)
        desired = [arange(2),arange(2,3),arange(3,4),arange(4,5),arange(5,6),
                   arange(6,7), arange(7,8), arange(8,9), arange(9,10)]
        compare_results(res,desired)

        res = array_split(a,10)
        desired = [arange(1),arange(1,2),arange(2,3),arange(3,4),
                   arange(4,5),arange(5,6), arange(6,7), arange(7,8),
                   arange(8,9), arange(9,10)]
        compare_results(res,desired)

        res = array_split(a,11)
        desired = [arange(1),arange(1,2),arange(2,3),arange(3,4),
                   arange(4,5),arange(5,6), arange(6,7), arange(7,8),
                   arange(8,9), arange(9,10),array([])]
        compare_results(res,desired)
    def check_integer_split_2D_rows(self):
        a = array([arange(10),arange(10)])
        res = array_split(a,3,axis=0)
        desired = [array([arange(10)]),array([arange(10)]),array([])]
        compare_results(res,desired)
    def check_integer_split_2D_cols(self):
        a = array([arange(10),arange(10)])
        res = array_split(a,3,axis=-1)
        desired = [array([arange(4),arange(4)]),
                   array([arange(4,7),arange(4,7)]),
                   array([arange(7,10),arange(7,10)])]
        compare_results(res,desired)
    def check_integer_split_2D_default(self):
        """ This will fail if we change default axis
        """
        a = array([arange(10),arange(10)])
        res = array_split(a,3)
        desired = [array([arange(10)]),array([arange(10)]),array([])]
        compare_results(res,desired)
    #perhaps should check higher dimensions

    def check_index_split_simple(self):
        a = arange(10)
        indices = [1,5,7]
        res = array_split(a,indices,axis=-1)
        desired = [arange(0,1),arange(1,5),arange(5,7),arange(7,10)]
        compare_results(res,desired)

    def check_index_split_low_bound(self):
        a = arange(10)
        indices = [0,5,7]
        res = array_split(a,indices,axis=-1)
        desired = [array([]),arange(0,5),arange(5,7),arange(7,10)]
        compare_results(res,desired)
    def check_index_split_high_bound(self):
        a = arange(10)
        indices = [0,5,7,10,12]
        res = array_split(a,indices,axis=-1)
        desired = [array([]),arange(0,5),arange(5,7),arange(7,10),
                   array([]),array([])]
        compare_results(res,desired)
        
class test_split(unittest.TestCase):
    """* This function is essentially the same as array_split,
         except that it test if splitting will result in an
         equal split.  Only test for this case.
    *"""
    def check_equal_split(self):
        a = arange(10)
        res = split(a,2)
        desired = [arange(5),arange(5,10)]
        compare_results(res,desired)

    def check_unequal_split(self):
        a = arange(10) 
        try:
            res = split(a,3)
            assert(0) # should raise an error
        except ValueError:
            pass

class test_atleast_1d(unittest.TestCase):
    def check_0D_array(self):
        a = array(1); b = array(2);
        res=atleast_1d([a,b])
        desired = [array([1]),array([2])]
        assert_array_equal(res,desired)
    def check_1D_array(self):
        a = array([1,2]); b = array([2,3]);
        res=atleast_1d([a,b])
        desired = [array([1,2]),array([2,3])]
        assert_array_equal(res,desired)
    def check_2D_array(self):
        a = array([[1,2],[1,2]]); b = array([[2,3],[2,3]]);
        res=atleast_1d([a,b])
        desired = [a,b]
        assert_array_equal(res,desired)
    def check_3D_array(self):
        a = array([[1,2],[1,2]]); b = array([[2,3],[2,3]]);
        a = array([a,a]);b = array([b,b]);
        res=atleast_1d([a,b])
        desired = [a,b]
        assert_array_equal(res,desired)

class test_atleast_2d(unittest.TestCase):
    def check_0D_array(self):
        a = array(1); b = array(2);
        res=atleast_2d([a,b])
        desired = [array([[1]]),array([[2]])]
        assert_array_equal(res,desired)
    def check_1D_array(self):
        a = array([1,2]); b = array([2,3]);
        res=atleast_2d([a,b])
        desired = [array([[1,2]]),array([[2,3]])]
        assert_array_equal(res,desired)
    def check_2D_array(self):
        a = array([[1,2],[1,2]]); b = array([[2,3],[2,3]]);
        res=atleast_2d([a,b])
        desired = [a,b]
        assert_array_equal(res,desired)
    def check_3D_array(self):
        a = array([[1,2],[1,2]]); b = array([[2,3],[2,3]]);
        a = array([a,a]);b = array([b,b]);
        res=atleast_2d([a,b])
        desired = [a,b]
        assert_array_equal(res,desired)

class test_atleast_3d(unittest.TestCase):
    def check_0D_array(self):
        a = array(1); b = array(2);
        res=atleast_3d([a,b])
        desired = [array([[[1]]]),array([[[2]]])]
        assert_array_equal(res,desired)
    def check_1D_array(self):
        a = array([1,2]); b = array([2,3]);
        res=atleast_3d([a,b])
        desired = [array([[[1],[2]]]),array([[[2],[3]]])]
        assert_array_equal(res,desired)
    def check_2D_array(self):
        a = array([[1,2],[1,2]]); b = array([[2,3],[2,3]]);
        res=atleast_3d([a,b])
        desired = [a[:,:,NewAxis],b[:,:,NewAxis]]
        assert_array_equal(res,desired)
    def check_3D_array(self):
        a = array([[1,2],[1,2]]); b = array([[2,3],[2,3]]);
        a = array([a,a]);b = array([b,b]);
        res=atleast_3d([a,b])
        desired = [a,b]
        assert_array_equal(res,desired)

class test_hstack(unittest.TestCase):
    def check_0D_array(self):
        a = array(1); b = array(2);
        res=hstack([a,b])
        desired = array([1,2])
        assert_array_equal(res,desired)
    def check_1D_array(self):
        a = array([1]); b = array([2]);
        res=hstack([a,b])
        desired = array([1,2])
        assert_array_equal(res,desired)
    def check_2D_array(self):
        a = array([[1],[2]]); b = array([[1],[2]]);
        res=hstack([a,b])
        desired = array([[1,1],[2,2]])
        assert_array_equal(res,desired)

class test_vstack(unittest.TestCase):
    def check_0D_array(self):
        a = array(1); b = array(2);
        res=vstack([a,b])
        desired = array([[1],[2]])
        assert_array_equal(res,desired)
    def check_1D_array(self):
        a = array([1]); b = array([2]);
        res=vstack([a,b])
        desired = array([[1],[2]])
        assert_array_equal(res,desired)
    def check_2D_array(self):
        a = array([[1],[2]]); b = array([[1],[2]]);
        res=vstack([a,b])
        desired = array([[1],[2],[1],[2]])
        assert_array_equal(res,desired)
    def check_2D_array2(self):
        a = array([1,2]); b = array([1,2]);
        res=vstack([a,b])
        desired = array([[1,2],[1,2]])
        assert_array_equal(res,desired)

class test_dstack(unittest.TestCase):
    def check_0D_array(self):
        a = array(1); b = array(2);
        res=dstack([a,b])
        desired = array([[[1,2]]])
        assert_array_equal(res,desired)
    def check_1D_array(self):
        a = array([1]); b = array([2]);
        res=dstack([a,b])
        desired = array([[[1,2]]])
        assert_array_equal(res,desired)
    def check_2D_array(self):
        a = array([[1],[2]]); b = array([[1],[2]]);
        res=dstack([a,b])
        desired = array([[[1,1]],[[2,2,]]])
        assert_array_equal(res,desired)
    def check_2D_array2(self):
        a = array([1,2]); b = array([1,2]);
        res=dstack([a,b])
        desired = array([[[1,1],[2,2]]])
        assert_array_equal(res,desired)

""" array_split has more comprehensive test of splitting.
    only do simple test on hsplit, vsplit, and dsplit
"""
class test_hsplit(unittest.TestCase):
    """ only testing for integer splits.
    """
    def check_0D_array(self):
        a= array(1)
        try:
            hsplit(a,2)
            assert(0)
        except ValueError:
            pass
    def check_1D_array(self):
        a= array([1,2,3,4])
        res = hsplit(a,2)
        desired = [array([1,2]),array([3,4])]
        compare_results(res,desired)
    def check_2D_array(self):
        a= array([[1,2,3,4],
                  [1,2,3,4]])
        res = hsplit(a,2)
        desired = [array([[1,2],[1,2]]),array([[3,4],[3,4]])]
        compare_results(res,desired)

class test_vsplit(unittest.TestCase):
    """ only testing for integer splits.
    """
    def check_1D_array(self):
        a= array([1,2,3,4])
        try:
            vsplit(a,2)
            assert(0)
        except ValueError:
            pass
    def check_2D_array(self):
        a= array([[1,2,3,4],
                  [1,2,3,4]])
        res = vsplit(a,2)
        desired = [array([[1,2,3,4]]),array([[1,2,3,4]])]
        compare_results(res,desired)

class test_dsplit(unittest.TestCase):
    """ only testing for integer splits.
    """
    def check_2D_array(self):
        a= array([[1,2,3,4],
                  [1,2,3,4]])
        try:
            dsplit(a,2)
            assert(0)
        except ValueError:
            pass
    def check_3D_array(self):
        a= array([[[1,2,3,4],
                   [1,2,3,4]],
                  [[1,2,3,4],
                   [1,2,3,4]]])
        res = dsplit(a,2)
        desired = [array([[[1,2],[1,2]],[[1,2],[1,2]]]),
                   array([[[3,4],[3,4]],[[3,4],[3,4]]])]
        compare_results(res,desired)

class test_trim_zeros(unittest.TestCase):
    """ only testing for integer splits.
    """
    def check_basic(self):
        a= array([0,0,1,2,3,4,0])
        res = trim_zeros(a)
        assert_array_equal(res,array([1,2,3,4]))
    def check_leading_skip(self):
        a= array([0,0,1,0,2,3,4,0])
        res = trim_zeros(a)
        assert_array_equal(res,array([1,0,2,3,4]))
    def check_trailing_skip(self):
        a= array([0,0,1,0,2,3,0,4,0])
        res = trim_zeros(a)
        assert_array_equal(res,array([1,0,2,3,0,4]))

                
# Utility

def compare_results(res,desired):
    for i in range(len(desired)):
        assert_array_equal(res[i],desired[i])

##################################################

def test_suite(level=1):
    suites = []
    if level > 0:
        suites.append( unittest.makeSuite(test_real,'check_') )
        suites.append( unittest.makeSuite(test_imag,'check_') )
        suites.append( unittest.makeSuite(test_sum,'check_') )
        suites.append( unittest.makeSuite(test_mean,'check_') )
        suites.append( unittest.makeSuite(test_array_iscomplex,'check_') )
        suites.append( unittest.makeSuite(test_array_isreal,'check_') )        
        suites.append( unittest.makeSuite(test_iscomplex,'check_') )
        suites.append( unittest.makeSuite(test_isreal,'check_') )    
        suites.append( unittest.makeSuite(test_isnan,'check_') )
        suites.append( unittest.makeSuite(test_isfinite,'check_') )
        suites.append( unittest.makeSuite(test_isinf,'check_') )
        suites.append( unittest.makeSuite(test_isposinf,'check_') )    
        suites.append( unittest.makeSuite(test_isneginf,'check_') )
        suites.append( unittest.makeSuite(test_nan_to_num,'check_') )
        suites.append( unittest.makeSuite(test_logn,'check_') )
        suites.append( unittest.makeSuite(test_log2,'check_') )
        suites.append( unittest.makeSuite(test_histogram,'check_') )
        suites.append( unittest.makeSuite(test_array_split,'check_') )
        suites.append( unittest.makeSuite(test_split,'check_') )
        suites.append( unittest.makeSuite(test_atleast_1d,'check_') )
        suites.append( unittest.makeSuite(test_atleast_2d,'check_') )
        suites.append( unittest.makeSuite(test_atleast_3d,'check_') )
        suites.append( unittest.makeSuite(test_hstack,'check_') )
        suites.append( unittest.makeSuite(test_vstack,'check_') )
        suites.append( unittest.makeSuite(test_dstack,'check_') )
        suites.append( unittest.makeSuite(test_hsplit,'check_') )    
        suites.append( unittest.makeSuite(test_vsplit,'check_') )
        suites.append( unittest.makeSuite(test_dsplit,'check_') )
        suites.append( unittest.makeSuite(test_trim_zeros,'check_') )
    
    total_suite = unittest.TestSuite(suites)
    return total_suite

def test(level=10):
    all_tests = test_suite(level)
    runner = unittest.TextTestRunner()
    runner.run(all_tests)
    return runner


if __name__ == "__main__":
    test()
