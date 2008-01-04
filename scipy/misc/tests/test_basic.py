''' Test functions for factorial, comb functions '''

from scipy.testing import *
import numpy as np
from scipy.misc import factorial, comb


class test_factorial(TestCase):
    def test_basic(self):
        for k in range(0,13):
            assert_equal(factorial(k),
                         np.product(np.arange(1,k+1),axis=0))
            assert_equal(factorial(-10), 0)
        
    def test_exact(self):
        resdict = {4:24L,10:3628800L,15:1307674368000L,
                   19:121645100408832000L,
                   100:93326215443944152681699238856266700490715968264381621468592963895217599993229915608941463976156518286253697920827223758251185210916864000000000000000000000000L}
        for key in resdict.keys():
            assert_equal(factorial(key,exact=1),resdict[key])


class test_comb(TestCase):
    
    def test_basic(self):
        for N in range(0,11):
            for k in range(0,N+1):
                ans = np.product(np.arange(N-k+1,N+1),axis=0) \
                      / np.product(np.arange(1,k+1),axis=0)
                assert_almost_equal(comb(N,k),ans,9)
        assert(comb(-10,1) == 0)
        assert(comb(10,-1) == 0)
        assert(comb(-10,-3) == 0)
        assert(comb(10,11) == 0)

    def test_exact(self):
        resdict = {(10,2):45L, (10,5):252L,
                   (1000,20):339482811302457603895512614793686020778700L,
                   (1000,975):47641862536236518640933948075167736642053976275040L
                   }
        for key in resdict.keys():
            assert_equal(comb(key[0],key[1],exact=1),resdict[key])


if __name__ == '__main__':
    unittest.main()
