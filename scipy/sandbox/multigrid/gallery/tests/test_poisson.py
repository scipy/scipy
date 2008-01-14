from scipy.testing import *

from scipy import matrix

from scipy.sandbox.multigrid.gallery.poisson import *

class TestPoisson(TestCase):
    def test_3pt(self):
        cases = []

        cases.append( (1,matrix([[2]])) )
        cases.append( (2,matrix([[ 2,-1],
                                 [-1, 2]])) )
        cases.append( (4,matrix([[ 2,-1, 0, 0],
                                 [-1, 2,-1, 0],
                                 [ 0,-1, 2,-1],
                                 [ 0, 0,-1, 2]])) )


        for N,expected in cases:
            assert_equal( poisson(N,stencil='3pt').todense(), expected )



if __name__ == '__main__':
    unittest.main()
