from scipy.testing import *

from scipy import matrix

from scipy.sandbox.multigrid.gallery.laplacian import *

class TestPoisson(TestCase):
    def test_poisson(self):
        cases = []

        # 1D
        cases.append( ((1,),matrix([[2]])) )
        cases.append( ((2,),matrix([[ 2,-1],
                                    [-1, 2]])) )
        cases.append( ((4,),matrix([[ 2,-1, 0, 0],
                                    [-1, 2,-1, 0],
                                    [ 0,-1, 2,-1],
                                    [ 0, 0,-1, 2]])) )
        
        # 2D        
        cases.append( ((1,1), matrix([[4]])) )
        cases.append( ((2,1), matrix([[ 4,-1],
                                      [-1, 4]])) )
        cases.append( ((1,2), matrix([[ 4,-1],
                                      [-1, 4]])) )
        cases.append( ((1,3), matrix([[ 4,-1, 0],
                                      [-1, 4,-1],
                                      [ 0,-1, 4]])) )
        cases.append( ((2,2), matrix([[ 4,-1,-1, 0],
                                      [-1, 4, 0,-1],
                                      [-1, 0, 4,-1],
                                      [ 0,-1,-1, 4]])) )
        # 3D
        cases.append( ((2,2,1), matrix([[ 6,-1,-1, 0],
                                        [-1, 6, 0,-1],
                                        [-1, 0, 6,-1],
                                        [ 0,-1,-1, 6]])) )
        cases.append( ((2,2,2), matrix([[ 6,-1,-1, 0,-1, 0, 0, 0],
                                        [-1, 6, 0,-1, 0,-1, 0, 0],
                                        [-1, 0, 6,-1, 0, 0,-1, 0],
                                        [ 0,-1,-1, 6, 0, 0, 0,-1],
                                        [-1, 0, 0, 0, 6,-1,-1, 0],
                                        [ 0,-1, 0, 0,-1, 6, 0,-1],
                                        [ 0, 0,-1, 0,-1, 0, 6,-1],
                                        [ 0, 0, 0,-1, 0,-1,-1, 6]])) )

        for grid,expected in cases:
            assert_equal( poisson(grid).todense(), expected )

if __name__ == '__main__':
    nose.run(argv=['', __file__])
