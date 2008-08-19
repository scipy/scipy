# P for pathological
import unittest
import numpy as np
from numpy import array
import dewall as dw
import interpolateSNd as SNd
reload(dw)

class Test(unittest.TestCase):
    
    def compare_array(self, a, b):
        return np.allclose(a,b)
    
    ## test Delaunay triangulation itself
    
    # testing pathological cases with
    
    def test_square(self):
        P = [array([0.,1.]), array([0.,0.]), array([1.,0.]), array([1.,1.])]
        tri = dw.dewall(P)
        self.assert_( len(tri)==2 )
        self.assert_( len(tri[0])==3 )
        self.assert_( len(tri[1])==3 )
        print "square triangulation:\n", tri
        
    def test_linear(self):
        P = [array([0.,1.]), array([0.,0.]), array([0.,-1.])]
        tri = dw.dewall(P)
        print "line triang:\n", tri
        
    # testing general case using random data
        
    ## test interpolation, and thus also triangulation by extension
    def test_linear_on_cube(self):
        x = array([0., 1, 0, 1, 0, 1, 0, 1])
        y = array([0., 0, 1, 1, 0, 0, 1, 1])
        z = array([0., 0, 0, 0, 1, 1, 1, 1])
        points = array([x,y,z]).reshape((3,8))
        fvals = x+y-z
        
        interp = SNd.InterpolateSNd(points, fvals)
        
        newdata = np.random.random_sample((3,20))
        interpvals = interp(newdata)
        
        realvals = newdata[0,:]+newdata[1,:]-newdata[2,:]
        
        self.assert_(compare_array(np.ravel(interpvals), np.ravel(realvals)))
        
        
        
if __name__ == "__main__":
    unittest.main()