# P for pathological
import unittest
import numpy as np
from numpy import array
import dewall as dw
import interpolateSNd as SNd
reload(dw)

class Test(unittest.TestCase):
    
    def compare_arrays(self, a, b):
        return np.allclose(a,b)
    
    ## test Delaunay triangulation itself
    
    # testing pathological cases with
    
    def test_square(self):
        P = [array([0.,1.]), array([0.,0.]), array([1.,0.]), array([1.,1.])]
        tri = dw.dewall(P)
        self.assert_( len(tri)==2 )
        self.assert_( len(tri[0])==3 )
        self.assert_( len(tri[1])==3 )
        
    def test_linear(self):
        P = [array([0.,1.]), array([0.,0.]), array([0.,-1.])]
        tri = dw.dewall(P)
        
    # testing general case using random data
    
    def test_2d(self):
        print "TESTING 2D"
        P = [np.random.random_sample(2) for i in range(15)]
        tri = dw.dewall(P)
        # not checking if its correct, just if it runs.
    
    def test_3d(self):
        print "TESTING 3D"
        P = [np.random.random_sample(3) for i in range(9)]
        tri = dw.dewall(P)
        # not checking if its correct, just if it runs.
        
    def test_4d(self):
        print "TESTING 4D"
        P = [np.random.random_sample(4) for i in range(9)]
        tri = dw.dewall(P)
        # not checking if its correct, just if it runs.
        
    def test_5d(self):
        P = [np.random.random_sample(5) for i in range(9)]
        tri = dw.dewall(P)
        # not checking if its correct, just if it runs.
        
    ## test interpolation, and thus also triangulation by extension
    
    def test_2d(self):
        points = np.array([[ 1.,0.,1.,0.],[0.,0.,1.,1.]])
        z = np.array([1.,0.,2.,1.])
        interp = SNd.InterpolateSNd(points,z)

        X=np.array([[.4,.1,.55],[.4,.1,.3]])

        output = interp(X)
        self.assert_(self.compare_arrays(output, array([.8,.2,.85])))
    
    def test_linear_on_cube(self):
        x = array([0., 1, 0, 1, 0, 1, 0, 1])
        y = array([0., 0, 1, 1, 0, 0, 1, 1])
        z = array([0., 0, 0, 0, 1, 1, 1, 1])
        points = array([x,y,z]).reshape((3,8))
        fvals = x+y-z
        
        interp = SNd.InterpolateSNd(points, fvals)
        
        newdata = np.random.random_sample((3,8))
        interpvals = interp(newdata)
        realvals = newdata[0,:]+newdata[1,:]-newdata[2,:]
        
        self.assert_(self.compare_arrays(np.ravel(interpvals), np.ravel(realvals)))
    
    def test_linear_5d(self):
        P = [np.random.random_sample(5) for i in range(9)]
        points = array(P).reshape((5,9))
        fvals = points[:,0]+points[:,1]+points[:,2]+points[:,3]+points[:,4]
        
        interp = SNd.InterpolateSNd(points, fvals)
        
        newdata = np.random.random_sample((5,8))
        interpvals = interp(newdata)
        realvals = np.sum(newdata, axis=0)
        
        self.assert_(self.compare_arrays(np.ravel(interpvals), np.ravel(realvals)))
        
if __name__ == "__main__":
    unittest.main()