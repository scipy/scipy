# P for pathological
import unittest
import numpy as np
from numpy import array
import dewall as dw
import interpolateSNd as SNd
reload(dw)

class Test(unittest.TestCase):
    
    def compare_arrays(self, a, b):
        return np.allclose(a,b,rtol=1e-3) or (np.isnan(a)&np.isnan(b)).all()
    
    ## test Delaunay triangulation itself
    
    # testing pathological cases with
    
    def _test_square(self):
        P = [array([0.,1.]), array([0.,0.]), array([1.,0.]), array([1.,1.])]
        tri = dw.dewall(P)
        self.assert_( len(tri)==2 )
        self.assert_( len(tri[0])==3 )
        self.assert_( len(tri[1])==3 )
        
    def _test_linear(self):
        P = [array([0.,1.]), array([0.,0.]), array([0.,-1.])]
        tri = dw.dewall(P)
        
    # testing general case using random data
    def _test_2d(self):
        ndim = 2
        nold = 15
        print "TESTING %iD"%ndim
        corner_points_matrix = corners(ndim)
        corner_points = [corner_points_matrix[:,i] for i in range(2**ndim)]
        interior_points = [np.random.random_sample(ndim) for i in range(nold)]
        P = corner_points+interior_points
        
        # not checking if its correct, just if it runs.
        tri = dw.dewall(P)
        
        # making sure it's correct
    
    def _test_3d(self):
        print "TESTING 3D"
        P = [np.random.random_sample(3) for i in range(9)]
        tri = dw.dewall(P)
        # not checking if its correct, just if it runs.
        
    def _test_4d(self):
        print "TESTING 4D"
        P = [np.random.random_sample(4) for i in range(9)]
        tri = dw.dewall(P)
        # not checking if its correct, just if it runs.
        
    def _test_5d(self):
        P = [np.random.random_sample(5) for i in range(9)]
        tri = dw.dewall(P)
        # not checking if its correct, just if it runs.
        
    ## test interpolation, and thus also triangulation by extension
    def _test_linear_on_cube(self):
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
    
    def test_linear_2d(self):
        ndim = 2 # num dimensions
        nold = 5 # num known data points
        nnew = 5 # num points at which to interpolate
        
        print "%iD Interpolation"%ndim
        
        P = [np.random.random_sample(ndim) for i in range(nold)]
        # points at corners of hypercube and radnimly scattered in the interior
        points = np.concatenate((corners(ndim) , array(P).reshape((ndim,nold))), axis=1)
        fvals = np.zeros((1,points.shape[1]))
        for i in range(ndim):
            fvals = fvals+points[i,:]
        fvals = fvals.reshape((points.shape[1]))
            
        print "points:\n",points
        print "fvals:\n",fvals
        
        interp = SNd.InterpolateSNd(points, fvals)
        
        print "\ntriang:"
        for x in interp._triangulation: print x
        
        newdata = np.random.random_sample((ndim,nnew))
        print "\nnewdata:\n",newdata
        interpvals = interp(newdata)
        realvals = np.sum(newdata, axis=0)
        
        print "%iD interpvals: "%ndim, interpvals
        print "%iD realvals:   "%ndim, realvals
        
        #self.assert_(self.compare_arrays(np.ravel(interpvals), np.ravel(realvals)))
        assert self.compare_arrays(np.ravel(interpvals), np.ravel(realvals)), "wrong data"

    def test_linear_3d(self):
        ndim = 3 # num dimensions
        nold = 1 # num known data points
        nnew = 5 # num points at which to interpolate
        
        print "%iD Interpolation"%ndim
        
        P = [np.random.random_sample(ndim) for i in range(nold)]
        # points at corners of hypercube and radnimly scattered in the interior
        points = np.concatenate((corners(ndim) , array(P).reshape((ndim,nold))), axis=1)
        fvals = np.zeros((1,points.shape[1]))
        for i in range(ndim):
            fvals = fvals+points[i,:]
        fvals = fvals.reshape((points.shape[1]))
            
        #print "points:\n",points
        #print "fvals:\n",fvals
        
        interp = SNd.InterpolateSNd(points, fvals)
        
        #print "\ntriang:"
        #for x in interp._triangulation: print x
        
        newdata = np.random.random_sample((ndim,nnew))
        #print "\nnewdata:\n",newdata
        interpvals = interp(newdata)
        realvals = np.sum(newdata, axis=0)
        
        print "%iD interpvals: "%ndim, interpvals
        print "%iD realvals:   "%ndim, realvals
        
        #self.assert_(self.compare_arrays(np.ravel(interpvals), np.ravel(realvals)))
        assert self.compare_arrays(np.ravel(interpvals), np.ravel(realvals)), "wrong data"

def corners(ndim):
    # returns matrix indicating corners of unit cube
    result = np.zeros((ndim,2**ndim))
    for i in range(ndim):
        # each block is 2**i spaces wide
        for j in range(2**ndim):
            if ( j%(2**(i+1)) )/(2**i) == 1: result[i,j]=1.
    return result
        
        
if __name__ == "__main__":
    unittest.main()