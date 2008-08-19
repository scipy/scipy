""" ND Scattered interpolation
"""

import dewall as dw
import numpy as np

class InterpolateSNd:
    """ Interpolation of scatter data in arbitrarily many dimensions
        
        *******
        WARNING : this code is currently extremely slow.  It is very much still
        *******
                        in development, but because it functions correctly it is being
                        presented.
                        
        Parameters:
        ------------
        
            points : 2D numpy array
                    dxn array, where d is dimensionality of the data and
                    there are n data points.  Each column of the array is
                    a data point.
                    
            fvals : 1D numpy array
                    Must have length n.  fvals[j] is the value of the function
                    at point point[:,j]
                    
            newX (when calling) : 2D numpy array
                    shape dxm.  Each column is a point at which to interpolate
                    the function value.
                    
        The valid interpolation range is the convex hull of the points.  All in-range
        points will be linearly interolated, and all out-of-bounds points will return
        NaN.
    """
    
    def __init__(self, P, fvals):
        # P = array of points, where each column is a point
        #       points must be distinct
        # fvals = array of known function values
        
        assert P.ndim == 2, "P must be 2-dimensional"
        d, n = P.shape
        assert len(fvals)==n, \
            "fvals must have length n, where n is number of points"
        
        # remember dimensionality of space
        self.ndim = d
                    
        # store list of points
        self.points = P.copy()
        
        # store function values
        self.fvals = fvals.copy()
        
        # calculate and store delaunay triangulation
        list_of_points = [P[:,i].reshape(d) for i in range(n)]
        self._triangulation = dw.dewall(list_of_points)
        
        # for each simplex, indices of its points in self.points
        # lists of tuples(lists) of integers
        indices_list = []
        for simplex in self._triangulation:
            indices_for_current_point = []
            for point in simplex:
                for i in range(n):
                    # if ith col = point, append i to list
                    if (point.reshape((d)) == np.reshape(P[:,i],(d))).all():
                        indices_for_current_point.append(i)
            indices_list.append(indices_for_current_point)
        self._indices_list = indices_list
        
        # FIXME : store list of (simplex, vals) tuples ????
            
    def __call__(self, X):
        # X = dxm array of points at which to interpolate
        d, m = X.shape
        assert d==self.ndim, "input must have correct dimensionality"
        
        # break input into a list of points for list comprehensions.
        # yeah, it's very slow, but it's elegant for now
        list_of_points = [X[:,i].reshape((d,1)) for i in range(m)]
                                        
        #  Tuples of form (simplex , vals_at_vertices, point)
        simplices_vals_points = [(self._get_simplex(pt), self._get_fvals(pt), pt) \
                                        for pt in list_of_points]
        
        result = [self.linear_interp( simplex, vals, point ) \
                        for simplex, vals, point in simplices_vals_points]
                            
        return result
        
    def _get_fvals(self, pt):
        assert pt.shape == (self.ndim, 1), "point shape: "+str(pt.shape)
        indices = self._get_simplex_indices(pt)
        if indices==None:
            return None
        return self.fvals[indices]
                
    def _get_simplex(self, pt):
        assert pt.shape == (self.ndim, 1), "point shape: "+str(pt.shape)
        indices = self._get_simplex_indices(pt)
        if indices==None:
            return None
        return self.points[:,indices]
        
    def _get_simplex_indices(self, pt):
        # returns tuple indicating vertices of simplex containing pt
        assert pt.shape == (self.ndim, 1), "point shape: "+str(pt.shape)
        indices = None
        for i in range(len(self._indices_list)):
            if self.point_is_in_simplex(pt.reshape((self.ndim,1)), self.points[:,self._indices_list[i]]):
                indices = self._indices_list[i]
                break
        return indices
        
    def point_is_in_simplex(self, point, simplex):
        # point = array
        # simplex = matrix
        assert point.shape == (self.ndim, 1), "wrong shape of point: "+str(point.shape)
        assert simplex.shape == (self.ndim, self.ndim+1), "wrong shape of simplex: "+str(simplex.shape)
        weights_vec = self.calculate_weights(simplex, point)
        weight_in_range = (0<= weights_vec) & (weights_vec <= 1)
        return np.alltrue(weight_in_range, axis=0)
    
    def linear_interp(self, simplex, vals, point):
        if simplex == None or vals == None: # point is out of range
            return np.NaN
        
        # testing
        assert point.shape == (self.ndim, 1), \
                "wrong shape of point: "+str(point.shape)
        assert simplex.shape == (self.ndim, self.ndim+1), \
                "wrong shape of simplex: "+str(simplex.shape)
        assert vals.shape == (self.ndim+1,), \
                "vals wrong shape: "+str(vals.shape)+"need: "+str((self.ndim,))
        
        weights = self.calculate_weights(simplex, point)
        return np.dot(np.ravel(weights), np.ravel(vals))
    
    def calculate_weights(self, simplex, points):
        """ Each column in points is a weighted average
            of columns in simplex.  Returns matrix where
            jth column gives these weights for the jth column
            of points.
        """
        assert simplex.shape == (self.ndim, self.ndim+1), "simplex shape: "+str(simplex.shape)
        d, V = simplex.shape
        d, P = points.shape
        
        matrix_to_solve = np.concatenate((simplex, \
                                                        np.ones((1,V))
                                                        ),
                                                        axis = 0
                                                        )
        vec_to_solve = np.concatenate((points, \
                                                        np.ones((1,P))
                                                        ),
                                                        axis = 0
                                                        )
        
        weights_vecs = np.linalg.solve(matrix_to_solve, vec_to_solve)
        
        return weights_vecs
        
def point_is_in_simplex(point, simplex):
    # point = array
    # simplex = matrix
    #assert point.shape == (self.ndim, 1), "wrong shape of point: "+str(point.shape)
    #assert simplex.shape == (self.ndim, self.ndim+1), "wrong shape of simplex: "+str(simplex.shape)
    weights_vec = calculate_weights(simplex, point)
    print "weights:\n", weights_vec
    weight_in_range = (0.<= weights_vec) & \
                                    (weights_vec <= 1.)
    print "in_range:\n", weight_in_range
    return np.alltrue(weight_in_range, axis=0)
    
def calculate_weights(simplex, points):
        """ Each column in points is a weighted average
            of columns in simplex.  Returns matrix where
            jth column gives these weights for the jth column
            of points.
        """
        N, V = simplex.shape
        N, P = points.shape
        
        matrix_to_solve = np.concatenate((simplex, \
                                                        np.ones((1,V))
                                                        ),
                                                        axis = 0
                                                        )
        vec_to_solve = np.concatenate((points, \
                                                        np.ones((1,P))
                                                        ),
                                                        axis = 0
                                                        )
        
        weights_vecs = np.linalg.solve(matrix_to_solve, vec_to_solve)
        
        return weights_vecs