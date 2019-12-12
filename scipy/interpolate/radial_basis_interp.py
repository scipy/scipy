import sys
if sys.version_info[0] > 2:
    unicode = str

import numpy as np
import scipy.spatial
from scipy.special import xlogy
import scipy.linalg

class RadialBasisFunction(object):
    """
    Augmented Radial Basis function interpolation with (optional) polynomial 
    augmentation. The construction follows [1] with code borrowed from [2] 
    for the kernel (which they call the "function")
    
    Parameters
    ----------
    X : (Npoints, Ndims) ndarray of floats
        Data point coordinates.

    f : (Npoints,) ndarray of float or complex
        Data values.
    
    kernel : string, callable, or int: optional
        Default: gaussian
        The radial basis kernel, based on the radius, r
        If string:
            'multiquadric': sqrt((r/epsilon)**2 + 1)
            'inverse': 1.0/sqrt((r/epsilon)**2 + 1)
            'gaussian','squared_exponential','rbf': exp(-(r/epsilon)**2)
            'linear': r
            'cubic': r**3
            'quintic': r**5
            'thin_plate': r**2 * log(r)
        If callable
            must take the arguments (r,epsilon) and ignore epsilon if 
            not needed. 
        If integer:
            Spline of the form: r**kernel if kernel is odd else r**kernel*log(r)
    
    degree : int: optional
        Use an augmented RBF. Specify the *total* degree fit if given as an 
        integer. If given as a matrix, must be (P+1,ndim) shape. This is the 
        degree of the augmenting polynomial. Default: 0            

    epsilon: float or None: optional
        Adjustable parameter for some kernels. Default is the average
        interpoint distance for X (it is computed anyway).
        
        See [3] and [4] for alternative methods to select epsilon
        
    smooth : float: optional
        Values greater than zero increase the smoothness of the
        approximation.  0 is for interpolation (default). Actually set to
        max(smooth,1e-10 * ||f||) for numerical precision.

    Methods:
    --------
    __call__(X)
        Evaluate the interpolant at X

    Tips:
    -----
    While not required, it is *highly* suggested to scale all inputs to a unit
    hypercube! This will improve the quality of the interpolant; especially for
    cases with large differences in scales.

    Theory:
    ------
    The RBF is based on a radial distance 
    
        r_j = ||x - x_j||                                               (1)
    
    with the kernel phi(r) and an optional "Agumenting" polynomial of set 
    degree. The interpolant is of the form: [1]
    
        s(x) = sum_{j=1}^N c_j * phi(r_j) + sum_{k \in P} g_k * p_k     (2)
    
    where `P` is a polynomial index and `p_k` is the corresponding
    polynomial (e.g. k = {1,2,0} then p_k = x*y^2*z^0).
    
    In order to solve for c_j and g_k we enforce 
        (a) s(x_j) = f_j for all j (i.e. exact interpolation)
        (b) sum_{j=1}^N c_j*p_k(x_j) = 0 for all k=1,...,N 
    
    That results in a matrix system:
        ____________________________________________________   ___     ___
       |                            |                       | |   |   |   |
       |   phi(||x_i - x_j||) (N,N) |   p_k(x_j) (N,len(K)) | | c |   | f |
       |                            |                       | |   | = |   | (3)
       |----------------------------|-----------------------| |---|   |---|
       |       p_k(x_j).T           |           0           | | g |   | 0 |
        ----------------------------------------------------   ---     ---
      
    which is then solved. Given c and g, (2) can be used to solve for
    any given x
    
    References:
    ----------  
    [1] G. B. Wright. Radial Basis Function Interpolation: Numerical and 
        Analytical Developments. PhD thesis, University of Colorado, 2003.

    [2] SciPy version 1.1.0
        https://github.com/scipy/scipy/blob/v1.1.0/scipy/interpolate/rbf.py
    
    [3] S. Rippa. An algorithm for selecting a good value for the parameter 
        c in radial basis function interpolation. Advances in Computational 
        Mathematics, 11(2-3):193–210, 1999.
        
    [4] J. D. Martin and T. W. Simpson. Use of Kriging Models to Approximate 
        Deterministic Computer Models. AIAA journal, 43(4):853–863, 2005.

    """
    def __init__(self,X,f,
                 kernel='gaussian',
                 degree=0,
                 epsilon=1,
                 smooth=0):
        self.X = X = self._checkX(X)
             
        self.N,self.ndim = X.shape
        self.f = f = np.ravel(f)
        self.degree = degree
        self.kernel = kernel
        
        r = scipy.spatial.distance.cdist(X,X)
        
        if epsilon is None:
            epsilon = np.mean(r[np.triu_indices(self.N,1)])
        self.epsilon = epsilon
        
        self.smooth = max(float(smooth),1e-10*scipy.linalg.norm(f)/np.sqrt(len(f))) # Scale by ||f||
        
        K = self._kernel(r) + np.eye(self.N)*self.smooth
        P = RadialBasisFunction.vandermond(X,degree=self.degree)
        
        # Build the matrix
        Z = np.zeros([P.shape[1]]*2)
        z = np.zeros(P.shape[1])
        KP = np.block([[K   , P],
                       [P.T , Z]])
        b = np.hstack([f,z])
        coef = scipy.linalg.solve(KP,b)
        self.rbf_coef = coef[:self.N]
        self.poly_coef = coef[self.N:]
        
    def _checkX(self,X):
        X = np.asarray(X)
        if X.ndim == 1:
            X = np.c_[X]
        elif X.ndim != 2:
            raise ValueError('X must be 1D or 2D')
        return X

    def __call__(self,X):
        X = self._checkX(X)
        K = self._kernel( scipy.spatial.distance.cdist(X,self.X) )
        P = RadialBasisFunction.vandermond(X,degree=self.degree)
        return K.dot(self.rbf_coef) + P.dot(self.poly_coef)

    def _kernel(self,r):
        r = np.asarray(r)
        if callable(self.kernel):
            return self.kernel(r,self.epsilon)
        elif isinstance(self.kernel,(int,np.integer)):
            if np.mod(self.kernel,2) == 0:
                return xlogy(r**self.kernel,r)
            else:
                return r**self.kernel
        elif not isinstance(self.kernel,(str,unicode)):
            raise ValueError('Kernel must be a callable with signature (r,epsilon), an integer (spline) or a valid string')
        
        kernel = self.kernel.lower().replace(' ','_').replace('-','_')
        
        if kernel in ['multiquadric']:
            return np.sqrt((r/self.epsilon)**2 + 1)
        elif kernel in ['inverse_multiquadric','inverse']:
            return 1.0/np.sqrt((1.0/self.epsilon*r)**2 + 1)
        elif kernel in ['gaussian','squared_exponential','rbf']:
            return np.exp(-(1.0/self.epsilon*r)**2)
        elif kernel in ['linear']:
            return r
        elif kernel in ['cubic']:
            return r**3
        elif kernel in ['quintic']:
            return r**5
        elif kernel in ['thin_plate']:
            return xlogy(r**2, r)
        else:
            raise ValueError('Not valid kernel name')
    
    @staticmethod
    def vandermond(X,degree):
        """
        Return a Vandermond matrix of X up to the
        *total* order degree
        """
        X = np.atleast_2d(X)
        ndim = X.shape[1]
        index = np.asarray(RadialBasisFunction.total_index(degree,ndim),dtype=int)
        V = np.ones([X.shape[0],len(index)])
        for d in range(ndim):
            v = np.fliplr(np.vander(X[:,d],N=degree+1))
            V *= v[:,index[:,d]]
        return V
    
    @staticmethod
    def total_index(P,ndim,sort=True):
        P = tuple([P]*ndim)
        curr = [0]*ndim
        ind = []
        ind.append(tuple(curr))
        while True:
            for d in range(ndim):
                if sum(curr) < P[0]:
                    curr[d] += 1
                    break
                else:
                    curr[d] = 0
            else:
                break
            ind.append(tuple(curr))
        if sort:
            ind.sort(key=lambda a:(sum(a),(np.array(a)**2).sum(),a[::-1]))
        return ind
