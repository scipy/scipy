#!/usr/bin/env python
"""
rbf - Radial basis functions for interpolation/smoothing scattered Nd data.

Written by John Travers <jtravs@gmail.com>, February 2007
Based closely on Matlab code by Alex Chirokov

Permission to use, modify, and distribute this software is given under the
terms of the SciPy (BSD style) license.  See LICENSE.txt that came with
this distribution for specifics.

NO WARRANTY IS EXPRESSED OR IMPLIED.  USE AT YOUR OWN RISK.

"""

import scipy as s
import scipy.linalg

class Rbf(object):
    """ A class for radial basis function approximation/interpolation of
        n-dimensional scattered data.
    """
    def __init__(self,x,y, function='multiquadrics', constant=None, smooth=0):
        """ Constructor for Rbf class.

            Inputs:
                x   (dim, n) array of coordinates for the nodes
                y   (n,) array of values at the nodes
                function    the radial basis function
                            'linear', 'cubic' 'thinplate', 'multiquadrics'
                            or 'gaussian', default is 'multiquadrics'
                constant    adjustable constant for gaussian or multiquadrics
                            functions - defaults to approximate average distance
                            between nodes (which is a good start)
                smooth      values greater than zero increase the smoothness
                            of the approximation. 
                            0 is for interpolation (default), the function will
                            always go through the nodal points in this case.

            Outputs: None
        """
        if len(x.shape) == 1:
            nxdim = 1
            nx = x.shape[0]
        else:
            (nxdim, nx)=x.shape
        if len(y.shape) == 1:
            nydim = 1
            ny = y.shape[0]
        else:
            (nydim, ny)=y.shape
        x.shape = (nxdim, nx)
        y.shape = (nydim, ny)
        if nx != ny:
            raise ValueError, 'x and y should have the same number of points'
        if nydim != 1:
            raise ValueError, 'y should be a length n vector'
        self.x = x
        self.y = y
        self.function = function
        if (constant==None 
            and ((function == 'multiquadrics') or (function == 'gaussian'))):
            # approx. average distance between the nodes
            constant = (s.product(x.T.max(0)-x.T.min(0),axis=0)/nx)**(1/nxdim)
        self.constant = constant
        self.smooth = smooth
        if self.function == 'linear':
            self.phi = lambda r: r
        elif self.function == 'cubic':
            self.phi = lambda r: r*r*r
        elif self.function == 'multiquadrics':
            self.phi = lambda r: s.sqrt(1.0+r*r/(self.constant*self.constant))
        elif self.function == 'thinplate':
            self.phi = lambda r: r*r*s.log(r+1)
        elif self.function == 'gaussian':
            self.phi = lambda r: s.exp(-0.5*r*r/(self.rbfconst*self.constant))
        else:
            raise ValueError, 'unkown function'
        A = self._rbf_assemble()
        b=s.r_[y.T, s.zeros((nxdim+1, 1), float)]
        self.coeff = s.linalg.solve(A,b)

    def __call__(self, xi):
        """ Evaluate the radial basis function approximation at points xi.

            Inputs:
                xi  (dim, n) array of coordinates for the points to evaluate at

            Outputs:
                y   (n,) array of values at the points xi
        """
        if len(xi.shape) == 1:
            nxidim = 1
            nxi = xi.shape[0]
        else:
            (nxidim, nxi)=xi.shape
        xi.shape = (nxidim, nxi)
        (nxdim, nx) = self.x.shape
        if nxdim != nxidim:
            raise ValueError, 'xi should have the same number of rows as an' \
                              ' array used to create RBF interpolation'
        f = s.zeros(nxi, float)
        r = s.zeros(nx, float)
        for i in range(nxi):
            st=0.0
            r = s.dot(xi[:,i,s.newaxis],s.ones((1,nx))) - self.x
            r = s.sqrt(sum(r*r))
            st = self.coeff[nx,:] + s.sum(self.coeff[0:nx,:].flatten()*self.phi(r))
            for k in range(nxdim):
                st=st+self.coeff[k+nx+1,:]*xi[k,i]
            f[i] = st
        return f

    def _rbf_assemble(self):
        (nxdim, nx)=self.x.shape
        A=s.zeros((nx,nx), float)
        for i in range(nx):
            for j in range(i+1):
                r=s.linalg.norm(self.x[:,i]-self.x[:,j])
                temp=self.phi(r)
                A[i,j]=temp
                A[j,i]=temp
            A[i,i] = A[i,i] - self.smooth
        P = s.c_[s.ones((nx,1), float), self.x.T]
        A = s.r_[s.c_[A, P], s.c_[P.T, s.zeros((nxdim+1,nxdim+1), float)]]
        return A
