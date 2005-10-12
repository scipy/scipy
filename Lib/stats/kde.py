#-------------------------------------------------------------------------------
#
#  Define classes for (uni/multi)-variate kernel density estimation. 
#
#  Currently, only Gaussian kernels are implemented.
#
#  Written by: Robert Kern
#
#  Date: 2004-08-09
#
#  Modified: 2005-02-10 by Robert Kern.
#            Contributed to Scipy.
#
#  Copyright 2004-2005 by Enthought, Inc.
#
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
# Imports:
#-------------------------------------------------------------------------------

import warnings

from scipy import dot, linalg, special
from scipy_base import *

from distributions import randint
import stats
from rv import multivariate_normal
import mvn

#-------------------------------------------------------------------------------
# 'gaussian_kde' class:
#-------------------------------------------------------------------------------

class gaussian_kde(object):
    """
    Parameters
    ----------
    dataset : (# of dims, # of data)-array
      datapoints to estimate from

    Members
    -------
    d : int
      number of dimensions
    n : int
      number of datapoints

    Methods
    -------
    kde.evaluate(points) : array
      evaluate the estimated pdf on a provided set of points
    kde(points) : array
      same as kde.evaluate(points)
    kde.integrate_gaussian(mean, cov) : float
      multiply pdf with a specified gaussian and integrate over the whole domain
    kde.integrate_box_1d(low, high) : float
      integrate pdf (1D only) between two bounds
    kde.integrate_box(low_bounds, high_bounds) : float
      integrate pdf over a rectangular space between low_bounds and high_bounds
    kde.integrate_kde(other_kde) : float
      integrate two kernel density estimates multiplied together

   Internal Methods
   ----------------
    kde.covariance_factor() : float
      computes the coefficient that multiplies the data covariance matrix to
      obtain the kernel covariance matrix. Set this method to kde.scotts_factor
      or kde.silverman_factor (or subclass to provide your own). The default is
      scotts_factor.
    """

    def __init__(self, dataset):
        self.dataset = atleast_2d(dataset)

        self.d, self.n = self.dataset.shape

        self._compute_covariance()
   
    def evaluate(self, points):
        """Evaluate the estimated pdf on a set of points.

        points.shape == (# of dimensions, # of points)
        """

        points = atleast_2d(points).astype(self.dataset.typecode())
        
        d, m = points.shape
        if d != self.d:
            if d == 1 and m == self.d:
                # points was passed in as a row vector
                points = reshape(points, (self.d, 1))
            else:
                msg = "points have dimension %s, dataset has dimension %s" % (d,
                    self.d)
                raise ValueError, msg

        result = zeros((m,), points.typecode())
    
        try:
            inv_cov = linalg.inv(self.covariance)
        except linalg.LinAlgError:
            warnings.warn('dataset is degenerate; cannot evaluate; returning zeros')
            return zeros(m, Float)

        if m >= self.n:
            # there are more points than data, so loop over data
            for i in range(self.n):
                diff = self.dataset[:,i,NewAxis] - points
                tdiff = dot(inv_cov, diff)
                energy = sum(diff*tdiff)/2.0
                result += exp(-energy)
        else:
            # loop over points
            for i in range(m):
                diff = self.dataset - points[:,i,NewAxis]
                tdiff = dot(self.inv_cov, diff)
                energy = sum(diff*tdiff)/2.0
                result[i] = sum(exp(-energy))

        det_cov = linalg.det(2*pi*self.covariance)
        result /= sqrt(det_cov)*self.n

        return result

    __call__ = evaluate

    def integrate_gaussian(self, mean, cov):
        """Multiply estimated density by a multivariate gaussian and integrate
        over the wholespace.

        Parameters
        ----------
        mean : vector
          the mean of the gaussian
        cov : matrix
          the covariance matrix of the gaussian

        Returns
        -------
        result : scalar
          the value of the integral

        """

        mean = atleast_1d(squeeze(mean))
        cov = atleast_2d(cov)
        
        if mean.shape != (self.d,):
            raise ValueError, "mean does not have dimension %s" % self.d
        if cov.shape != (self.d, self.d):
            raise ValueError, "covariance does not have dimension %s" % self.d

        # make mean a column vector
        mean = mean[:,NewAxis]

        sum_cov = self.covariance + cov
        sum_cov_inv = linalg.inv(sum_cov)

        diff = self.dataset - mean
        tdiff = dot(sum_cov_inv, diff)

        energies = sum(diff*tdiff)/2.0
        result = sum(exp(-energies))/sqrt(linalg.det(2*pi*sum_cov))/self.n

        return result

    def integrate_box_1d(self, low, high):
        """Computes the integral of a 1D pdf between two bounds.

        Parameters
        ----------
        low : scalar
          lower bound of integration
        high : scalar
          upper bound of integration

        Returns
        -------
        value : scalar
          the result of the integral
        """
        if self.d != 1:
            raise ValueError("integrate_box_1d() only handles 1D pdfs")

        stdev = ravel(sqrt(self.covariance))[0]

        normalized_low = ravel((low - self.dataset)/stdev)
        normalized_high = ravel((high - self.dataset)/stdev)

        value = stats.mean(special.ndtr(normalized_high) -
                           special.ndtr(normalized_low))
        return value


    def integrate_box(self, low_bounds, high_bounds, maxpts=None):
        """Computes the integral of a pdf over a rectangular interval.

        Parameters
        ----------
        low_bounds : vector
          lower bounds of integration
        high_bounds : vector
          upper bounds of integration
        maxpts=None : int
          maximum number of points to use for integration

        Returns
        -------
        value : scalar
          the result of the integral
        """
        if maxpts is not None:
            extra_kwds = {'maxpts': maxpts}
        else:
            extra_kwds = {}

        value, inform = mvn.mvnun(low_bounds, high_bounds, self.dataset,
            self.covariance, **extra_kwds)
        if inform:
            msg = ('an integral in mvn.mvnun requires more points than %s' %
                (self.d*1000))
            warnings.warn(msg)

        return value

    def integrate_kde(self, other):
        """Computes the integral of the product of this  kernel density estimate
        with another.

        Parameters
        ----------
        other : GaussianKDE instance
          the other kde

        Returns
        -------
        value : scalar
          the result of the integral
        """

        if other.d != self.d:
            raise ValueError("KDEs are not the same dimensionality")

        # we want to iterate over the smallest number of points
        if other.n < self.n:
            small = other
            large = self
        else:
            small = self
            large = other

        sum_cov = small.covariance + large.covariance
        result = 0.0
        for i in range(small.n):
            mean = small.dataset[:,i,NewAxis]
            diff = large.dataset - mean
            tdiff = dot(linalg.inv(sum_cov), diff)

            energies = sum(diff*tdiff)/2.0
            result += sum(exp(-energies))
        
        result /= sqrt(linalg.det(2*pi*sum_cov))*large.n*small.n

        return result

    def resample(self, size=None):
        """Randomly sample a dataset from the estimated pdf.

        Parameters
        ----------
        size : int=None
          if None, then size = self.n; otherwise, the number of samples to draw
        
        Returns
        -------
        dataset : (self.d, size)-array
          sampled dataset
        """

        if size is None:
            size = self.n

        norm = transpose(multivariate_normal(zeros((self.d,), Float),
            self.covariance, size=size))
        indices = randint.rvs(0, self.n, size=size)
        means = take(self.dataset, indices, axis=1)

        return means + norm


    def scotts_factor(self):
        return power(self.n, -1./(self.d+4))

    def silverman_factor(self):
        return power(self.n*(self.d+2.0)/4.0, -1./(self.d+4))

    covariance_factor = scotts_factor

    def _compute_covariance(self):
        """Computes the covariance matrix for each gaussian kernel using
        covariance_factor
        """
        self.factor = self.covariance_factor()
        self.covariance = atleast_2d(stats.cov(self.dataset, rowvar=1) * 
            self.factor * self.factor)
