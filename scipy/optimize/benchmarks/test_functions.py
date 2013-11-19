import numpy as np

class SimpleQuadratic(object):
    def fun(self, x):
        return np.dot(x, x)
    
    def der(self, x):
        return 2. * x

    def hess(self, x):
        return 2. * np.eye(x.size)

class AsymmetricQuadratic(object):
    def fun(self, x):
        return np.dot(x, x) + x[0]
    
    def der(self, x):
        d = 2. * x
        d[0] += 1
        return d

    def hess(self, x):
        return 2. * np.eye(x.size)

class LJ(object):
    """
    Lennard-Jones pairwise potential energy
    """
    def __init__(self, eps=1.0, sig=1.0):
        """ simple lennard jones potential"""
        self.sig = sig
        self.eps = eps

    def vij(self, r):
        return 4.*self.eps * ( (self.sig/r)**12 - (self.sig/r)**6 )

    def dvij(self, r):
        return 4.*self.eps * ( -12./self.sig*(self.sig/r)**13 + 6./self.sig*(self.sig/r)**7 )

    def get_energy(self, coords):
        coords = np.reshape(coords, [-1,3])
        natoms = coords.shape[0]
        energy=0.
        for i in xrange(natoms):
            for j in xrange(i+1,natoms):
                dr = coords[j,:]- coords[i,:]
                r = np.linalg.norm(dr)
                energy += self.vij(r)
        return energy

    def get_energy_gradient(self, coords):
        coords = np.reshape(coords, [-1,3])
        natoms = coords.shape[0]
        energy=0.
        V = np.zeros([natoms,3])
        for i in xrange(natoms):
            for j in xrange(i+1,natoms):
                dr = coords[j,:]- coords[i,:]
                r = np.linalg.norm(dr)
                energy += self.vij(r)
                g = self.dvij(r)
                V[i,:] += -g * dr/r
                V[j,:] += g * dr/r
        V = V.reshape([natoms*3])
        return energy,V

    def get_gradient(self, coords):
        e, g = self.get_energy_gradient(coords)
        return g

class Booth(object):
#    target_E = 0.
#    target_coords = np.array([1., 3.])
#    xmin = np.array([-10., -10.])
##    xmin = np.array([0., 0.])
#    xmax = np.array([10., 10.])
    def fun(self, coords):
        x, y = coords
        return (x + 2.*y - 7.)**2 + (2.*x + y - 5.)**2

    def der(self, coords):
        x, y = coords
        dx = 2.*(x + 2.*y - 7.) + 4.*(2.*x + y - 5.)
        dy = 4.*(x + 2.*y - 7.) + 2.*(2.*x + y - 5.)
        return np.array([dx, dy])

