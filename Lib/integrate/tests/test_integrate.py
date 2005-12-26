#!/usr/bin/env python

# Test provided by Nils Wagner.
# File created by Ed Schofield on Nov 16. 

""" Tests for numerical integration.
"""

import scipy.base
from scipy.base import arange, zeros, array, dot, sqrt, cos, sin
from scipy.linalg import norm
from scipy.testing import *
set_package_path()
from scipy.integrate import odeint
restore_path()

class test_odeint(ScipyTestCase):
    """ Test odeint: free vibration of a simple oscillator
        m \ddot{u} + k u = 0, u(0) = u_0 \dot{u}(0) \dot{u}_0

    Solution:
        u(t) = u_0*cos(sqrt(k/m)*t)+\dot{u}_0*sin(sqrt(k/m)*t)/sqrt(k/m)
    """

    def setUp(self):
        self.k = 4.0
        self.m = 1.0

    def F(self, z, t):
        tmp = zeros((2,2), float)
        tmp[0,1] = 1.0
        tmp[1,0] = -self.k / self.m 
        return dot(tmp,z)

    def check_odeint1(self):
        omega = sqrt(self.k / self.m)
        z0 = zeros(2, float)
        z0[0] = 1.0     # initial displacement
        z0[1] = 0.1     # initial velocity
        t = arange(0.0, 1+0.09, 0.1)

        # Analytical solution
        #
        u = z0[0]*cos(omega*t)+z0[1]*sin(omega*t)/omega

        # Numerical solution
        z, infodict = odeint(self.F, z0, t, full_output=True)

        res = norm(u - z[:,0])
        print 'Residual:', res
        assert res < 1.0e-6

if __name__ == "__main__":
    ScipyTest().run()
