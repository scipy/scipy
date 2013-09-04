""" Test Hessian estimation using finite differences.

"""
from __future__ import division, print_function, absolute_import

import numpy as np
from scipy.diff import Hessian
from scipy import linalg, optimize, constants 

from numpy.testing import run_module_suite, TestCase


_TINY = np.finfo(float).machar.tiny

#-----------------------------------------------------------------------------------------
#       Hamiltonian      H=sum_i(p_i2/(2m)+ 1/2 * m * w2 x_i2)+ sum_(i!=j)(a/|x_i-x_j|)
#-----------------------------------------------------------------------------------------

class classicalHamiltonian:
    def __init__(self):
       
        self.N = 2                            #N is a scalar, it's the number of ions in the chain
        f = 1000000                            #f is a scalar, it's the trap frequency
        self.w = 2 * np.pi * f                         #w is a scalar, it's the angular velocity corresponding to the trap frequency
        self.C = (4 * np.pi * constants.epsilon_0) ** (-1) * constants.e ** 2    #C is a scalar, it's the Coulomb constant times the electronic charge in SI
        self.m = 39.96 * 1.66e-27                        #m is the mass of a single trapped ion in the chain
       
    def potential(self, positionvector):                     #Defines the potential that is going to be minimized
       
        x = positionvector                         #x is an 1-d array (vector) of length N that contains the positions of the N ions
        w = self.w
        C = self.C
        m = self.m
       
       
        #First we consider the potential of the harmonic oscillator
        Vx = 0.5 * m * (w ** 2) * sum(x ** 2)               
       
        for i in range(len(x)):
            for j in range(i + 1, len(x)):
                Vx += C / (abs(x[i] - x[j]))    #then we add the coulomb interaction
               
        return Vx
    
    def initialposition(self):        #Defines the initial position as an estimate for the minimize process
       
        N = self.N
        x_0 = np.r_[-(N - 1) / 2:(N - 1) / 2:N * 1j]
        return x_0
       
    def normal_modes(self, eigenvalues):    #the computed eigenvalues of the matrix Vx are of the form (normal_modes)2*m.
        m = self.m
        normal_modes = np.sqrt(eigenvalues / m)
        return normal_modes

def _run_hamiltonian(verbose=True):
    #C=(4*np.pi*constants.epsilon_0)**(-1)*constants.e**2
    c = classicalHamiltonian()
    if verbose:
        print(c.potential(np.array([-0.5, 0.5])))
        print(c.potential(np.array([-0.5, 0.0])))
        print(c.potential(np.array([0.0, 0.0])))
        
    xopt = optimize.fmin(c.potential, c.initialposition(), xtol=1e-10, disp=False)
    # Important to restrict the step in order to avoid the discontinutiy at x=[0,0]
    hessian = Hessian(c.potential, step_max=1.0, stepNom=np.abs(xopt))
    #hessian = Hessian(c.potential)
    #hessian = algopy.Hessian(c.potential) # Does not work
    #hessian = scientific.Hessian(c.potential) # does not work
    H = hessian(xopt)
    true_H = np.array([[  5.23748385e-12,  -2.61873829e-12],
                       [ -2.61873829e-12,   5.23748385e-12]])
    if verbose:
        print(xopt)
        print('H', H)
        print('H-true_H', np.abs(H-true_H))
        print('error_estimate', hessian.error_estimate)
        
        eigenvalues = linalg.eigvals(H)
        normal_modes = c.normal_modes(eigenvalues)
        
        print('eigenvalues', eigenvalues)
        print('normal_modes', normal_modes)
    return H, hessian.error_estimate, true_H
    
class TestHessian(TestCase):
    def test_hessian(self):
        H, error_estimate, true_H =_run_hamiltonian(verbose=False) #@UnusedVariable
        self.assertTrue((np.abs(H-true_H)<1e-18).all())
       

if __name__ == "__main__":
    run_module_suite()
