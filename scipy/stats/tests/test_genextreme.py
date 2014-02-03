from __future__ import division, print_function, absolute_import

import numpy as np
from numpy import inf
from numpy.testing import assert_allclose, assert_equal, run_module_suite
from scipy.stats.distributions import genextreme, frechet, weibull_max, gumbel_r

"""
Test both the generalized extreme value distribution (genextreme) and related
special cases weibull_max (c<0), gumbel_r (c=0), and frechet (c>0).
"""

def test_genextreme_support(mu=1,sigma=2):
    """The support of this distribution depends on c,mu,sigma."""

    # Gumbel case; unbounded support
    c=0
    gev=genextreme(c,loc=mu,scale=sigma)
    assert_allclose(gev.interval(1),(-inf,inf))

    gum=gumbel_r(loc=mu,scale=sigma)
    assert_allclose(gum.interval(1),(-inf,inf))
    assert_allclose(gum.cdf([-1,-.1,0,.1,1]),gev.cdf([-1,-.1,0,.1,1]))

    # Frechet case
    c=1
    gev=genextreme(c,loc=mu,scale=sigma)
    assert_allclose(gev.interval(1),(mu-sigma/c,inf))

    F=frechet(1./abs(c),loc=(sigma/abs(c)+mu),scale=sigma/abs(c))
    assert_allclose(F.interval(1),mu-sigma/c,inf)

    # Weibull case
    c=-1
    gev=genextreme(c,loc=mu,scale=sigma)
    assert_allclose(gev.interval(1),(-inf,mu-sigma/c))

    W=weibull_max(1./abs(c),loc=(sigma/abs(c)+mu),scale=sigma/abs(c))
    assert_allclose(W.interval(1),-inf,mu-sigma/c)


def test_genextreme_known_exact():
    """
    Compare results with some known exact formulas.

    Some exact values of the generalized extreme value cdf with (mu,sigma)=(0,1)

    == == =============
     c  x cdf
    == == =============
    -1 -1 exp(-2)
    -1  0 exp(-1)
     0 -1 exp(-exp(1))
     0  0 exp(-1)
     0  1 exp(-exp(-1))
     1  0 exp(-1)
     1  1 exp(-(1/2))
    == == =============
    """
    # c = -1
    gev=genextreme(-1)
    assert_allclose(gev.cdf([-1,0]), np.exp([-2,-1]))

    # c = 0
    gev=genextreme(0)
    assert_allclose(gev.cdf([-1,0,1]), np.exp([-np.exp(1),-1,-np.exp(-1)]))

    # c = 1
    gev=genextreme(1)
    assert_allclose(gev.cdf([0,1]), np.exp([-1,-1./2]))

def test_genextreme_equivalents(mu=1,sigma=2):
    """
    Test equivalency between GEV and Frechet and Weibull distributions.
    """
    c=-1 # EVI parameter
    G=genextreme(c,loc=mu,scale=sigma).cdf
    W=weibull_max(-1./c,loc=(mu-sigma/c),scale=sigma/abs(c)).cdf

    assert_allclose(G(0),W(0))

    assert_allclose(G(-1),W(-1))

    c=1 # EVI parameter
    G=genextreme(c,loc=mu,scale=sigma).cdf
    F=frechet(1./abs(c),loc=mu-sigma/c,scale=sigma/abs(c)).cdf

    assert_allclose(G(0),F(0))

    assert_allclose(G(1),F(1))




if __name__ == "__main__":
    run_module_suite()
