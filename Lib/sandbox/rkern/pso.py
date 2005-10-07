"""Particle Swarm Optimization

:Author: Robert Kern

Copyright 2005 by Robert Kern.
"""

import scipy as sp

def pso(func, nswarm, lbound, ubound, vmax, args=(), maxiter=1000, cp=2.0, cg=2.0):
    ndim = len(lbound)
    lbound = sp.asarray(lbound)
    ubound = sp.asarray(ubound)
    vmax = sp.asarray(vmax)

    # initialize the swarm
    swarm = lbound + sp.rand(nswarm, ndim)*(ubound-lbound)

    # initialize the "personal best" values
    pbestv = sp.zeros(nswarm, sp.Float)
    for i in sp.arange(nswarm):
        pbestv[i] = func(swarm[i])
    pbest = sp.array(swarm)

    # initialize the "global best" values
    gbesti = sp.argmin(pbestv)
    gbestv = pbestv[gbesti]
    gbest = pbest[gbesti]

    # initialize velocities
    velocities = 2*vmax*sp.randn(nswarm, ndim) - vmax

    for i in sp.arange(maxiter):

        values = sp.zeros(nswarm, sp.Float)
        for j in sp.arange(nswarm):
            values[j] = func(swarm[j])

        mask = values < pbestv
        mask2d = sp.repeat(mask, ndim)
        mask2d.shape = (nswarm, ndim)
        pbestv = sp.where(mask, values, pbestv)
        pbest = sp.where(mask2d, swarm, pbest)

        if sp.minimum.reduce(pbestv) < gbestv:
            gbesti = sp.argmin(pbestv)
            gbestv = pbestv[gbesti]
            gbest = pbest[gbesti]

        velocities += (cp*sp.rand()*(pbest - swarm) + 
                       cg*sp.rand()*(gbest - swarm))
        velocities = sp.clip(velocities, -vmax, vmax)
        swarm += velocities
        swarm = sp.clip(swarm, lbound, ubound)
        yield gbest

##    return gbest

def _testfunc(x):
    return -sp.multiply.reduce(sp.cos(x/2.0/sp.pi))*sp.exp(-sp.sum(x*x))


