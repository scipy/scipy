"""
 Methods for Integrating Functions given function object.

   odeint        -- Integrate ordinary differential equations.
   quad          -- General purpose integration.
   dblquad       -- General purpose double integration.
   tplquad       -- General purpose triple integration.
   fixed_quad    -- Integrate func(x) using Gaussian quadrature of order n.
   quadrature    -- Integrate with given tolerance using Gaussian quadrature.
   romberg       -- Integrate func using Romberg integration.

 Methods for Integrating Functions given fixed samples.

   trapz         -- Use trapezoidal rule to compute integral of sampled function
   simps         -- Use Simpson's rule to compute integral from samples.
   romb          -- Use Romberg Integration to compute integral from
                    (2**k + 1) evenly-spaced samples.

   See the orthogonal module (scipy.integrate.orthogonal) for Gaussian
      quadrature roots and weights for other weighting factors and regions.

 Interface to numerical integrators of ODE systems.

   ode           -- Integrate ODE using vode routine.

"""
_moddict = {'quadrature': ['fixed_quad', 'quadrature','romberg', 'trapz',
                           'simps', 'romb'],
            'odepack': ['odeint'],
            'quadpack': ['quad', 'dblquad', 'tplquad', 'quad_explain', 'Inf',
                         'inf'],
            'ode':['ode'],
            }

__all__ = []
import scipy
scipy.somenames2all(__all__, _moddict, globals())
scipy.modules2all(__all__, ['orthogonal'], globals())
del scipy
