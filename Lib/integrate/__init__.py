"""
 Methods for Integrating Functions

   gauss_quad    -- Integrate func(x) using Gaussian quadrature of order n.
   gauss_quadtol -- Integrate with given tolerance using Gaussian quadrature.
   odeint        -- Integrate ordinary differential equations.
   quad          -- General purpose integration.
   dblquad       -- General purpose double integration.
   tplquad       -- General purpose triple integration.

   See the orthogonal module (scipy.integrate.orthogonal) for Gaussian
      quadrature roots and weights.
"""

from quadrature import gauss_quad, gauss_quadtol
from odepack import odeint
from quadpack import quad, dblquad, tplquad
import orthogonal

