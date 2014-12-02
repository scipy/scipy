# -*- coding: utf-8 -*-
from __future__ import division, print_function, absolute_import
"""
==============================================================================
`go_benchmark_functions` --  Problems for testing global optimization routines
==============================================================================

This module provides a comprehensive set of problems for benchmarking global
optimization routines, such as scipy.optimize.basinhopping, or
scipy.optimize.differential_evolution.  The purpose is to see whether a given
optimization routine can find the global minimum, and how many function
evaluations it requires to do so.
The range of problems is extensive, with a range of difficulty. The problems are
multivariate, with N=2 to N=17 provided.

References
----------
.. [1] Momin Jamil and Xin-She Yang, A literature survey of benchmark
    functions for global optimization problems, Int. Journal of Mathematical
    Modelling and Numerical Optimisation, Vol. 4, No. 2, pp. 150--194 (2013).
    http://arxiv.org/pdf/1308.4008v1.pdf
    (and references contained within)
.. [2] http://infinity77.net/global_optimization/index.html
.. [3] S. K. Mishra, Global Optimization By Differential Evolution and
    Particle Swarm Methods: Evaluation On Some Benchmark Functions, Munich
    Research Papers in Economics
.. [4] E. P. Adorio, U. P. Dilman, MVF - Multivariate Test Function Library
    in C for Unconstrained Global Optimization Methods, [Available Online]:
    http://www.geocities.ws/eadorio/mvf.pdf
.. [5] S. K. Mishra, Some New Test Functions For Global Optimization And
    Performance of Repulsive Particle Swarm Method, [Available Online]:
    http://mpra.ub.uni-muenchen.de/2718/
.. [6] NIST StRD Nonlinear Regression Problems, retrieved on 1 Oct, 2014
    http://www.itl.nist.gov/div898/strd/nls/nls_main.shtml

"""

"""
Copyright 2013 Andrea Gavana
Author: <andrea.gavana@gmail.com>

Modifications 2014 Andrew Nelson
<andyfaff@gmail.com>
"""

from .go_funcs_A import *
from .go_funcs_B import *
from .go_funcs_C import *
from .go_funcs_D import *
from .go_funcs_E import *
from .go_funcs_F import *
from .go_funcs_G import *
from .go_funcs_H import *
from .go_funcs_I import *
from .go_funcs_J import *
from .go_funcs_K import *
from .go_funcs_L import *
from .go_funcs_M import *
from .go_funcs_N import *
from .go_funcs_O import *
from .go_funcs_P import *
from .go_funcs_Q import *
from .go_funcs_R import *
from .go_funcs_S import *
from .go_funcs_T import *
from .go_funcs_U import *
from .go_funcs_V import *
from .go_funcs_W import *
from .go_funcs_X import *
from .go_funcs_Y import *
from .go_funcs_Z import *

__all__ = [s for s in dir() if not s.startswith('_')]
