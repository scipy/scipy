"""
Statistical models
==================

This module contains a several linear statistical models
- model formulae as in R (to some extent)
- OLS (ordinary least square regression)
- WLS (weighted least square regression)
- ARp regression
- GLMS (generalized linear models)
- robust linear models using M estimators (with a number of standard default robust norms as in R's rlm)
- robust scale estimates (MAD, Huber's proposal 2).
- mixed effects models
- generalized additive models (gam)
"""

depends = ['weave',
           'special.orthogonal',
           'integrate',
           'optimize',
           'linalg']

postpone_import = True
