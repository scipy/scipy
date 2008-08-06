"""
Statistical models

 - model `formula`
 - standard `regression` models

  - `OLSModel` (ordinary least square regression)
  - `WLSModel` (weighted least square regression)
  - `ARModel` (autoregressive model)

 - `glm.Model` (generalized linear models)
 - robust statistical models

  - `rlm.Model` (robust linear models using M estimators)
  - `robust.norms` estimates
  - `robust.scale` estimates (MAD, Huber's proposal 2).

 - `mixed` effects models
 - `gam` (generalized additive models)
"""
__docformat__ = 'restructuredtext en'

depends = ['weave',
           'special.orthogonal',
           'integrate',
           'optimize',
           'linalg']

postpone_import = True
