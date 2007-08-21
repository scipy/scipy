"""
Statistical models

 - model `formula`
 - standard `regression` models

  - `ols_model` (ordinary least square regression)
  - `wls_model` (weighted least square regression)
  - `ar_model` (autoregression)

 - `glm.model` (generalized linear models)
 - robust statistical models

  - `rlm.model` (robust linear models using M estimators)
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
