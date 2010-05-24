========================================
Special functions (:mod:`scipy.special`)
========================================

.. module:: scipy.special

Nearly all of the functions below are universal functions and follow
broadcasting and automatic array-looping rules. Exceptions are noted.

Error handling
==============

Errors are handled by returning nans, or other appropriate values.
Some of the special function routines will print an error message
when an error occurs.  By default this printing
is disabled.  To enable such messages use errprint(1)
To disable such messages use errprint(0).

Example:
    >>> print scipy.special.bdtr(-1,10,0.3)
    >>> scipy.special.errprint(1)
    >>> print scipy.special.bdtr(-1,10,0.3)

.. autosummary::
   :toctree: generated/

   errprint
   errstate

Available functions
===================

Airy functions
--------------

.. autosummary::
   :toctree: generated/

   airy
   airye
   ai_zeros
   bi_zeros


Elliptic Functions and Integrals
--------------------------------

.. autosummary::
   :toctree: generated/

   ellipj
   ellipk
   ellipkinc
   ellipe
   ellipeinc

Bessel Functions
----------------

.. autosummary::
   :toctree: generated/

   jn
   jv
   jve
   yn
   yv
   yve
   kn
   kv
   kve
   iv
   ive
   hankel1
   hankel1e
   hankel2
   hankel2e

The following is not an universal function:

.. autosummary::
   :toctree: generated/

   lmbda

Zeros of Bessel Functions
^^^^^^^^^^^^^^^^^^^^^^^^^

These are not universal functions:

.. autosummary::
   :toctree: generated/

   jnjnp_zeros
   jnyn_zeros
   jn_zeros
   jnp_zeros
   yn_zeros
   ynp_zeros
   y0_zeros
   y1_zeros
   y1p_zeros

Faster versions of common Bessel Functions
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. autosummary::
   :toctree: generated/

   j0
   j1
   y0
   y1
   i0
   i0e
   i1
   i1e
   k0
   k0e
   k1
   k1e

Integrals of Bessel Functions
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. autosummary::
   :toctree: generated/

   itj0y0
   it2j0y0
   iti0k0
   it2i0k0
   besselpoly

Derivatives of Bessel Functions
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. autosummary::
   :toctree: generated/

   jvp
   yvp
   kvp
   ivp
   h1vp
   h2vp

Spherical Bessel Functions
^^^^^^^^^^^^^^^^^^^^^^^^^^

These are not universal functions:

.. autosummary::
   :toctree: generated/

   sph_jn
   sph_yn
   sph_jnyn
   sph_in
   sph_kn
   sph_inkn

Riccati-Bessel Functions
^^^^^^^^^^^^^^^^^^^^^^^^

These are not universal functions:

.. autosummary::
   :toctree: generated/

   riccati_jn
   riccati_yn

Struve Functions
----------------

.. autosummary::
   :toctree: generated/

   struve
   modstruve
   itstruve0
   it2struve0
   itmodstruve0


Raw Statistical Functions
-------------------------

.. seealso:: :mod:`scipy.stats`: Friendly versions of these functions.

.. autosummary::
   :toctree: generated/

   bdtr
   bdtrc
   bdtri
   btdtr
   btdtri
   fdtr
   fdtrc
   fdtri
   gdtr
   gdtrc
   gdtria
   gdtrib
   gdtrix
   nbdtr
   nbdtrc
   nbdtri
   pdtr
   pdtrc
   pdtri
   stdtr
   stdtridf
   stdtrit
   chdtr
   chdtrc
   chdtri
   ndtr
   ndtri
   smirnov
   smirnovi
   kolmogorov
   kolmogi
   tklmbda

Gamma and Related Functions
---------------------------

.. autosummary::
   :toctree: generated/

   gamma
   gammaln
   gammainc
   gammaincinv
   gammaincc
   gammainccinv
   beta
   betaln
   betainc
   betaincinv
   psi
   rgamma
   polygamma
   multigammaln


Error Function and Fresnel Integrals
------------------------------------

.. autosummary::
   :toctree: generated/

   erf
   erfc
   erfinv
   erfcinv
   erf_zeros
   fresnel
   fresnel_zeros
   modfresnelp
   modfresnelm

These are not universal functions:

.. autosummary::
   :toctree: generated/

   fresnelc_zeros
   fresnels_zeros

Legendre Functions
------------------

.. autosummary::
   :toctree: generated/

   lpmv
   sph_harm

These are not universal functions:

.. autosummary::
   :toctree: generated/

   lpn
   lqn
   lpmn
   lqmn

Orthogonal polynomials
----------------------

The following functions evaluate values of orthogonal polynomials:

.. autosummary::
   :toctree: generated/

   eval_legendre
   eval_chebyt
   eval_chebyu
   eval_chebyc
   eval_chebys
   eval_jacobi
   eval_laguerre
   eval_genlaguerre
   eval_hermite
   eval_hermitenorm
   eval_gegenbauer
   eval_sh_legendre
   eval_sh_chebyt
   eval_sh_chebyu
   eval_sh_jacobi

The functions below, in turn, return :ref:`orthopoly1d` objects, which
functions similarly as :ref:`numpy.poly1d`.  The :ref:`orthopoly1d`
class also has an attribute ``weights`` which returns the roots, weights,
and total weights for the appropriate form of Gaussian quadrature.
These are returned in an ``n x 3`` array with roots in the first column,
weights in the second column, and total weights in the final column.

.. autosummary::
   :toctree: generated/

   legendre
   chebyt
   chebyu
   chebyc
   chebys
   jacobi
   laguerre
   genlaguerre
   hermite
   hermitenorm
   gegenbauer
   sh_legendre
   sh_chebyt
   sh_chebyu
   sh_jacobi

.. warning::

   Large-order polynomials obtained from these functions
   are numerically unstable.

   ``orthopoly1d`` objects are converted to ``poly1d``, when doing
   arithmetic.  ``numpy.poly1d`` works in power basis and cannot
   represent high-order polynomials accurately, which can cause
   significant inaccuracy.


Hypergeometric Functions
------------------------

.. autosummary::
   :toctree: generated/

   hyp2f1
   hyp1f1
   hyperu
   hyp0f1
   hyp2f0
   hyp1f2
   hyp3f0


Parabolic Cylinder Functions
----------------------------

.. autosummary::
   :toctree: generated/

   pbdv
   pbvv
   pbwa

These are not universal functions:

.. autosummary::
   :toctree: generated/

   pbdv_seq
   pbvv_seq
   pbdn_seq

Mathieu and Related Functions
-----------------------------

.. autosummary::
   :toctree: generated/

   mathieu_a
   mathieu_b

These are not universal functions:

.. autosummary::
   :toctree: generated/

   mathieu_even_coef
   mathieu_odd_coef

The following return both function and first derivative:

.. autosummary::
   :toctree: generated/

   mathieu_cem
   mathieu_sem
   mathieu_modcem1
   mathieu_modcem2
   mathieu_modsem1
   mathieu_modsem2

Spheroidal Wave Functions
-------------------------

.. autosummary::
   :toctree: generated/

   pro_ang1
   pro_rad1
   pro_rad2
   obl_ang1
   obl_rad1
   obl_rad2
   pro_cv
   obl_cv
   pro_cv_seq
   obl_cv_seq

The following functions require pre-computed characteristic value:

.. autosummary::
   :toctree: generated/

   pro_ang1_cv
   pro_rad1_cv
   pro_rad2_cv
   obl_ang1_cv
   obl_rad1_cv
   obl_rad2_cv

Kelvin Functions
----------------

.. autosummary::
   :toctree: generated/

   kelvin
   kelvin_zeros
   ber
   bei
   berp
   beip
   ker
   kei
   kerp
   keip

These are not universal functions:

.. autosummary::
   :toctree: generated/

   ber_zeros
   bei_zeros
   berp_zeros
   beip_zeros
   ker_zeros
   kei_zeros
   kerp_zeros
   keip_zeros

Other Special Functions
-----------------------

.. autosummary::
   :toctree: generated/

   expn
   exp1
   expi
   wofz
   dawsn
   shichi
   sici
   spence
   lambertw
   zeta
   zetac

Convenience Functions
---------------------

.. autosummary::
   :toctree: generated/

   cbrt
   exp10
   exp2
   radian
   cosdg
   sindg
   tandg
   cotdg
   log1p
   expm1
   cosm1
   round
