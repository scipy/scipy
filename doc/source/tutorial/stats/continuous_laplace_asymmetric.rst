
.. _continuous-laplace:

Asymmetric Laplace Distribution
================================================================

.. math::
   :nowrap:

    \begin{eqnarray*}
    F(x, \kappa) & = & 1-\frac{\kappa^{-1}}{\kappa+\kappa^{-1}}\exp(-x\kappa),\quad x\ge0; \\
                 & = & \frac{\kappa}{\kappa+\kappa^{-1}}\exp(x/\kappa),\quad x<0. \\
    f(x, \kappa) & = & \frac{1}{\kappa+\kappa^{-1}}\exp(-x\kappa),\quad x\ge0; \\
                 & = & \frac{1}{\kappa+\kappa^{-1}}\exp(x/\kappa),\quad x<0.
    \end{eqnarray*}

.. math::
   :nowrap:

    \begin{eqnarray*}
    \mu & = & \kappa^{-1}-\kappa\\
    \sigma^2 & = & \kappa^{-2}-\kappa^2\\
    \tilde\mu_3 & = & \frac{2(1-\kappa^6)}{(1+\kappa^4)^{3/2}}
    \end{eqnarray*}


References
----------

-  "Asymmetric Laplace distribution", Wikipedia
   https://en.wikipedia.org/wiki/Asymmetric_Laplace_distribution

-  Kozubowski TJ and Podg\'orski K, "A Multivariate and Asymmetric
   Generalization of Laplace Distribution," *Computational Statistics*
   15, 531--540 (2000). :doi:`10.1007/PL00022717`


Implementation: `scipy.stats.laplace_asymmetric`
