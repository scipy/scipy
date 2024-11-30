
.. _continuous-dpareto_lognorm:

Double Pareto Lognormal Distribution
====================================

For real numbers :math:`x` and :math:`\mu`, :math:`\sigma > 0`,
:math:`\alpha > 0`, and :math:`\beta > 0`, the PDF of a double
Pareto lognormal distribution is:

.. math::
   :nowrap:

    \begin{eqnarray*}
        f(x, \mu, \sigma, \alpha, \beta) =
        \frac{\alpha \beta}{(\alpha + \beta) x}
        \phi\left( \frac{\log x - \mu}{\sigma} \right)
        \left( R(y_1) + R(y_2) \right)
    \end{eqnarray*}

where :math:`R(t) = \frac{1 - \Phi(t)}{\phi(t)}` is a Mills' ratio,
:math:`y_1 = \alpha \sigma - \frac{\log x - \mu}{\sigma}`,
and :math:`y_2 = \beta \sigma + \frac{\log x - \mu}{\sigma}`.
The CDF is:

.. math::
   :nowrap:

    \begin{eqnarray*}
        F(x, \mu, \sigma, \alpha, \beta) =
        \Phi \left(\frac{\log x - \mu}{\sigma} \right) -
        \phi \left(\frac{\log x - \mu}{\sigma} \right)
        \left(\frac{\beta R(x_1) - \alpha R(x_2)}{\alpha + \beta} \right)
    \end{eqnarray*}

Raw moment :math:`k > \alpha` is given by:

.. math::
   :nowrap:

    \begin{eqnarray*}
        \mu_k' = \frac{\alpha \beta}{(\alpha - k)(\beta + k)} 
                 \exp \left(k \mu + \frac{k^2 \sigma^2}{2} \right)
    \end{eqnarray*}

Implementation: `scipy.stats.dpareto_lognorm`
