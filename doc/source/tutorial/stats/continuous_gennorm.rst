
.. _continuous-gennorm:

Generalized Normal Distribution
===============================

This distribution is also known as the exponential power distribution. It has
a single shape parameter :math:`\beta>0`. It reduces to a number of common
distributions.

Functions
---------

.. math::
    :nowrap:

    \begin{eqnarray*} f\left(x; \beta\right) & = &\frac{\beta}{2\Gamma(1/\beta)} e^{-\left|x\right|^{\beta}} \end{eqnarray*}

.. math::
    :nowrap:

    \begin{eqnarray*} F\left(x; \beta\right) & = & \frac{1}{2} + \mathrm{sgn}\left(x\right) \frac{\gamma\left(1/\beta, x^{\beta}\right)}{2\Gamma\left(1/\beta\right)} \end{eqnarray*}
     
:math:`\gamma` is the lower incomplete gamma function. 
:math:`\gamma\left(s, x\right) = \int_0^x t^{s-1} e^{-t} dt`.

.. math::
    :nowrap:

    \begin{eqnarray*} h\left[X; \beta\right] = \frac{1}{\beta} - \log\left(\frac{\beta}{2\Gamma\left(1/\beta\right)}\right)\end{eqnarray*}

Moments
-------

.. math::
   :nowrap:

    \begin{eqnarray*}
    \mu & = & 0 \\
    m_{n} & = & 0 \\
    m_{d} & = & 0 \\
    \mu_2 &  = & \frac{\Gamma\left(3/\beta\right)}{\gamma\left(1/\beta\right)} \\
    \gamma_1 & = & 0 \\
    \gamma_2 & = & \frac{\Gamma\left(5/\beta\right) \Gamma\left(1/\beta\right)}{\Gamma\left(3/\beta\right)^2} - 3 \\
    \end{eqnarray*}


Special Cases
-------------
* Laplace distribution (:math:`\beta = 1`)
* Normal distribution with :math:`\mu_2 = 1/2` (:math:`\beta = 2`)
* Uniform distribution over the interval :math:`[-1, 1]`
  (:math:`\beta \rightarrow \infty`)

Sources
-------
* https://en.wikipedia.org/wiki/Generalized_normal_distribution#Version_1
* https://en.wikipedia.org/wiki/Incomplete_gamma_function#Lower_incomplete_Gamma_function

Implementation: `scipy.stats.gennorm`
