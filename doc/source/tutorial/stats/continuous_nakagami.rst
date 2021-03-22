
.. _continuous-nakagami:

Nakagami Distribution
=====================

Generalization of the chi distribution. Shape parameter is :math:`\nu>0.` The support is :math:`x\geq0.`

.. math::
   :nowrap:

    \begin{eqnarray*}
    f\left(x;\nu\right) & = & \frac{2\nu^{\nu}}{\Gamma\left(\nu\right)}x^{2\nu-1}\exp\left(-\nu x^{2}\right)\\
    F\left(x;\nu\right) & = & \frac{\gamma\left(\nu,\nu x^{2}\right)}{\Gamma\left(\nu\right)}\\
    G\left(q;\nu\right) & = & \sqrt{\frac{1}{\nu}\gamma^{-1}\left(\nu,q{\Gamma\left(\nu\right)}\right)}
    \end{eqnarray*}

where :math:`\gamma` is the lower incomplete gamma function, :math:`\gamma\left(\nu, x\right) = \int_0^x t^{\nu-1} e^{-t} dt`.

.. math::
   :nowrap:

    \begin{eqnarray*}
    \mu & = & \frac{\Gamma\left(\nu+\frac{1}{2}\right)}{\sqrt{\nu}\Gamma\left(\nu\right)}\\
    \mu_{2} & = & \left[1-\mu^{2}\right]\\
    \gamma_{1} & = & \frac{\mu\left(1-4v\mu_{2}\right)}{2\nu\mu_{2}^{3/2}}\\
    \gamma_{2} & = & \frac{-6\mu^{4}\nu+\left(8\nu-2\right)\mu^{2}-2\nu+1}{\nu\mu_{2}^{2}}
    \end{eqnarray*}

Implementation: `scipy.stats.nakagami`


MLE of the Nakagami Distribution in SciPy (:code:`nakagami.fit`)
----------------------------------------------------------------

The probability density function of the :code:`nakagami` distribution in SciPy is

.. math::
   :nowrap:

    \begin{equation}
    f(x; \nu, \mu, \sigma) = 2 \frac{\nu^\nu}{ \sigma \Gamma(\nu)}\left(\frac{x-\mu}{\sigma}\right)^{2\nu - 1} \exp\left(-\nu \left(\frac{x-\mu}{\sigma}\right)^2 \right),\tag{1}
    \end{equation}

for :math:`x` such that :math:`\frac{x-\mu}{\sigma} \geq 0`, where :math:`\nu \geq \frac{1}{2}` is the shape parameter,
:math:`\mu` is the location, and :math:`\sigma` is the scale.

The log-likelihood function is therefore

.. math::
   :nowrap:

    \begin{equation}
    l(\nu, \mu, \sigma) = \sum_{i = 1}^{N} \log \left( 2 \frac{\nu^\nu}{ \sigma\Gamma(\nu)}\left(\frac{x_i-\mu}{\sigma}\right)^{2\nu - 1} \exp\left(-\nu \left(\frac{x_i-\mu}{\sigma}\right)^2 \right) \right),\tag{2}
    \end{equation}

which can be expanded as

.. math::
   :nowrap:

    \begin{equation}
    l(\nu, \mu, \sigma) = N \log(2) + N\nu \log(\nu) - N\log\left(\Gamma(\nu)\right) - 2N \nu \log(\sigma) + \left(2 \nu - 1 \right) \sum_{i=1}^N  \log(x_i - \mu) - \nu \sigma^{-2} \sum_{i=1}^N \left(x_i-\mu\right)^2, \tag{3}
    \end{equation}

Leaving supports constraints out, the first-order condition for optimality on the likelihood derivatives gives estimates of parameters:

.. math::
   :nowrap:

    \begin{align}
    \frac{\partial l}{\partial \nu}(\nu, \mu, \sigma) &= N\left(1 + \log(\nu) - \psi^{(0)}(\nu)\right) + 2 \sum_{i=1}^N \log \left( \frac{x_i - \mu}{\sigma} \right) - \sum_{i=1}^N \left( \frac{x_i - \mu}{\sigma} \right)^2  = 0
    \text{,} \tag{4}\\
    \frac{\partial l}{\partial \mu}(\nu, \mu, \sigma) &= (1 - 2 \nu) \sum_{i=1}^N \frac{1}{x_i-\mu} + \frac{2\nu}{\sigma^2} \sum_{i=1}^N x_i-\mu = 0
    \text{, and} \tag{5}\\
    \frac{\partial l}{\partial \sigma}(\nu, \mu, \sigma) &= -2 N \nu \frac{1}{\sigma} + 2 \nu \sigma^{-3} \sum_{i=1}^N \left(x_i-\mu\right)^2 = 0
    \text{,}\tag{6}
    \end{align}

where :math:`\psi^{(0)}` is the polygamma function of order :math:`0`; i.e. :math:`\psi^{(0)}(\nu) = \frac{d}{d\nu} \log \Gamma(\nu)`.

However, the support of the distribution is the values of :math:`x` for which :math:`\frac{x-\mu}{\sigma} \geq 0`, and this provides an additional constraint that

.. math::
   :nowrap:

    \begin{equation}
    \mu \leq \min_i{x_i}.\tag{7}
    \end{equation}

For :math:`\nu = \frac{1}{2}`, the partial derivative of the log-likelihood with respect to :math:`\mu` reduces to:

.. math::
   :nowrap:

    \begin{equation}
    \frac{\partial l}{\partial \mu}(\nu, \mu, \sigma) = {\sigma^2} \sum_{i=1}^N (x_i-\mu),
    \end{equation}

which is positive when the support constraint is satisfied. Because the partial derivative with respect to :math:`\mu`
is positive, increasing :math:`\mu` increases the log-likelihood, and therefore the constraint is active at the maximum likelihood estimate for :math:`\mu`

.. math::
   :nowrap:

    \begin{equation}
    \mu = \min_i{x_i}, \quad \nu = \frac{1}{2}. \tag{8}
    \end{equation}

For :math:`\nu` sufficiently greater than :math:`\frac{1}{2}`, the likelihood equation :math:`\frac{\partial l}{\partial \mu}(\nu, \mu, \sigma)=0` has a solution, and this solution provides the maximum likelihood estimate for :math:`\mu`. In either case, however, the condition :math:`\mu = \min_i{x_i}` provides a reasonable initial guess for numerical optimization.

Furthermore, the likelihood equation for :math:`\sigma` can be solved explicitly, and it provides the maximum likelihood estimate

.. math::
   :nowrap:

    \begin{equation}
    \sigma = \sqrt{ \frac{\sum_{i=1}^N \left(x_i-\mu\right)^2}{N}}. \tag{9}
    \end{equation}

Hence, the :code:`_fitstart` method for :code:`nakagami` uses

.. math::
   :nowrap:

    \begin{align}
    \mu_0 &= \min_i{x_i} \,
    \text{and} \\
    \sigma_0 &= \sqrt{ \frac{\sum_{i=1}^N \left(x_i-\mu_0\right)^2}{N}}
    \end{align}

as initial guesses for numerical optimization accordingly.