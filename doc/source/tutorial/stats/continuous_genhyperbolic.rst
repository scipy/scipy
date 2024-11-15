
.. _continuous-genhyperbolic:

Generalized Hyperbolic Distribution
===================================

The Generalized Hyperbolic Distribution is defined as the normal variance-mean mixture with Generalized Inverse Gaussian distribution as the mixing distribution.
The "hyperbolic" characterization refers to the fact that the shape of the log-probability distribution can be described as a hyperbola. Hyperbolic distributions are sometime referred to as semi-fat tailed because their probability density decrease slower than "sub-hyperbolic" distributions (e.g. normal distribution, whose log-probability decreases quadratically), but faster than other "extreme value" distributions (e.g. `pareto` distribution, whose log-probability decreases logarithmically).

Functions
---------

Different parameterizations exist in the literature; SciPy implements the "4th parametrization" in Prause (1999).

.. math::
   :nowrap:

    \begin{eqnarray*}
        f(x, p, a, b) & = &
        \frac{(a^2 - b^2)^{p/2}}
        {\sqrt{2\pi}a^{p-0.5}
        K_p\Big(\sqrt{a^2 - b^2}\Big)}
        e^{bx} \times \frac{K_{p - 1/2}
        (a \sqrt{1 + x^2})}
        {(\sqrt{1 + x^2})^{1/2 - p}}
    \end{eqnarray*}

for:

-  :math:`x, p \in ( - \infty; \infty)`
-  :math:`|b| < a` if :math:`p \ge 0`
-  :math:`|b| \le a` if :math:`p < 0`
-  :math:`K_{p}(.)` denotes the modified Bessel function of the second kind and order :math:`p` (`scipy.special.kn`)

The probability density above is defined in the "standardized" form. To shift and/or scale the distribution use the :math:`\text{loc}` and :math:`\text{scale}` parameters. Specifically, :math:`f(x, p, a, b, \text{loc}, \text{scale})` is identically equivalent to :math:`\frac{1}{\text{scale}}f(y, p, a, b)` with :math:`y = \frac{1}{\text{scale}}(x - \text{loc})`.

This parameterization derives from the original :math:`(\lambda, \alpha, \beta, \delta, \mu)` parameterization in  Barndorff (1978) by setting:

-  :math:`\lambda = p`
-  :math:`\alpha = \frac{a}{\delta} = \frac{\hat{\alpha}}{\delta}`
-  :math:`\beta = \frac{b}{\delta} = \frac{\hat{\beta}}{\delta}`
-  :math:`\delta = \text{scale}`
-  :math:`\mu = \text{location}`


Random variates for the `scipy.stats.genhyperbolic` can be efficiently sampled from the above-mentioned normal variance-mean mixture where `scipy.stats.geninvgauss` is parametrized as :math:`GIG\Big(p = p, b = \sqrt{\hat{\alpha}^2 - \hat{\beta}^2}, \text{loc} = \text{location}, \text{scale} = \frac{1}{\sqrt{\hat{\alpha}^2 - \hat{\beta}^2}}\Big)` so that: :math:`GH(p, \hat{\alpha}, \hat{\beta}) = \hat{\beta} \cdot GIG + \sqrt{GIG} \cdot N(0,1)`


The "generalized" characterization suggests the fact that this distribution is a superclass of several other probability distribution, for instance:

-  :math:`f(p = -\nu/2,  a = 0, b = 0, \text{loc} = 0, \text{scale} = \sqrt{\nu})` has a Student's t-distribution (`scipy.stats.t`) with :math:`\nu` degrees of freedom.
-  :math:`f(p = 1, a = \hat{\alpha}, b = \hat{\beta}, \text{loc} = \mu, \text{scale} = \delta)` has a Hyperbolic Distribution.
-  :math:`f(p = - 1/2, a = \hat{\alpha}, b = \hat{\beta}, \text{loc} = \mu, \text{scale} = \delta)` has a Normal Inverse Gaussian Distribution (`scipy.stats.norminvgauss`).
-  :math:`f(p = 1, a = \delta, b = 0, loc = \mu, \text{scale} = \delta)` has a Laplace Distribution (`scipy.stats.laplace`) for :math:`\delta \rightarrow 0`


Examples
--------

It is useful to understand how the parameters affect the shape of the distribution. While it is fairly straightforward to interpret the meaning of :math:`b` as the skewness, understanding the difference between :math:`a` and :math:`p` is not as obvious, as both affect the kurtosis of the distribution. :math:`a` can be interpreted as the speed of the decay of the probability density (where :math:`a > 1` the asymptotic decay is faster than :math:`log_e` and vice versa) or - equivalently - as the slope of the log-probability hyperbola asymptote (where :math:`a > 1` decay is faster than :math:`|1|` and vice versa). :math:`p` can be seen as the width of the shoulders of the probability density distribution (where :math:`p < 1` results in narrow shoulders and vice versa) or - equivalently - as the shape of the log-probability hyperbola, which is convex for :math:`p < 1` and concave otherwise.

.. code-block:: python

    import numpy as np
    from matplotlib import pyplot as plt
    from scipy import stats
    
    p, a, b, loc, scale = 1, 1, 0, 0, 1
    x = np.linspace(-10, 10, 100)
    
    # plot GH for different values of p
    plt.figure(0)
    plt.title("Generalized Hyperbolic | -10 < p < 10")
    plt.plot(x, stats.genhyperbolic.pdf(x, p, a, b, loc, scale),
            label = 'GH(p=1, a=1, b=0, loc=0, scale=1)')
    plt.plot(x, stats.genhyperbolic.pdf(x, p, a, b, loc, scale),
            color = 'red', alpha = 0.5, label='GH(p>1, a=1, b=0, loc=0, scale=1)')
    [plt.plot(x, stats.genhyperbolic.pdf(x, p, a, b, loc, scale),
            color = 'red', alpha = 0.1) for p in np.linspace(1, 10, 10)]
    plt.plot(x, stats.genhyperbolic.pdf(x, p, a, b, loc, scale),
            color = 'blue', alpha = 0.5, label='GH(p<1, a=1, b=0, loc=0, scale=1)')
    [plt.plot(x, stats.genhyperbolic.pdf(x, p, a, b, loc, scale),
            color = 'blue', alpha = 0.1) for p in np.linspace(-10, 1, 10)]
    plt.plot(x, stats.norm.pdf(x, loc, scale), label = 'N(loc=0, scale=1)')
    plt.plot(x, stats.laplace.pdf(x, loc, scale), label = 'Laplace(loc=0, scale=1)')
    plt.plot(x, stats.pareto.pdf(x+1, 1, loc, scale), label = 'Pareto(a=1, loc=0, scale=1)')
    plt.ylim(1e-15, 1e2)
    plt.yscale('log')
    plt.legend(bbox_to_anchor=(1.1, 1))
    plt.subplots_adjust(right=0.5)
    
    # plot GH for different values of a
    plt.figure(1)
    plt.title("Generalized Hyperbolic | 0 < a < 10")
    plt.plot(x, stats.genhyperbolic.pdf(x, p, a, b, loc, scale),
            label = 'GH(p=1, a=1, b=0, loc=0, scale=1)')
    plt.plot(x, stats.genhyperbolic.pdf(x, p, a, b, loc, scale),
            color = 'blue', alpha = 0.5, label='GH(p=1, a>1, b=0, loc=0, scale=1)')
    [plt.plot(x, stats.genhyperbolic.pdf(x, p, a, b, loc, scale),
            color = 'blue', alpha = 0.1) for a in np.linspace(1, 10, 10)]
    plt.plot(x, stats.genhyperbolic.pdf(x, p, a, b, loc, scale),
            color = 'red', alpha = 0.5, label='GH(p=1, 0<a<1, b=0, loc=0, scale=1)')
    [plt.plot(x, stats.genhyperbolic.pdf(x, p, a, b, loc, scale),
            color = 'red', alpha = 0.1) for a in np.linspace(0, 1, 10)]
    plt.plot(x, stats.norm.pdf(x, loc, scale),  label = 'N(loc=0, scale=1)')
    plt.plot(x, stats.laplace.pdf(x, loc, scale), label = 'Laplace(loc=0, scale=1)')
    plt.plot(x, stats.pareto.pdf(x+1, 1, loc, scale), label = 'Pareto(a=1, loc=0, scale=1)')
    plt.ylim(1e-15, 1e2)
    plt.yscale('log')
    plt.legend(bbox_to_anchor=(1.1, 1))
    plt.subplots_adjust(right=0.5)
    
    plt.show()

References
----------

-  Normal Variance-Mean Mixture
   https://en.wikipedia.org/wiki/Normal_variance-mean_mixture

-  Generalized Hyperbolic Distribution
   https://en.wikipedia.org/wiki/Generalised_hyperbolic_distribution

-  O. Barndorff-Nielsen, "Hyperbolic Distributions and Distributions
   on Hyperbolae", Scandinavian Journal of Statistics, Vol. 5(3),
   pp. 151-157, 1978. https://www.jstor.org/stable/4615705

-  Eberlein E., Prause K. (2002) The Generalized Hyperbolic Model:
   Financial Derivatives and Risk Measures. In: Geman H., Madan D.,
   Pliska S.R., Vorst T. (eds) Mathematical Finance - Bachelier
   Congress 2000. Springer Finance. Springer, Berlin, Heidelberg.
   https://doi.org/10.1007/978-3-662-12429-1_12

-  Scott, David J, Würtz, Diethelm, Dong, Christine and Tran,
   Thanh Tam, (2009), Moments of the generalized hyperbolic
   distribution, MPRA Paper, University Library of Munich, Germany,
   https://EconPapers.repec.org/RePEc:pra:mprapa:19081.

Implementation: `scipy.stats.genhyperbolic`
