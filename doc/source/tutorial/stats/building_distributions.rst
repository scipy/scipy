Building specific distributions
-------------------------------

The next examples shows how to build your own distributions. Further
examples show the usage of the distributions and some statistical
tests.


Making a continuous distribution, i.e., subclassing ``rv_continuous``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Making continuous distributions is fairly simple.

    >>> from scipy import stats
    >>> class deterministic_gen(stats.rv_continuous):
    ...     def _cdf(self, x):
    ...         return np.where(x < 0, 0., 1.)
    ...     def _stats(self):
    ...         return 0., 0., 0., 0.

    >>> deterministic = deterministic_gen(name="deterministic")
    >>> deterministic.cdf(np.arange(-3, 3, 0.5))
    array([ 0.,  0.,  0.,  0.,  0.,  0.,  1.,  1.,  1.,  1.,  1.,  1.])

Interestingly,  the ``pdf`` is now computed automatically:

    >>> deterministic.pdf(np.arange(-3, 3, 0.5))
    array([  0.00000000e+00,   0.00000000e+00,   0.00000000e+00,
             0.00000000e+00,   0.00000000e+00,   0.00000000e+00,
             5.83333333e+04,   4.16333634e-12,   4.16333634e-12,
             4.16333634e-12,   4.16333634e-12,   4.16333634e-12])


Be aware of the performance issues mentioned in
:ref:`performance_issues_label`. The computation of unspecified
common methods can become very slow, since only general methods are
called, which, by their very nature, cannot use any specific
information about the distribution. Thus, as a cautionary example:

    >>> from scipy.integrate import quad
    >>> quad(deterministic.pdf, -1e-1, 1e-1)
    (4.163336342344337e-13, 0.0)

But this is not correct: the integral over this pdf should be 1. Let's make the
integration interval smaller:

    >>> quad(deterministic.pdf, -1e-3, 1e-3)  # warning removed
    (1.000076872229173, 0.0010625571718182458)

This looks better. However, the problem originated from the fact that
the pdf is not specified in the class definition of the deterministic
distribution.


Subclassing ``rv_discrete``
^^^^^^^^^^^^^^^^^^^^^^^^^^^

In the following, we use `stats.rv_discrete` to generate a discrete
distribution that has the probabilities of the truncated normal for the
intervals centered around the integers.

**General info**

From the docstring of rv_discrete, ``help(stats.rv_discrete)``,

  "You can construct an arbitrary discrete rv where P{X=xk} = pk by
  passing to the rv_discrete initialization method (through the values=
  keyword) a tuple of sequences (xk, pk) which describes only those
  values of X (xk) that occur with nonzero probability (pk)."

Next to this, there are some further requirements for this approach to
work:

* The keyword `name` is required.
* The support points of the distribution xk have to be integers.
* The number of significant digits (decimals) needs to be specified.

In fact, if the last two requirements are not satisfied, an exception
may be raised or the resulting numbers may be incorrect.

**An example**

Let's do the work. First:

    >>> npoints = 20   # number of integer support points of the distribution minus 1
    >>> npointsh = npoints // 2
    >>> npointsf = float(npoints)
    >>> nbound = 4   # bounds for the truncated normal
    >>> normbound = (1+1/npointsf) * nbound   # actual bounds of truncated normal
    >>> grid = np.arange(-npointsh, npointsh+2, 1)   # integer grid
    >>> gridlimitsnorm = (grid-0.5) / npointsh * nbound   # bin limits for the truncnorm
    >>> gridlimits = grid - 0.5   # used later in the analysis
    >>> grid = grid[:-1]
    >>> probs = np.diff(stats.truncnorm.cdf(gridlimitsnorm, -normbound, normbound))
    >>> gridint = grid

And, finally, we can subclass ``rv_discrete``:

    >>> normdiscrete = stats.rv_discrete(values=(gridint,
    ...              np.round(probs, decimals=7)), name='normdiscrete')

Now that we have defined the distribution, we have access to all
common methods of discrete distributions.

    >>> print('mean = %6.4f, variance = %6.4f, skew = %6.4f, kurtosis = %6.4f' %
    ...       normdiscrete.stats(moments='mvsk'))
    mean = -0.0000, variance = 6.3302, skew = 0.0000, kurtosis = -0.0076

    >>> nd_std = np.sqrt(normdiscrete.stats(moments='v'))

**Testing the implementation**

Let's generate a random sample and compare observed frequencies with
the probabilities.

    >>> n_sample = 500
    >>> rvs = normdiscrete.rvs(size=n_sample)
    >>> f, l = np.histogram(rvs, bins=gridlimits)
    >>> sfreq = np.vstack([gridint, f, probs*n_sample]).T
    >>> print(sfreq)
    [[-1.00000000e+01  0.00000000e+00  2.95019349e-02]  # random
     [-9.00000000e+00  0.00000000e+00  1.32294142e-01]
     [-8.00000000e+00  0.00000000e+00  5.06497902e-01]
     [-7.00000000e+00  2.00000000e+00  1.65568919e+00]
     [-6.00000000e+00  1.00000000e+00  4.62125309e+00]
     [-5.00000000e+00  9.00000000e+00  1.10137298e+01]
     [-4.00000000e+00  2.60000000e+01  2.24137683e+01]
     [-3.00000000e+00  3.70000000e+01  3.89503370e+01]
     [-2.00000000e+00  5.10000000e+01  5.78004747e+01]
     [-1.00000000e+00  7.10000000e+01  7.32455414e+01]
     [ 0.00000000e+00  7.40000000e+01  7.92618251e+01]
     [ 1.00000000e+00  8.90000000e+01  7.32455414e+01]
     [ 2.00000000e+00  5.50000000e+01  5.78004747e+01]
     [ 3.00000000e+00  5.00000000e+01  3.89503370e+01]
     [ 4.00000000e+00  1.70000000e+01  2.24137683e+01]
     [ 5.00000000e+00  1.10000000e+01  1.10137298e+01]
     [ 6.00000000e+00  4.00000000e+00  4.62125309e+00]
     [ 7.00000000e+00  3.00000000e+00  1.65568919e+00]
     [ 8.00000000e+00  0.00000000e+00  5.06497902e-01]
     [ 9.00000000e+00  0.00000000e+00  1.32294142e-01]
     [ 1.00000000e+01  0.00000000e+00  2.95019349e-02]]


.. plot:: tutorial/examples/normdiscr_plot1.py
   :align: center
   :alt: "An X-Y histogram plot showing the distribution of random variates. A blue trace shows a normal bell curve. A blue bar chart perfectly approximates the curve showing the true distribution. A red bar chart representing the sample is well described by the blue trace but not exact."
   :include-source: 0


.. plot:: tutorial/examples/normdiscr_plot2.py
   :align: center
   :alt: "An X-Y histogram plot showing the cumulative distribution of random variates. A blue trace shows a CDF for a typical normal distribution. A blue bar chart perfectly approximates the curve showing the true distribution. A red bar chart representing the sample is well described by the blue trace but not exact."
   :include-source: 0


Next, we can test whether our sample was generated by our norm-discrete
distribution. This also verifies whether the random numbers were generated
correctly.

The chisquare test requires that there are a minimum number of observations
in each bin. We combine the tail bins into larger bins so that they contain
enough observations.

    >>> f2 = np.hstack([f[:5].sum(), f[5:-5], f[-5:].sum()])
    >>> p2 = np.hstack([probs[:5].sum(), probs[5:-5], probs[-5:].sum()])
    >>> ch2, pval = stats.chisquare(f2, p2*n_sample)

    >>> print('chisquare for normdiscrete: chi2 = %6.3f pvalue = %6.4f' % (ch2, pval))
    chisquare for normdiscrete: chi2 = 12.466 pvalue = 0.4090  # random

The pvalue in this case is high, so we can be quite confident that
our random sample was actually generated by the distribution.
