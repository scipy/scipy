.. _resampling-and-monte-carlo-introduction:

Resampling and Monte Carlo Methods
==================================

.. currentmodule:: scipy.stats.resampling

Introduction
------------

Resampling and Monte Carlo methods are statistical techniques that
replace mathematical analysis with lots of computation.

For example, suppose you and your brother Kyle find yourselves
hitchhiking down a long and lonesome road. All of a sudden, there shines
a shiny demon... in the middle... of the road. And he says:

   If you flip a coin with probability of heads :math:`p=0.5` exactly
   :math:`n=100` times, what is the probability that the number of heads
   will be less than or equal to :math:`x=45`? Answer correctly, or I'll
   eat your souls.

    >>> import math
    >>> import numpy as np
    >>> p = 0.5  # probability of flipping heads each flip
    >>> n = 100  # number of coin flips per trial
    >>> x = 45  # we want to know the probability that the number of heads per trial will be less than or equal to this

Your brother Kyle is the analytical one. He answers:

   As the number of coin tosses increases, the distribution of the
   number of heads will tend towards normality with mean
   :math:`\mu = p n` and standard deviation
   :math:`\sigma = \sqrt{n p (1 - p)}`, where :math:`p = 0.5` is the
   probability of heads and :math:`n=100` is the number of flips. The
   probability of :math:`x=45` heads can be approximated as the
   cumulative density function :math:`F(x)` of this normal distribution.
   Specifically:

.. math::

    F(x; \mu, \sigma) = \frac{1}{2} \left[ 1 + \mbox{erf} \left( \frac{x-\mu}{\sigma \sqrt{2}} \right) \right]

.. code-block::

    >>> # Kyle's Analytical Approach
    >>> mean = p*n
    >>> std = math.sqrt(n*p*(1-p))
    >>> # CDF of the normal distribution. (Unfortunately, Kyle forgets a continuity correction that would produce a more accurate answer.)
    >>> prob = 0.5 * (1 + math.erf((x - mean) / (std * math.sqrt(2))))
    >>> print(f"The normal approximation estimates the probability as {prob}")
    The normal approximation estimates the probability as 0.15865525393145713

You are a little more practical, so you decide to take a computational
approach (or more precisely, a Monte Carlo approach): just simulate many
sequences of coin tosses, count the number of heads in each toss,
and estimate the probability as the fraction of sequences in which the
count does not exceed 45.

    >>> # Your Monte Carlo Approach
    >>> N = 100000  # We'll do 100000 trials, each with 100 flips
    >>> rng = np.random.default_rng()  # use the "new" Generator interface
    >>> simulation = rng.random(size=(n, N)) < p  # False for tails, True for heads
    >>> counts = np.sum(simulation, axis=0)  # count the number of heads each trial
    >>> prob = np.sum(counts <= x) / N  # estimate the probability as the observed proportion of cases in which the count did not exceed 45
    >>> print(f"The Monte Carlo approach estimates the probability as {prob}")
    The Monte Carlo approach estimates the probability as 0.18348

The demon replies:

   You are both incorrect. The probability is given by the binomial
   distribution. Specifically.

.. math::

    \sum_{i=0}^{x} {n \choose i} p^i (1-p)^{n-i}

.. code-block::

    >>> # The Demon's Exact Probability
    >>> from scipy.stats import binom
    >>> prob = binom.cdf(x, n, p)
    >>> print(f"The correct answer is approximately {prob}")
    The correct answer is approximately 0.18410080866334788

As your soul is being eaten, you take solace in the knowledge that your
simple Monte Carlo approach was more accurate than the normal
approximation. This is not uncommon: when an exact answer is unknown,
often a computational approximation is more accurate than an analytical
approximation. Also, it's easy for demons to invent questions for which
analytical approximations (let alone exact answers) are unavailable. In
such cases, a computational approach is the only way to go.

Resampling and Monte Carlo methods tutorials
--------------------------------------------

Although it's best to use an exact approach when it's available,
learning to use computational statistics techniques can improve the
accuracy of `scipy.stats` features that rely on analytical
approximations, dramatically extend your statistical analysis
capabilities, and even improve your understanding of statistics.
The following tutorials will help you get started with the resampling
and Monte Carlo methods in `scipy.stats`.

1. `Monte Carlo Hypothesis Tests <https://nbviewer.org/github/scipy/scipy-cookbook/blob/main/ipython/ResamplingAndMonteCarloMethods/resampling_tutorial_1.ipynb>`_
2. `Permutation Tests <https://nbviewer.org/github/scipy/scipy-cookbook/blob/main/ipython/ResamplingAndMonteCarloMethods/resampling_tutorial_2.ipynb>`_

   a. `Independent-Sample Permutation Tests <https://nbviewer.org/github/scipy/scipy-cookbook/blob/main/ipython/ResamplingAndMonteCarloMethods/resampling_tutorial_2a.ipynb>`_
   b. `Paired-Sample Permutation Tests <https://nbviewer.org/github/scipy/scipy-cookbook/blob/main/ipython/ResamplingAndMonteCarloMethods/resampling_tutorial_2b.ipynb>`_
   c. `Correlated-Sample Permutation Tests <https://nbviewer.org/github/scipy/scipy-cookbook/blob/main/ipython/ResamplingAndMonteCarloMethods/resampling_tutorial_2c.ipynb>`_

3. `The Bootstrap <https://nbviewer.org/github/scipy/scipy-cookbook/blob/main/ipython/ResamplingAndMonteCarloMethods/resampling_tutorial_3.ipynb>`_
