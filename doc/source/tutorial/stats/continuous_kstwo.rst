
.. _continuous-kstwo:

KStwo Distribution
==================

This is the distribution of the maximum absolute differences between an
empirical distribution function, computed from :math:`n` samples or observations,
and a comparison (or target) cumulative distribution function.
(The "two" in the name is because this is the two-sided difference.
``ksone`` is the distribution of the positive differences, :math:`D_n^+`,
hence it concerns one-sided differences.
``kstwobign`` is the limiting
distribution of the *normalized* maximum absolute differences :math:`\sqrt{n} D_n`.)


Writing :math:`D_n = \sup_t \left|F_{empirical,n}(t)-F_{target}(t)\right|`,
``kstwo`` is the distribution of the :math:`D_n` values.


``kstwo`` can also be used with the differences between two empirical distribution functions,
for sets of observations with :math:`m` and :math:`n` samples respectively.
Writing :math:`D_{m,n} = \sup_t \left|F_{1,m}(t)-F_{2,n}(t)\right|`,  where
:math:`F_{1,m}` and :math:`F_{2,n}` are the two empirical distribution functions, then
:math:`Pr(D_{m,n} \le x) \approx Pr(D_N \le x)` under appropriate conditions,
where :math:`N = \sqrt{\left(\frac{mn}{m+n}\right)}`.


There is one shape parameter :math:`n`, a positive integer, and the support is :math:`x\in\left[0,1\right]`.

The implementation follows Simard & L'Ecuyer, which combines exact algorithms of Durbin and Pomeranz
with asymptotic estimates of Li-Chien, Pelz and Good to compute the CDF with 5-15 accurate digits.

Examples
========

>>> from scipy.stats import kstwo

Show the probability of a gap at least as big as 0, 0.5 and 1.0 for a sample of size 5

>>> kstwo.sf([0, 0.5, 1.0], 5)
array([1.   , 0.112, 0.   ])

Compare a sample of size 5 drawn from a source N(0.5, 1) distribution against
a target N(0, 1) CDF.

>>> from scipy.stats import norm
>>> n = 5
>>> gendist = norm(0.5, 1)       # Normal distribution, mean 0.5, stddev 1
>>> np.random.seed(seed=233423)  # Set the seed for reproducibility
>>> x = np.sort(gendist.rvs(size=n))
>>> x
array([-0.20946287,  0.71688765,  0.95164151,  1.44590852,  3.08880533])
>>> target = norm(0, 1)
>>> cdfs = target.cdf(x)
>>> cdfs
array([ 0.41704346,  0.76327829,  0.82936059,  0.92589857,  0.99899518])
# Construct the Empirical CDF and the K-S statistics (Dn+, Dn-, Dn)
>>> ecdfs = np.arange(n+1, dtype=float)/n
>>> cols = np.column_stack([x, ecdfs[1:], cdfs, cdfs - ecdfs[:n], ecdfs[1:] - cdfs])
>>> np.set_printoptions(precision=3)
>>> cols
array([[ -2.095e-01,   2.000e-01,   4.170e-01,   4.170e-01,  -2.170e-01],
       [  7.169e-01,   4.000e-01,   7.633e-01,   5.633e-01,  -3.633e-01],
       [  9.516e-01,   6.000e-01,   8.294e-01,   4.294e-01,  -2.294e-01],
       [  1.446e+00,   8.000e-01,   9.259e-01,   3.259e-01,  -1.259e-01],
       [  3.089e+00,   1.000e+00,   9.990e-01,   1.990e-01,   1.005e-03]])
>>> gaps = cols[:, -2:]
>>> Dnpm = np.max(gaps, axis=0)
>>> Dn = np.max(Dnpm)
>>> print('Dn-=%f, Dn+=%f, Dn=%f' % (Dnpm[0], Dnpm[1], Dn))
Dn-=0.563278, Dn+=0.001005, Dn=0.563278
>>> probs = kstwo.sf(Dn, n)
>>> print(chr(10).join(['For a sample of size %d drawn from a N(0, 1) distribution:' % n,
...      ' Kolmogorov-Smirnov 2-sided n=%d: Prob(Dn >= %f) = %.4f' % (n, Dn, probs)]))
For a sample of size 5 drawn from a N(0, 1) distribution:
 Kolmogorov-Smirnov 2-sided n=5: Prob(Dn >= 0.563278) = 0.0500

Plot the Empirical CDF against the target N(0, 1) CDF

>>> import matplotlib.pyplot as plt
>>> plt.step(np.concatenate([[-3], x]), ecdfs, where='post', label='Empirical CDF')
>>> x3 = np.linspace(-3, 3, 100)
>>> plt.plot(x3, target.cdf(x3), label='CDF for N(0, 1)')
>>> plt.ylim([0, 1]); plt.grid(True); plt.legend();
# Add vertical lines marking Dn+ and Dn-
>>> iminus, iplus = np.argmax(gaps, axis=0)
>>> plt.vlines([x[iminus]], ecdfs[iminus], cdfs[iminus], color='r', linestyle='dashed', lw=4)
>>> plt.vlines([x[iplus]], cdfs[iplus], ecdfs[iplus+1], color='m', linestyle='dashed', lw=4)
>>> plt.show()


References
----------

-  "Kolmogorov-Smirnov test", Wikipedia
   https://en.wikipedia.org/wiki/Kolmogorov-Smirnov_test

-  Durbin J. "The Probability that the Sample Distribution Function Lies Between Two
   Parallel Straight Lines." *Ann. Math. Statist*., 39 (1968) 39, 398-411.

-  Pomeranz J.  "Exact Cumulative Distribution of the Kolmogorov-Smirnov Statistic for
   Small Samples (Algorithm 487)."  *Communications of the ACM*, 17(12), (1974) 703-704.

-  Li-Chien, C.  "On the exact distribution of the statistics of A. N. Kolmogorov and
   their asymptotic expansion."  *Acta Matematica Sinica*, 6, (1956) 55-81.

-  Pelz W, Good IJ. "Approximating the Lower Tail-areas of the Kolmogorov-Smirnov One-sample
   Statistic." *Journal of the Royal Statistical Society*, Series B, (1976) 38(2), 152-156.

-  Simard, R., L'Ecuyer, P. "Computing the Two-Sided Kolmogorov-Smirnov Distribution",
   *Journal of Statistical Software*, Vol 39, (2011) 11.

Implementation: `scipy.stats.kstwo`
