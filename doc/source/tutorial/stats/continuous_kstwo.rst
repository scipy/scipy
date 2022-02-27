
.. _continuous-kstwo:

KStwo Distribution
==================

This is the distribution of the maximum absolute differences between an
empirical distribution function, computed from :math:`n` samples or observations,
and a comparison (or target) cumulative distribution function, which is
assumed to be continuous.
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
--------

>>> from scipy.stats import kstwo

Show the probability of a gap at least as big as 0, 0.5 and 1.0 for a sample of size 5

>>> kstwo.sf([0, 0.5, 1.0], 5)
array([1.   , 0.112, 0.   ])

Compare a sample of size 5 drawn from a source N(0.5, 1) distribution against
a target N(0, 1) CDF.

>>> from scipy.stats import norm
>>> n = 5
>>> gendist = norm(0.5, 1)       # Normal distribution, mean 0.5, stddev 1
>>> x = np.sort(gendist.rvs(size=n, random_state=np.random.default_rng()))
>>> x
array([-1.59113056, -0.66335147,  0.54791569,  0.78009321,  1.27641365])
>>> target = norm(0, 1)
>>> cdfs = target.cdf(x)
>>> cdfs
array([0.0557901 , 0.25355274, 0.7081251 , 0.78233199, 0.89909533])
# Construct the Empirical CDF and the K-S statistics (Dn+, Dn-, Dn)
>>> ecdfs = np.arange(n+1, dtype=float)/n
>>> cols = np.column_stack([x, ecdfs[1:], cdfs, cdfs - ecdfs[:n], ecdfs[1:] - cdfs])
>>> np.set_printoptions(precision=3)
>>> cols
array([[-1.591,  0.2  ,  0.056,  0.056,  0.144],
       [-0.663,  0.4  ,  0.254,  0.054,  0.146],
       [ 0.548,  0.6  ,  0.708,  0.308, -0.108],
       [ 0.78 ,  0.8  ,  0.782,  0.182,  0.018],
       [ 1.276,  1.   ,  0.899,  0.099,  0.101]])
>>> gaps = cols[:, -2:]
>>> Dnpm = np.max(gaps, axis=0)
>>> Dn = np.max(Dnpm)
>>> iminus, iplus = np.argmax(gaps, axis=0)
>>> print('Dn- = %f (at x=%.2f)' % (Dnpm[0], x[iminus]))
Dn- = 0.308125 (at x=0.55)
>>> print('Dn+ = %f (at x=%.2f)' % (Dnpm[1], x[iplus]))
Dn+ = 0.146447 (at x=-0.66)
>>> print('Dn  = %f' % (Dn))
Dn  = 0.308125

>>> probs = kstwo.sf(Dn, n)
>>> print(chr(10).join(['For a sample of size %d drawn from a N(0, 1) distribution:' % n,
...      ' Kolmogorov-Smirnov 2-sided n=%d: Prob(Dn >= %f) = %.4f' % (n, Dn, probs)]))
For a sample of size 5 drawn from a N(0, 1) distribution:
 Kolmogorov-Smirnov 2-sided n=5: Prob(Dn >= 0.308125) = 0.6319

Plot the Empirical CDF against the target N(0, 1) CDF

>>> import matplotlib.pyplot as plt
>>> plt.step(np.concatenate([[-3], x]), ecdfs, where='post', label='Empirical CDF')
>>> x3 = np.linspace(-3, 3, 100)
>>> plt.plot(x3, target.cdf(x3), label='CDF for N(0, 1)')
>>> plt.ylim([0, 1]); plt.grid(True); plt.legend();
>>> plt.vlines([x[iminus]], ecdfs[iminus], cdfs[iminus], color='r', linestyle='solid', lw=4)
>>> plt.vlines([x[iplus]], cdfs[iplus], ecdfs[iplus+1], color='m', linestyle='solid', lw=4)
>>> plt.annotate('Dn-', xy=(x[iminus], (ecdfs[iminus]+ cdfs[iminus])/2),
...              xytext=(x[iminus]+1, (ecdfs[iminus]+ cdfs[iminus])/2 - 0.02),
...              arrowprops=dict(facecolor='white', edgecolor='r', shrink=0.05), size=15, color='r');
>>> plt.annotate('Dn+', xy=(x[iplus], (ecdfs[iplus+1]+ cdfs[iplus])/2),
...             xytext=(x[iplus]-2, (ecdfs[iplus+1]+ cdfs[iplus])/2 - 0.02),
...             arrowprops=dict(facecolor='white', edgecolor='m', shrink=0.05), size=15, color='m');
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
