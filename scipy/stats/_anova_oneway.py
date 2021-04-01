
from dataclasses import dataclass
from collections import namedtuple
import numpy as np
from scipy import special
from scipy.special import stdtrit
from ._anova_util import (_make_table, _encode, _proj)


__all__ = [
    'AnovaOnewayResult',
    'oneway', 'oneway_from_labels_values',
]


ConfidenceInterval = namedtuple('ConfidenceInterval', ['low', 'high'])
Mean = namedtuple('Mean', ['mean', 'se', 'ci'])


def _fmt(value):
    return np.array2string(np.array(value))


def _afmt(x):
    # This code assumes x is a 1d sequence of floats or ints
    s = [_fmt(t) for t in x]
    s = [t + ('.' if '.' not in t else '') for t in s]
    intpart, fracpart = zip(*[t.split('.') for t in s])
    intwidth = max(len(t) for t in intpart)
    fracwidth = max(len(t) for t in fracpart)
    jintpart = [t.rjust(intwidth) for t in intpart]
    jfracpart = [t.ljust(fracwidth) for t in fracpart]
    js = ['.'.join(parts) for parts in zip(jintpart, jfracpart)]
    return js


@dataclass
class AnovaOnewayResult:

    mean: float
    group_means: np.ndarray
    group_sizes: np.ndarray
    SSb: float
    SSw: float
    DFb: int
    DFw: int
    MSb: float
    MSw: float
    F: float
    pvalue: float
    levels: np.ndarray

    def __str__(self):
        vartot = self.SSb + self.SSw
        dftot = self.DFb + self.DFw

        header = "ANOVA one-way"
        rowlabels = ['Between groups', 'Within groups', 'Total']
        return _make_table(header, rowlabels,
                           [self.SSb, self.SSw, vartot],
                           [self.DFb, self.DFw, dftot],
                           [self.MSb, self.MSw],
                           [self.F],
                           [self.pvalue])

    def means_ci(self, confidence_level=0.95):
        """
        Confidence intervals for the means of the groups.
        """
        means = self.group_means
        ngroups = len(means)
        t = stdtrit(self.DFw, 0.5 + 0.5*confidence_level)
        cis = []
        for i in range(ngroups):
            se = np.sqrt(self.MSw * (1/self.group_sizes[i]))
            delta = t*se
            cis.append(Mean(mean=means[i],
                            se=se,
                            ci=ConfidenceInterval(low=means[i] - delta,
                                                  high=means[i] + delta)))
        return cis

    def deltas(self, confidence_level=0.95):
        """
        Confidence intervals for the differences of the means.
        """
        means = self.group_means
        ngroups = len(means)
        result = {}
        t = stdtrit(self.DFw, 0.5 + 0.5*confidence_level)
        for i in range(ngroups - 1):
            for j in range(i + 1, ngroups):
                se = np.sqrt(self.MSw * (1/self.group_sizes[i]
                                         + 1/self.group_sizes[j]))
                diff = means[i] - means[j]
                result[(i, j)] = Mean(mean=diff, se=se,
                                      ci=ConfidenceInterval(low=diff - t*se,
                                                            high=diff + t*se))
        return result


def oneway(*args, **kwds):
    """
    One-way analyis of variance.

    This is also know as single factor analysis of variance.

    Parameters
    ----------
    group1, group2, ... : 1-d sequences
        The groups of values to be tested.  Each group must
        be a one-dimensional sequence of values.

    Returns
    -------
    aov : AnovaOnewayResult
        An object whose attributes contain the information normally
        presented in an ANOVA table.  The object also has the methods
        `means_ci()` that computes the confidence intervals for the
        means of each group, and `deltas()`, the computes the confidence
        intervals for the pairwise differences in the means of the groups.

    See Also
    --------
    scipy.stats.anova.oneway_from_labels_values

    Notes
    -----
    .. versionadded:: 1.7.0

    References
    ----------
    .. [1] Jerrold H. Zar, Biostatistical Analysis (fifth edition),
           Prentice Hall, Upper Saddle River, New Jersey USA (2010)

    Examples
    --------
    >>> from scip.stats import anova

    This is Example 10.1 from Zar [1]_.

    >>> feed1 = [60.8, 67.0, 65.0, 68.6, 61.7]
    >>> feed2 = [68.7, 67.7, 75.0, 73.3, 71.8]
    >>> feed3 = [69.6, 77.1, 75.2, 71.5]
    >>> feed4 = [61.9, 64.2, 63.1, 66.7, 60.3]
    >>> aov = anova.oneway(feed1, feed2, feed3, feed4)

    When the `aov` object is printed, it is formatted as an
    ANOVA table.

    >>> print(aov)
    ANOVA one-way
    Source                   SS  DF            MS            F           p
    Between groups 338.93736842   3  112.97912281  12.04040385  0.00028301
    Within groups  140.75        15  9.38333333
    Total          479.68736842  18  26.64929825

    The p-value is 0.000283, which suggests that the means of the four
    populations from which the samples were drawn are not all the same.
    """
    num_groups = len(args)
    groups = [np.asarray(arg, dtype=np.float64) for arg in args]
    means = [group.mean() for group in groups]
    n = 0
    grand_total = 0
    for group in groups:
        n += len(group)
        grand_total += group.sum()
    grand_mean = grand_total / n

    v = sum(((group - grand_mean)**2).sum() for group in groups)
    vb = sum(len(group)*(group.mean() - grand_mean)**2 for group in groups)

    # Check for the edge case where the values in each group are constant.
    # When this happens, vw is 0.  If we don't handle this explicitly, vw
    # might contain numerical noise, and then F will be nonsense.
    if all([np.all(group[0] == group) for group in groups]):
        vw = 0.0
    else:
        vw = v - vb

    dfb = num_groups - 1
    dfw = n - num_groups
    msb = vb / dfb
    msw = vw / dfw
    if msw > 0:
        F = msb / msw
    else:
        F = np.inf
    p = special.fdtrc(dfb, dfw, F)
    result = AnovaOnewayResult(
                mean=grand_mean,
                group_means=means,
                group_sizes=[len(g) for g in groups],
                SSb=vb, SSw=vw,
                DFb=dfb, DFw=dfw,
                MSb=msb, MSw=msw,
                F=F, pvalue=p,
                levels=np.arange(num_groups))
    return result


# Old version that works by reorganizing the data and passing it to oneway().
def _oneway_from_labels_values(labels, values):
    """
    One-way analysis of variance.

    This is also known as single factor analysis of variance.

    This does the same calculation as `anova.oneway`.  The difference
    is in how the data is given to the function.

    Parameters
    ----------
    labels : array_like, 1-d
        Group labels of the data in `values`.
    values : array_like, 1-d
        Values associated with the labels in `labels`.

    Returns
    -------
    aov : AnovaOnewayResult
        An object whose attributes contain the information normally
        presented in an ANOVA table.  The object also has the methods
        `means_ci()` that computes the confidence intervals for the
        means of each group, and `deltas()`, the computes the confidence
        intervals for the pairwise differences in the means of the groups.
        This function addes the attribute `levels` to `aov`, containing
        the unique values in `labels` the define the groups.

    See Also
    --------
    scipy.stats.anova.oneway

    Notes
    -----
    .. versionadded:: 1.7.0

    Examples
    --------
    >>> from scipy.stats import anova
    >>> labels = ['A', 'B', 'B', 'C', 'B', 'A', 'A', 'C', 'A',
    ...           'C', 'C', 'A', 'B', 'B', 'C', 'A', 'B', 'B']
    >>> values = [4.5, 3.0, 3.9, 4.8, 4.1, 4.2, 4.0, 4.7, 3.9,
    ...           4.1, 4.7, 4.2, 4.1, 3.8, 4.5, 4.1, 4.0, 4.1]
    >>> aov = anova.oneway_from_labels_values(labels, values)
    >>> print(aov)
    ANOVA one-way
    Source                 SS  DF          MS           F          p
    Between groups 1.44085714   2  0.72042857  7.38072007  0.0058652
    Within groups  1.46414286  15  0.09760952
    Total          2.905       17  0.17088235

    Take a look at the means of the groups and their confidence intervals.

    >>> for lab, m in zip(aov.levels, aov.means_ci()):
    ...     print(f'{lab}: mean = {m.mean:5.3f}', end='')
    ...     print(f'  ci = ({m.ci.low:5.3f}, {m.ci.high:5.3f})')
    ...
    A: mean = 4.150  ci = (3.878, 4.422)
    B: mean = 3.857  ci = (3.605, 4.109)
    C: mean = 4.560  ci = (4.262, 4.858)

    Also inspect the pairwise differences of the means and their
    confidence intervals.

    >>> labs = aov.levels
    >>> for (i, j), m in aov.deltas().items():
    ...     print(f'{labs[i]} - {labs[j]}: mean diff = {m.mean:7.4f}', end='')
    ...     print(f'  ci = ({m.ci.low:7.4f}, {m.ci.high:7.4f})')
    ...
    A - B: mean diff =  0.2929  ci = (-0.0776,  0.6633)
    A - C: mean diff = -0.4100  ci = (-0.8132, -0.0068)
    B - C: mean diff = -0.7029  ci = (-1.0928, -0.3129)

    """
    values = np.asarray(values)
    levels, idx = np.unique(labels, return_inverse=True)
    groups = [values[idx == i] for i in range(len(levels))]
    result = oneway(*groups)
    result.levels = levels
    return result


def oneway_from_labels_values(labels, values):
    """
    One-way analysis of variance.

    This is also known as single factor analysis of variance.

    This does the same calculation as `anova.oneway`.  The difference
    is in how the data is given to the function.

    Parameters
    ----------
    labels : array_like, 1-d
        Group labels of the data in `values`.
    values : array_like, 1-d
        Values associated with the labels in `labels`.

    Returns
    -------
    aov : AnovaOnewayResult
        An object whose attributes contain the information normally
        presented in an ANOVA table.  The object also has the methods
        `means_ci()` that computes the confidence intervals for the
        means of each group, and `deltas()`, the computes the confidence
        intervals for the pairwise differences in the means of the groups.
        This function addes the attribute `levels` to `aov`, containing
        the unique values in `labels` the define the groups.

    See Also
    --------
    scipy.stats.anova.oneway

    Notes
    -----
    .. versionadded:: 1.7.0

    Examples
    --------
    >>> from scipy.stats import anova
    >>> labels = ['A', 'B', 'B', 'C', 'B', 'A', 'A', 'C', 'A',
    ...           'C', 'C', 'A', 'B', 'B', 'C', 'A', 'B', 'B']
    >>> values = [4.5, 3.0, 3.9, 4.8, 4.1, 4.2, 4.0, 4.7, 3.9,
    ...           4.1, 4.7, 4.2, 4.1, 3.8, 4.5, 4.1, 4.0, 4.1]
    >>> aov = anova.oneway_from_labels_values(labels, values)
    >>> print(aov)
    ANOVA one-way
    Source                 SS  DF          MS           F          p
    Between groups 1.44085714   2  0.72042857  7.38072007  0.0058652
    Within groups  1.46414286  15  0.09760952
    Total          2.905       17  0.17088235

    Take a look at the means of the groups and their confidence intervals.

    >>> for lab, m in zip(aov.levels, aov.means_ci()):
    ...     print(f'{lab}: mean = {m.mean:5.3f}', end='')
    ...     print(f'  ci = ({m.ci.low:5.3f}, {m.ci.high:5.3f})')
    ...
    A: mean = 4.150  ci = (3.878, 4.422)
    B: mean = 3.857  ci = (3.605, 4.109)
    C: mean = 4.560  ci = (4.262, 4.858)

    Also inspect the pairwise differences of the means and their
    confidence intervals.

    >>> labs = aov.levels
    >>> for (i, j), m in aov.deltas().items():
    ...     print(f'{labs[i]} - {labs[j]}: mean diff = {m.mean:7.4f}', end='')
    ...     print(f'  ci = ({m.ci.low:7.4f}, {m.ci.high:7.4f})')
    ...
    A - B: mean diff =  0.2929  ci = (-0.0776,  0.6633)
    A - C: mean diff = -0.4100  ci = (-0.8132, -0.0068)
    B - C: mean diff = -0.7029  ci = (-1.0928, -0.3129)

    """
    values = np.asarray(values)
    A, levels, inv = _encode(labels)

    grand_mean = np.mean(values)
    num_groups = len(levels)
    groups = [values[inv == k] for k in range(num_groups)]
    group_means = [np.mean(group) for group in groups]
    group_sizes = [len(group) for group in groups]

    n = len(values)
    I = np.eye(n)
    c = np.ones(n)

    # Full model
    X1 = np.column_stack((c, A))
    P1 = _proj(X1)

    SSbetween = ((P1 @ values - np.mean(values))**2).sum()

    # Compute SSwithin.  Note: SSwithin is also known as SSerror.
    # Check for the edge case where the values in each group are constant.
    # When this happens, SSwithin is 0.  If we don't handle this explicitly,
    # SSwithin might contain numerical noise, and then F will be nonsense.
    if all([np.all(group[0] == group) for group in groups]):
        SSwithin = 0.0
    else:
        SSwithin = values @ (I - P1) @ values

    DFbetween = num_groups - 1
    DFwithin = n - num_groups
    MSbetween = SSbetween / DFbetween
    MSwithin = SSwithin / DFwithin
    if MSwithin > 0:
        F = MSbetween / MSwithin
    else:
        F = np.inf
    p = special.fdtrc(DFbetween, DFwithin, F)

    result = AnovaOnewayResult(
                mean=grand_mean,
                group_means=group_means,
                group_sizes=group_sizes,
                SSb=SSbetween, SSw=SSwithin,
                DFb=DFbetween, DFw=DFwithin,
                MSb=MSbetween, MSw=MSwithin,
                F=F, pvalue=p,
                levels=levels)
    return result
