from dataclasses import dataclass
import itertools
import numpy as np
from scipy import special
from ._anova_util import (_nway_groups, _make_table, _unpack_2d_data_grid,
                          _encode, _interaction, _proj)


__all__ = [
    'AnovaTwowayResult',
    'twoway_from_a_b_values',
    'twoway_from_data_grid',
]

_EMPTY_CELL_MSG = 'There is an input combination that has no data.'
_UNBALANCED_MSG = 'the data is unbalanced, so `sstype` must be given'


def _validate_sstype(sstype):
    if sstype not in [None, 1, 2, 3]:
        raise ValueError('if given, `sstype` must be 1, 2 or 3.')


@dataclass
class AnovaTwowayResult:

    SSA: float
    SSB: float
    SSAB: float
    SSerror: float

    DFA: int
    DFB: int
    DFAB: int
    DFerror: int

    MSA: float
    MSB: float
    MSAB: float
    MSerror: float

    FA: float
    FB: float
    FAB: float

    pA: float
    pB: float
    pAB: float

    sstype: int

    levels: list

    def __str__(self):
        header = f'ANOVA two-way (SS type {self.sstype})'
        rowlabels = ['Factor A', 'Factor B', 'Interaction', 'Error']
        return _make_table(header, rowlabels,
                           [self.SSA, self.SSB, self.SSAB, self.SSerror],
                           [self.DFA, self.DFB, self.DFAB, self.DFerror],
                           [self.MSA, self.MSB, self.MSAB, self.MSerror],
                           [self.FA, self.FB, self.FAB],
                           [self.pA, self.pB, self.pAB])


@dataclass
class AnovaTwoway1Result:

    SSA: float
    SSB: float
    SSerror: float

    DFA: int
    DFB: int
    DFerror: int

    MSA: float
    MSB: float
    MSerror: float

    FA: float
    FB: float

    pA: float
    pB: float

    levels: list

    def __str__(self):
        header = 'ANOVA two-way'
        rowlabels = ['Factor A', 'Factor B', 'Error']
        return _make_table(header, rowlabels,
                           [self.SSA, self.SSB, self.SSerror],
                           [self.DFA, self.DFB, self.DFerror],
                           [self.MSA, self.MSB, self.MSerror],
                           [self.FA, self.FB],
                           [self.pA, self.pB])


def _twoway_balanced(data):
    """
    Two-way ANOVA with balanced replication.

    The number of replications for each pair of factors must be the same.

    Parameters
    ----------
    data : array-like, either (m, n) or (m, n, r)
        `r` is the number of replications.  If `data` has shape (m, n),
        it implies `r` is 1.

    """
    # In the following, the two factors are "labeled" A and B.

    data = np.asarray(data)
    shp = data.shape
    if data.ndim == 2:
        raise ValueError("ndim = 2 not implemented yet, use anova.twoway1 "
                         "instead.")

    grand_mean = data.mean()

    mean2 = data.mean(axis=2, keepdims=True)

    meanB = data.mean(axis=(0, 2), keepdims=True)
    meanA = data.mean(axis=(1, 2), keepdims=True)

    ssB = shp[0]*shp[2]*((meanB - grand_mean)**2).sum()
    dofB = shp[1] - 1
    msB = ssB / dofB

    ssA = shp[1]*shp[2]*((meanA - grand_mean)**2).sum()
    dofA = shp[0] - 1
    msA = ssA / dofA

    ss_inter = shp[2]*((mean2 - meanA - meanB + grand_mean)**2).sum()
    dof_inter = (shp[0] - 1)*(shp[1] - 1)
    ms_inter = ss_inter / dof_inter

    # These are from R. Johnson "Miller & Freund's Prob. & Stats for Engineers"
    # ss_error  = ((data - mean2 - mean01 + grand_mean)**2).sum()
    # dof_error = (shp[0]*shp[1] - 1)*(shp[2] - 1)
    #
    # These are from Zar (fifth ed.)
    ss_error = ((data - mean2)**2).sum()
    dof_error = (shp[0]*shp[1])*(shp[2] - 1)
    ms_error = ss_error / dof_error

    FB = msB / ms_error
    FA = msA / ms_error
    F_inter = ms_inter / ms_error

    pA = special.fdtrc(dofA, dof_error, FA)
    pB = special.fdtrc(dofB, dof_error, FB)
    p_inter = special.fdtrc(dof_inter, dof_error, F_inter)

    result = AnovaTwowayResult(SSB=ssB, SSA=ssA,
                               SSAB=ss_inter, SSerror=ss_error,
                               DFB=dofB, DFA=dofA,
                               DFAB=dof_inter, DFerror=dof_error,
                               MSB=msB, MSA=msA,
                               MSAB=ms_inter, MSerror=ms_error,
                               FB=FB, FA=FA, FAB=F_inter,
                               pB=pB, pA=pA, pAB=p_inter,
                               sstype=1,
                               levels=[np.arange(m) for m in shp[:2]])
    return result


def _twoway1(data):
    """
    Two-way anova without replication.

    Parameters
    ----------
    data : array-like with shape (m, n)

    """
    data = np.asarray(data)
    shp = data.shape
    if data.ndim != 2:
        raise ValueError("This function is for two-way ANOVA with no "
                         "replication.")
    r, c = shp

    grand_mean = data.mean()
    mean0 = data.mean(axis=0, keepdims=True)
    mean1 = data.mean(axis=1, keepdims=True)

    ss_total = ((data - grand_mean)**2).sum()

    ss0 = r*((mean0 - grand_mean)**2).sum()
    ss1 = c*((mean1 - grand_mean)**2).sum()

    df0 = c - 1
    df1 = r - 1
    ms0 = ss0 / df0
    ms1 = ss1 / df1

    sse = ss_total - ss0 - ss1
    dfe = (c - 1)*(r - 1)
    mse = sse / dfe

    F1 = ms1 / mse
    F0 = ms0 / mse

    p1 = special.fdtrc(df1, dfe, F1)
    p0 = special.fdtrc(df0, dfe, F0)

    result = AnovaTwoway1Result(SSB=ss0, SSA=ss1, SSerror=sse,
                                DFB=df0, DFA=df1, DFerror=dfe,
                                MSB=ms0, MSA=ms1, MSerror=mse,
                                FB=F0, FA=F1,
                                pB=p0, pA=p1,
                                levels=[np.arange(m) for m in shp])
    return result


def twoway_from_a_b_values(a, b, values, sstype=None):
    """
    Two-way analyis of variance.

    Parameters
    ----------
    a, b : array_like, 1-d sequences
        These give the levels or labels of the A and B factors associated
        with each value in `values`.
    values : array_like, 1-d
        The array of observations or measurements.
    sstype : int
        Type of sum-of-squares.  Must be 1, 2, or 3.  If the data is
        unbalanced, this parameter is required.  If the data is balanced,
        the ANOVA result is the same for any of the three types, so the
        parameter is not required.

    Returns
    -------
    aov : `AnovaTwoWayResult` or `AnovaTwoway1Result`
        An object whose attributes contain the information normally
        presented in an ANOVA table.  When the object is printed, the
        output is formatted as an ANOVA table.

        If there is exactly one value per grid cell, an instance of
        `AnovaTwoway1Result` is returned.  This object does not include
        results for the interaction of the factors.

    Notes
    -----
    .. versionadded:: 1.7.0

    Examples
    --------
    The response times of three subjects are measured multiple times
    for three different tests, with the following results::

            Subject  Test    Time    Subject   Test   Time
               A       1      850       B        2     749
               A       1      912       B        3     801
               A       1      799       B        3     676
               A       2      651       B        3     703
               A       2      777       B        3     652
               A       3      545       C        1     788
               A       3      629       C        1     695
               A       3      633       C        1     513
               B       1      788       C        2     560
               B       1      834       C        2     629
               B       1      874       C        2     521
               B       2      749       C        3     483
               B       2      661       C        3     551
               B       2      589       C        3     613

    >>> from scipy.stats import anova

    We'll enter the data from each column as a list, and use
    `twoway_from_a_b_values` to perform the two-way analysis of
    variance.

    >>> subj = ['A', 'A', 'A', 'A', 'A', 'A', 'A', 'A',
    ...         'B', 'B', 'B', 'B', 'B', 'B', 'B', 'B', 'B', 'B', 'B',
    ...         'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C']
    >>> test = [1, 1, 1, 2, 2, 3, 3, 3,
    ...         1, 1, 1, 2, 2, 2, 2, 3, 3, 3, 3,
    ...         1, 1, 1, 2, 2, 2, 3, 3, 3]
    >>> time = [850, 912, 799, 651, 777, 545, 629, 633,
    ...         788, 834, 874, 749, 661, 589, 749, 801, 676, 703, 652,
    ...         788, 695, 513, 560, 629, 521, 483, 551, 613]

    Because the number of measurements per subject and test is not the
    same, we must specify the type of the two-way ANOVA to perform.
    Here we will perform a type 3 ANOVA.

    >>> aov = anova.twoway_from_a_b_values(subj, test, values=time, sstype=3)
    >>> print(aov)
    ANOVA two-way (SS type 3)
    Source                   SS DF             MS           F          p
    Factor A     119712.0186916  2  59856.0093458 10.59319452 0.00081172
    Factor B    135886.85228153  2 67943.42614076 12.02448906 0.00042219
    Interaction  23066.37699044  4  5766.59424761   1.0205601 0.42186674
    Error       107357.99999999 19  5650.42105263

    """
    _validate_sstype(sstype)
    if len(a) != len(b):
        raise ValueError('`a` and `b` must have the same length.')
    if len(values) != len(a):
        raise ValueError('`values` must have the same length as `a`.')
    levels, groups = _nway_groups(a, b, values=values)
    groups1d = groups.ravel()
    if any(g is None for g in groups1d):
        raise ValueError(_EMPTY_CELL_MSG)
    lengths = [len(t) for t in groups1d]
    if all(n == 1 for n in lengths):
        # One replicate; do ANOVA without interaction.
        g = np.array(groups.tolist())
        result = _twoway1(g[:, :, 0])
    else:
        len0 = lengths[0]
        if any(lenk != len0 for lenk in lengths[1:]) and sstype is None:
            raise ValueError(_UNBALANCED_MSG)
        result = _twoway_lm(a, b, values, sstype)
    result.levels = levels
    return result


def _twoway_lm(a, b, values, sstype):
    """
    `sstype` must be one of 1, 2 or 3.
    """
    if len(a) != len(b):
        raise ValueError("'a' and 'b' must have the same length")
    n = len(values)
    if len(a) != n:
        raise ValueError("'a' and 'values' must have the same length")

    A, levelsA, _ = _encode(a)
    B, levelsB, _ = _encode(b)
    AB = _interaction(A, B)

    I = np.eye(n)
    c = np.ones(n)

    # Full model
    X1 = np.column_stack((c, A, B, AB))
    P1 = _proj(X1)

    SSerror = values @ (I - P1) @ values

    # No interaction
    X2 = np.column_stack((c, A, B))
    P2 = _proj(X2)

    SSab = values @ (P1 - P2) @ values

    if sstype == 1:
        X5 = np.column_stack((c, A))
        P5 = _proj(X5)
        SSa = ((P5 @ values - np.mean(values))**2).sum()
        SSb = values @ (P2 - P5) @ values
    elif sstype == 2:
        X5 = np.column_stack((c, A))
        P5 = _proj(X5)
        X6 = np.column_stack((c, B))
        P6 = _proj(X6)
        SSa = values @ (P2 - P6) @ values
        SSb = values @ (P2 - P5) @ values
    else:
        # sstype == 3
        X4 = np.column_stack((c, B, AB))
        P4 = _proj(X4)
        SSa = values @ (P1 - P4) @ values
        X3 = np.column_stack((c, A, AB))
        P3 = _proj(X3)
        SSb = values @ (P1 - P3) @ values

    dfa = A.shape[1]
    dfb = B.shape[1]
    dfab = dfa*dfb
    dferror = len(values) - (dfa + 1)*(dfb + 1)

    MSa = SSa / dfa
    MSb = SSb / dfb
    MSab = SSab / dfab
    MSerror = SSerror / dferror

    Fa = MSa / MSerror
    Fb = MSb / MSerror
    Fab = MSab / MSerror

    pa = special.fdtrc(dfa, dferror, Fa)
    pb = special.fdtrc(dfb, dferror, Fb)
    pab = special.fdtrc(dfab, dferror, Fab)

    result = AnovaTwowayResult(SSA=SSa, SSB=SSb, SSAB=SSab, SSerror=SSerror,
                               DFA=dfa, DFB=dfb, DFAB=dfab, DFerror=dferror,
                               MSA=MSa, MSB=MSb, MSAB=MSab, MSerror=MSerror,
                               FA=Fa, FB=Fb, FAB=Fab,
                               pA=pa, pB=pb, pAB=pab,
                               sstype=sstype,
                               levels=[levelsA, levelsB])
    return result


def twoway_from_data_grid(data, sstype=None):
    """
    Two-way analyis of variance.

    Perform the two-way (or two factor) analyis of variance calculation.

    Parameters
    ----------
    data : array_like with shape (m, n) or (m, n, r), OR nested sequences
           that looks like an array with shape (m, n, r(m, n)), where
           r(m, n) represents the number of replicates in row m and column n.
    sstype : int
        Type of sum-of-squares.  Must be 1, 2, or 3.  If the data is
        unbalanced, this parameter is required.  If the data is balanced,
        the ANOVA result is the same for any of the three types, so the
        parameter is not required.

    Returns
    -------
    aov : `AnovaTwoWayResult` or `AnovaTwoway1Result`
        An object whose attributes contain the information normally
        presented in an ANOVA table.  When the object is printed, the
        output is formatted as an ANOVA table.

        If there is exactly one value per grid cell, an instance of
        `AnovaTwoway1Result` is returned.  This object does not include
        results for the interaction of the factors.

    Notes
    -----
    .. versionadded:: 1.7.0

    Examples
    --------
    The performance of three brands of oil is tested in four machines.
    Each machine is run with each oil, and after a fixed amount of time,
    the quality of the oil is tested.  For each oil and each machine, the
    test is done three times, so we have a total of 36 measurements::

               |                     Machine                    |
          -----+------------+-----------+-----------+-----------+
          Oil  |      1     |     2     |     3     |     4     |
          -----+------------+-----------+-----------+-----------+
          A    |  12 12 13  | 10 19 13  | 14 12 12  | 15 16 14  |
          B    |  14 16 14  | 18 15 14  | 13 15 12  | 16 16 17  |
          C    |  10 12 13  | 11 15 16  | 12 17 13  | 17 18 19  |
          -----+------------+-----------+-----------+-----------+

    We perform a two-way analysis of variance on the measurements.

    >>> from scipy.stats import anova

    We enter the data as a nested list; the innermost lists hold the
    measurements.  (We could make this a NumPy array, but that is not
    required.)

    >>> data = [[[12, 12, 13], [10, 19, 13], [14, 12, 12], [15, 16, 14]],
    ...         [[14, 16, 14], [18, 15, 14], [13, 15, 12], [16, 16, 17]],
    ...         [[10, 12, 13], [11, 15, 18], [12, 17, 13], [17, 18, 19]]]

    Compute and print the two-way analysis of variance.

    >>> aov = anova.twoway_from_data_grid(data)
    >>> print(aov)
    ANOVA two-way (SS type 1)
    Source                SS DF          MS          F          p
    Factor A     14.38888889  2  7.19444444 1.57926829 0.22680808
    Factor B     69.63888889  3 23.21296296 5.09552846 0.00718026
    Interaction  20.94444444  6  3.49074074 0.76626016 0.60363198
    Error       109.33333333 24  4.55555556

    """
    _validate_sstype(sstype)

    row_lengths = [len(row) for row in data]
    if len(row_lengths) < 2:
        raise ValueError('at least two rows required')
    if any(len(t) < 2 for t in data):
        raise ValueError('at least two columns required')
    if any(t != row_lengths[0] for t in row_lengths[1:]):
        raise ValueError('each row in `data` must have the same length')

    some_scalar = False
    all_scalar = True
    for row in data:
        for cell in row:
            if np.isscalar(cell):
                some_scalar = True
            else:
                all_scalar = False
    if some_scalar and not all_scalar:
        raise ValueError('inner elements of `data` must be all scalars or'
                         ' all one-dimensional sequences')

    if all_scalar:
        return _twoway1(data)

    # Check for empty cells, and determine if the data is balanced.
    len00 = len(data[0][0])
    unbalanced = False
    for i, j in itertools.product(range(len(data)), range(len(data[0]))):
        cell = data[i][j]
        if len(cell) == 0:
            # We require at least one replicate in each cell.
            raise ValueError(_EMPTY_CELL_MSG)
        if len(cell) != len00:
            unbalanced = True

    if unbalanced:
        if sstype is None:
            raise ValueError(_UNBALANCED_MSG)
        # Unpack `data` into 1d sequences `a`, `b` and `values`,
        # and pass to _twoway_lm
        a, b, values = _unpack_2d_data_grid(data)
        return _twoway_lm(a, b, values, sstype=sstype)

    # Balanced data, more than 1 replicate.
    return _twoway_balanced(data)
