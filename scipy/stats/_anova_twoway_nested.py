from dataclasses import dataclass
import numpy as np
from scipy.special import fdtrc
from ._anova_util import _make_table, _encode, _encode_nested, _proj


__all__ = [
    'AnovaTwowayNestedResult',
    'twoway_nested_from_data_grid',
    'twoway_nested_from_a_b_values',
]


@dataclass
class AnovaTwowayNestedResult:

    SSA: float
    SSBA: float
    SSerror: float

    DFA: int
    DFBA: int
    DFerror: int

    MSA: float
    MSBA: float
    MSerror: float

    FA: float
    FBA: float

    pA: float
    pBA: float

    btype: str

    def __str__(self):
        vartot = self.SSA + self.SSBA + self.SSerror
        dftot = self.DFA + self.DFBA + self.DFerror

        header = ("ANOVA two-way, factor B nested in factor"
                  f" A, B is {self.btype}")
        rowlabels = ['Factor A', 'Factor B(A)', 'Error', 'Total']
        return _make_table(header, rowlabels,
                           [self.SSA, self.SSBA, self.SSerror, vartot],
                           [self.DFA, self.DFBA, self.DFerror, dftot],
                           [self.MSA, self.MSBA, self.MSerror],
                           [self.FA, self.FBA],
                           [self.pA, self.pBA])


def twoway_nested_from_data_grid(data, btype='random'):
    """
    Balanced, fully nested two-way ANOVA with factor B nested in factor A.

    The number of replications for each pair of factors must be the same.

    Parameters
    ----------
    data : array-like, either (m, n) or (m, n, r)
        `r` is the number of replications.  If `data` has shape (m, n),
        it implies `r` is 1.
    btype : str, optional
        Effect type of factor B. Must be one of ``'random'`` or ``'fixed'``.
        Default is ``'random'``.

    Returns
    -------
    aov : `AnovaTwoWayNestedResult`
        An object whose attributes contain the information normally
        presented in an ANOVA table.  When the object is printed, the
        output is formatted as an ANOVA table.

    Notes
    -----
    .. versionadded:: 1.7.0

    References
    ----------
    .. [1] NIST/SEMATECH e-Handbook of Statistical Methods,
           http://www.itl.nist.gov/div898/handbook/,
           https://doi.org/10.18434/M32189

    Examples
    --------
    This example is from section 3.2.3.3 of the NIST Engineering Statistics
    Handbook [1]_.

    There are five different machines, each with two operators, one for the
    day shift and one for the night shift.  Five samples are taken from each
    machine for each operator::

                                 Machine
                      1       2       3       4       5
                    0.125   0.118   0.123   0.126   0.118
          Day       0.127   0.122   0.125   0.128   0.129
          Operator  0.125   0.120   0.125   0.126   0.127
                    0.126   0.124   0.124   0.127   0.120
                    0.128   0.119   0.126   0.129   0.121

                    0.124   0.116   0.122   0.126   0.125
          Night     0.128   0.125   0.121   0.129   0.123
          Operator  0.127   0.119   0.124   0.125   0.114
                    0.126   0.125   0.126   0.130   0.124
                    0.129   0.120   0.125   0.124   0.117

    Because the operators are different for each machine, this is not
    a crossed data set.  The main factor is the machine, and the nested
    factor is the operator.  To arrange this data for the function
    `twoway_nested_from_data_grid`, we must put it in an array such that
    first dimension is the machine, the second is the operator, and the
    third is the set of measurements:

    >>> data = [[[0.125, 0.127, 0.125, 0.126, 0.128],
    ...          [0.124, 0.128, 0.127, 0.126, 0.129]],
    ...         [[0.118, 0.122, 0.120, 0.124, 0.119],
    ...          [0.116, 0.125, 0.119, 0.125, 0.120]],
    ...         [[0.123, 0.125, 0.125, 0.124, 0.126],
    ...          [0.122, 0.121, 0.124, 0.126, 0.125]],
    ...         [[0.126, 0.128, 0.126, 0.127, 0.129],
    ...          [0.126, 0.129, 0.125, 0.130, 0.124]],
    ...         [[0.118, 0.129, 0.127, 0.120, 0.121],
    ...          [0.125, 0.123, 0.114, 0.124, 0.117]]]

    Note that each horizontal row of data corresponds to a specific
    machine and operator, and comes from a column in the table.

    >>> from scipy.stats import anova
    >>> aov = anova.twoway_nested_from_data_grid(data)

    We can print `aov` to see the full ANOVA table.

    >>> print(aov)
    ANOVA two-way, factor B nested in factor A, B is random
    Source              SS DF        MS          F          p
    Factor A    0.00030332  4 7.583e-05 20.3844086 0.00269263
    Factor B(A)   1.86e-05  5  3.72e-06  0.4300578 0.82493038
    Error         0.000346 40  8.65e-06
    Total       0.00066792 49

    We can access attributes of `aov` for specific values, such as
    the p-values.

    >>> aov.pA
    0.002692630233071879
    >>> aov.pBA
    0.8249303776476935

    If our criterion for statisical significance is a p-value
    of 0.05 or less, than we conclude that the effect of the
    operator is not significant, but that of the machine is.

    """
    # In the following, the letters A and B in variable names
    # refer to the first and second factor, respectively.

    if btype not in ['random', 'fixed']:
        raise ValueError('btype must be "random" or "fixed".')

    data = np.asarray(data)
    shp = data.shape
    if data.ndim == 2:
        raise ValueError("ndim = 2 not implemented yet, use anova.twoway1 "
                         "instead.")
    a, b, s = shp

    grand_mean = data.mean()
    group_means = data.mean(axis=(1, 2), keepdims=True)
    subgroup_means = data.mean(axis=2, keepdims=True)

    SS_among_groups = b * s * ((group_means - grand_mean)**2).sum()
    df_among_groups = a - 1
    MS_among_groups = SS_among_groups / df_among_groups

    SS_subgroups_within_groups = s * ((subgroup_means - group_means)**2).sum()
    df_subgroups_within_groups = a * (b - 1)
    MS_subgroups_within_groups = (SS_subgroups_within_groups
                                  / df_subgroups_within_groups)

    SS_within_groups = ((data - subgroup_means)**2).sum()
    df_within_groups = a * b * (s - 1)
    MS_within_groups = SS_within_groups / df_within_groups

    if btype == 'random':
        F_among_groups_over_subgroups = (MS_among_groups
                                         / MS_subgroups_within_groups)
        DFAdenom = df_subgroups_within_groups
    else:
        # btype is 'fixed'
        F_among_groups_over_subgroups = (MS_among_groups
                                         / MS_within_groups)
        DFAdenom = df_within_groups

    p_among_groups_over_subgroups = fdtrc(df_among_groups,
                                          DFAdenom,
                                          F_among_groups_over_subgroups)

    F_subgroups_over_within_groups = (MS_subgroups_within_groups
                                      / MS_within_groups)

    p_subgroups_over_within_groups = fdtrc(df_subgroups_within_groups,
                                           df_within_groups,
                                           F_subgroups_over_within_groups)

    result = AnovaTwowayNestedResult(SSA=SS_among_groups,
                                     SSBA=SS_subgroups_within_groups,
                                     SSerror=SS_within_groups,
                                     DFA=df_among_groups,
                                     DFBA=df_subgroups_within_groups,
                                     DFerror=df_within_groups,
                                     MSA=MS_among_groups,
                                     MSBA=MS_subgroups_within_groups,
                                     MSerror=MS_within_groups,
                                     FA=F_among_groups_over_subgroups,
                                     FBA=F_subgroups_over_within_groups,
                                     pA=p_among_groups_over_subgroups,
                                     pBA=p_subgroups_over_within_groups,
                                     btype=btype)
    return result


def twoway_nested_from_a_b_values(a, b, values, sstype=None, btype='fixed'):
    """
    Nested two-way ANOVA with factor B nested in factor A.

    The number of replications for each pair of factors must be the same.

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
    btype : str, optional
        Effect type of factor B. Must be one of ``'random'`` or ``'fixed'``.
        Default is ``'fixed'``.

    Returns
    -------
    aov : `AnovaTwoWayNestedResult`
        An object whose attributes contain the information normally
        presented in an ANOVA table.  When the object is printed, the
        output is formatted as an ANOVA table.

    Notes
    -----
    .. versionadded:: 1.7.0

    Examples
    --------

    """
    if btype not in ['random', 'fixed']:
        raise ValueError('`btype` must be "random" or "fixed".')

    a = np.array(a)
    b = np.array(b)

    n = len(a)

    A, Alevels, _ = _encode(a)
    BinA, cell_sizes = _encode_nested(a, b)
    # print("cell_sizes:")
    # print(cell_sizes)
    # row0_sizes = cell_sizes[0]
    # for row_sizes in cell_sizes[1:]:
    #     if not np.array_equal(row_sizes, row0_sizes):
    #         raise ValueError('balanced data required')

    I = np.eye(n)
    c = np.ones(n)

    # Full model
    X1 = np.column_stack((c, A, BinA))
    P1 = _proj(X1)

    X2 = np.column_stack((c, BinA))
    P2 = _proj(X2)

    X3 = np.column_stack((c, A))
    P3 = _proj(X3)

    SSerror = values @ (I - P1) @ values
    SSa = values @ (P1 - P2) @ values
    SSba = values @ (P1 - P3) @ values

    dfa = A.shape[1]
    dfba = BinA.shape[1]
    dferror = sum(cs.sum() - len(cs) for cs in cell_sizes)

    msa = SSa/dfa
    msba = SSba/dfba
    mserror = SSerror/dferror

    if btype == 'random':
        FA = msa / msba
        DFAdenom = dfba
    else:
        # btype is 'fixed'
        FA = msa / mserror
        DFAdenom = dferror

    FBA = msba / mserror

    pA = fdtrc(dfa, DFAdenom, FA)
    pBA = fdtrc(dfba, dferror, FBA)

    result = AnovaTwowayNestedResult(SSA=SSa,
                                     SSBA=SSba,
                                     SSerror=SSerror,
                                     DFA=dfa,
                                     DFBA=dfba,
                                     DFerror=dferror,
                                     MSA=msa,
                                     MSBA=msba,
                                     MSerror=mserror,
                                     FA=FA,
                                     FBA=FBA,
                                     pA=pA,
                                     pBA=pBA,
                                     btype=btype)

    return result
