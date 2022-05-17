import numpy as np

#pythran export _mood_inner_lc(float[], float[])

def _mood_inner_lc(y, x) -> float:
    n = x.shape[0]
    m = y.shape[0]
    N = n + m
    xy = np.concatenate((x, y))

    # obtain the unique values and the counts of each.
    # "a_i, + b_i, = t_i, for j = 1, ... k", where `k` is the number of unique
    # classes, and "[t]he number of values associated with the x's and y's in
    # the jth class will be denoted by a_i, and b_i respectively."
    # (Mielke, 312)
    uniques, t = np.unique(sorted(xy), return_counts=1)
    k = len(uniques)
    js = list(np.arange(1, k + 1))

    # the `b` array mentioned in the paper is not used, outside of the
    # calculation of `t`, so we do not need to calculate it separately. Here
    # we calculate `a`. In plain language, `a[i]` is the number of values in
    # `x` that equal `uniques[i]`.
    _, xyx_counts = np.unique(sorted(np.concatenate((xy, x))),
                              return_counts=1)
    a = xyx_counts - t
    # "Define .. a_0 = b_0 = t_0 = S_0 = 0" (Mielke 312) so we shift  `a`
    # and `t` arrays over 1 to allow a first element of 0 to accommodate this
    # indexing.
    t = np.concatenate(([0], t))
    a = np.concatenate(([0], a))

    # S is built from `t`, so it does not need a preceding zero added on.
    S = np.cumsum(t)
    # define a copy of `S` with a prepending zero for later use to avoid
    # indexing
    S_i_m1 = np.concatenate(([0], S[:-1]))

    # Psi, as defined by the 6th unnumbered equation on page 312 (Mielke).
    # Note that in the paper there is an error where the denominator `2` is
    # squared when it should be the entire equation.
    psi = lambda I: (I - (N + 1)/2)**2

    # define summation range for use in calculation of phi, as seen in sum
    # in the unnumbered equation on the bottom of page 312 (Mielke).
    s_lower = [S[jsi - 1] + 1 for jsi in js]
    s_upper = [S[jsi] + 1 for jsi in js]
    phi_I = [range(s_lower[idx], s_upper[idx]) for idx in range(k)]

    # for every range in the above array, determine the sum of psi(I) for
    # every element in the range. Divide all the sums by `t`. Following the
    # last unnumbered equation on page 312.
    phis = []
    for i in range(len(phi_I)):
        psi_temp = [psi(phi_I_i) for phi_I_i in phi_I[i]]
        phis.append(sum(psi_temp) / t[(i + 1)])

    # `T` is equal to a[j] * phi[j], per the first unnumbered equation on
    # page 312. `phis` is already in the order based on `js`, so we index
    # into `a` with `js` as well.
    T = sum(phis * a[js])

    # The approximate statistic
    E_0_T = n * (N * N - 1) / 12

    varM = (m * n * (N + 1.0) * (N ** 2 - 4) / 180 -
            m * n / (180 * N * (N - 1)) * np.sum(
                t * (t**2 - 1) * (t**2 - 4 + (15 * (N - S - S_i_m1) ** 2))
            ))

    return (T - E_0_T) / np.sqrt(varM)
