import numpy as np

#pythran export _mood_inner_lc(int, int, int[], int[], float[], int)

def _mood_inner_lc(js, k, s_lower, s_upper, t, N):
    pass
    # js, k, s_lower, s_upper, t, N = 1, 2, [4,5],[4,5],  [4,5] ,5
    psi = lambda I: (I - (N + 1)/2)**2
    #
    phi_I = [np.arange(s_lower[idx], s_upper[idx]) for idx in range(k)]
    # # for every range in the above array, determine the sum of psi(I) for
    # # every element in the range. Divide all the sums by `t`. Following the
    # # last unnumbered equation on page 312.
    # # 132 µs
    phis = []
    for i in range(len(phi_I)):
        phis.append(sum(psi(phi_I[i])) / t[i + 1])

    # phis = [sum(psi(I_i)) for I_i in phi_I] / t[js]
    return phis
