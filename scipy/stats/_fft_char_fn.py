import numpy as np
from scipy.integrate._quadrature import _builtincoeffs

# cotes numbers - see sequence from http://oeis.org/A100642
C_table = [[], [1]] + [v[2] for v in _builtincoeffs.values()]
C_npa = np.array(
    [
        np.pad(r, (0, len(C_table) - 1 - len(r)), mode='constant')
        for r in C_table
    ]
)


def pdf_from_cf_with_fft(cf, h=0.01, q=9, level=3):
    """Calculates pdf from cf using fft using Simpsons as suggest by [WZ]
    when level=3. Also can calculate pdf using higher order cote rules that
    provide greater accuracy.
    """
    n = level
    N = 2**q
    j = l = np.arange(0, N)
    L = N * h / 2
    x_l = np.pi * (l - N / 2) / L
    if level > 1:
        k = np.arange(n).reshape(n, 1)
        s1 = np.sum(
            (-1)**l * C_npa[n, k] * np.fft.fft(
                (-1)**j * cf(-L + h * j + h * k / (n - 1))
            ) * np.exp(
                1j * np.pi * k / (n - 1) - 2 * 1j * np.pi * k * l /
                (N * (n - 1))
            ),
            axis=0
        )
    else:
        s1 = (-1)**l * C_npa[n, 0] * np.fft.fft((-1)**j * cf(-L + h * j))
    density = h * s1 / (2 * np.pi * np.sum(C_npa[n]))
    return (x_l, density)
