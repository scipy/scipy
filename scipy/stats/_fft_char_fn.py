import numpy as np
from scipy.integrate._quadrature import _builtincoeffs

# cotes numbers - see sequence from http://oeis.org/A100642
Cotes_table = [[], [1]] + [v[2] for v in _builtincoeffs.values()]
Cotes = np.array(
    [
        np.pad(r, (0, len(Cotes_table) - 1 - len(r)), mode='constant')
        for r in Cotes_table
    ]
)


def pdf_from_cf_with_fft(cf, h=0.01, q=9, level=3):
    """Calculates pdf from characteristic function.

    Uses fast Fourier transform with Newton-Cotes integration following [WZ].
    Defaults to using Simpson's method (3-point Newton-Cotes integration).

    Parameters
    ----------
    cf : callable
        Single argument function from float -> complex expressing a
        characteristic function for some distribution.
    h : Optional[float]
        Step size for Newton-Cotes integration. Default: 0.01
    q : Optional[int]
        Use 2**q steps when peforming Newton-Cotes integration.
        The infinite integral in the inverse Fourier transform will then
        be restricted to the interval [-2**q * h / 2, 2**q * h / 2]. Setting
        the number of steps equal to a power of 2 allows the fft to be
        calculated in O(n*log(n)) time rather than O(n**2).
        Default: 9
    level : Optional[int]
        Calculate integral using n-point Newton-Cotes integration for
        n = level. The 3-point Newton-Cotes formula corresponds to Simpson's
        rule. Default: 3

    Returns
    -------
    x_l : ndarray
        Array of points x at which pdf is estimated. 2**q equally spaced
        points from -pi/h up to but not including pi/h.
    density : ndarray
        Estimated values of pdf corresponding to cf at points in x_l.

    References
    ----------
    .. [WZ] Wang, Li and Zhang, Ji-Hong, 2008. Simpson's rule based FFT method
        to compute densities of stable distribution.
    """
    n = level
    N = 2**q
    steps = np.arange(0, N)
    L = N * h / 2
    x_l = np.pi * (steps - N / 2) / L
    if level > 1:
        indices = np.arange(n).reshape(n, 1)
        s1 = np.sum(
            (-1) ** steps * Cotes[n, indices] * np.fft.fft(
                (-1)**steps * cf(-L + h * steps + h * indices / (n - 1))
            ) * np.exp(
                1j * np.pi * indices / (n - 1)
                - 2 * 1j * np.pi * indices * steps /
                (N * (n - 1))
            ),
            axis=0
        )
    else:
        s1 = (-1) ** steps * Cotes[n, 0] * np.fft.fft(
            (-1) ** steps * cf(-L + h * steps)
        )
    density = h * s1 / (2 * np.pi * np.sum(Cotes[n]))
    return (x_l, density)
