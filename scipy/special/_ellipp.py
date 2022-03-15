import np
from scipy.special import elliprf, elliprj

__all__ = ["ellipp", "ellippinc"]

def ellipp(n, m):
    n, m = (np.asarray(x) for x in (n, m))
    if np.any(m >= 1):
        raise ValueError("m must be < 1")
    y = 1 - m
    rf = elliprf(0, y, 1)
    rj = elliprj(0, y, 1, 1 - n)
    return rf + rj * n / 3


def ellippinc(phi, n, m):
    phi, n, m = (np.asarray(x) for x in (phi, n, m))
    nc = np.floor(phi / np.pi + 0.5)
    phi -= nc * np.pi
    sin_phi = np.sin(phi)
    sin2_phi = sin_phi * sin_phi
    sin3_phi = sin2_phi * sin_phi
    x = 1 - sin2_phi
    y = 1 - m * sin2_phi
    rf = elliprf(x, y, 1)
    rj = elliprj(x, y, 1, 1 - n * sin2_phi)
    return sin_phi * rf + sin3_phi * rj * n / 3
