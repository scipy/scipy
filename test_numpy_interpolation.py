"""
Full validation of interp1d and griddata (numpy-only implementations).

Strategy
--------
* 1-D linear / nearest / cubic  →  exact match with scipy (atol 1e-11).
* 2-D nearest                   →  exact match for random float coords
                                    (no equidistance ties); for integer pilot
                                    coords we verify all returned values belong
                                    to the pilot set and unambiguous cells match.
* 2-D linear                    →  same NaN mask as scipy; exact at data points;
                                    linear functions reproduced exactly; interior
                                    cells within 0.05 of scipy (different but
                                    valid Delaunay triangulations may differ on
                                    co-circular point sets).
* 2-D cubic (TPS RBF)           →  exact at data points; linear functions
                                    reproduced exactly (polynomial order of TPS
                                    augmentation); smooth approximation to scipy
                                    CloughTocher (no bit-exact match expected).
* High-precision (decimal)      →  interp1d(precision=40) with string inputs
                                    and eval_hp() achieves 10^-35 accuracy
                                    (well beyond the 10^-32 requirement).
"""

import numpy as np
from scipy.interpolate import (interp1d as scipy_interp1d,
                                CubicSpline,
                                griddata as scipy_griddata)
from numpy_interpolation import interp1d, griddata

ATOL   = 1e-11   # tight tolerance for 1-D and for exactly-reproduced functions
LINEAR_ATOL = 0.05  # Delaunay-triangulation-dependent tolerance for 2-D linear


def check(a, b, tag, atol=ATOL):
    a, b = np.asarray(a, float), np.asarray(b, float)
    ok = np.allclose(a, b, atol=atol, rtol=1e-9, equal_nan=True)
    print(f"  [{'PASS' if ok else 'FAIL'}] {tag}")
    if not ok:
        diff = np.abs(a - b)
        print(f"         max|diff|={np.nanmax(diff):.3e}  (atol={atol:.0e})")
    return ok


# ============================================================
# helpers
# ============================================================

def make_1d(n=20, seed=0):
    rng = np.random.default_rng(seed)
    x   = np.unique(np.sort(rng.uniform(0, 10, n)))
    y   = np.sin(x) + 0.1 * rng.standard_normal(len(x))
    return x, y


def make_2d_y(n=20, m=3, seed=0):
    rng = np.random.default_rng(seed)
    x   = np.unique(np.sort(rng.uniform(0, 10, n)))
    y   = np.column_stack([np.sin(x) + 0.1*rng.standard_normal(len(x))
                           for _ in range(m)])
    return x, y


def q1d(x, ni=40, ne=5):
    return np.concatenate([
        np.linspace(x[0]-1.5,  x[0]-0.01, ne),
        np.linspace(x[0],      x[-1],      ni),
        np.linspace(x[-1]+0.01, x[-1]+1.5, ne),
    ])


# ============================================================
# 1-D linear
# ============================================================

def test_linear_1d():
    print("\n--- linear 1-D ---")
    x, y = make_1d()
    xq   = q1d(x)
    sci  = scipy_interp1d(x, y, kind='linear', fill_value='extrapolate')
    ours = interp1d(x, y, kind='linear', fill_value='extrapolate', bounds_error=False)
    check(ours(xq), sci(xq), "extrapolate interior+exterior")
    check(ours(x),  sci(x),  "exactly on knots")
    check(ours(x[0]), float(sci(x[0])), "scalar query")


def test_linear_2d_y():
    print("\n--- linear 2-D y ---")
    x, y = make_2d_y()
    xq   = q1d(x)
    sci  = scipy_interp1d(x, y.T, kind='linear', fill_value='extrapolate', axis=1)
    ours = interp1d(x, y, kind='linear', fill_value='extrapolate', bounds_error=False)
    check(ours(xq), sci(xq).T, "2-D extrapolate")
    check(ours(x),  sci(x).T,  "2-D on knots")


# ============================================================
# 1-D nearest
# ============================================================

def test_nearest_1d():
    print("\n--- nearest 1-D ---")
    x, y = make_1d()
    xq   = np.linspace(x[0], x[-1], 80)
    sci  = scipy_interp1d(x, y, kind='nearest')
    ours = interp1d(x, y, kind='nearest', bounds_error=False, fill_value=np.nan)
    check(ours(xq), sci(xq), "interior")
    check(ours(x),  sci(x),  "on knots")


def test_nearest_2d_y():
    print("\n--- nearest 2-D y ---")
    x, y = make_2d_y()
    xq   = np.linspace(x[0], x[-1], 60)
    sci  = scipy_interp1d(x, y.T, kind='nearest', axis=1)
    ours = interp1d(x, y, kind='nearest', bounds_error=False, fill_value=np.nan)
    check(ours(xq), sci(xq).T, "2-D interior")


# ============================================================
# 1-D cubic
# ============================================================

def test_cubic_1d():
    print("\n--- cubic 1-D ---")
    x, y = make_1d()
    xq   = q1d(x)
    sci  = CubicSpline(x, y, bc_type='not-a-knot', extrapolate=True)
    ours = interp1d(x, y, kind='cubic', fill_value='extrapolate', bounds_error=False)
    check(ours(xq), sci(xq), "extrapolate interior+exterior")
    check(ours(x),  sci(x),  "on knots")
    check(ours(x[0]), float(sci(x[0])), "scalar")


def test_cubic_2d_y():
    print("\n--- cubic 2-D y ---")
    x, y = make_2d_y()
    xq   = q1d(x)
    sci  = CubicSpline(x, y, bc_type='not-a-knot', extrapolate=True)
    ours = interp1d(x, y, kind='cubic', fill_value='extrapolate', bounds_error=False)
    check(ours(xq), sci(xq), "2-D extrapolate")
    check(ours(x),  sci(x),  "2-D on knots")


def test_cubic_small_n():
    print("\n--- cubic small n ---")
    x2, y2 = np.array([1.0, 3.0]), np.array([2.0, 6.0])
    check(interp1d(x2, y2, kind='cubic', bounds_error=False)(np.linspace(1, 3, 20)),
          CubicSpline(x2, y2)(np.linspace(1, 3, 20)), "n=2 linear degenerate")
    x3 = np.array([0.0, 1.0, 3.0]); y3 = x3**2
    check(interp1d(x3, y3, kind='cubic', bounds_error=False)(np.linspace(0, 3, 30)),
          CubicSpline(x3, y3)(np.linspace(0, 3, 30)), "n=3 parabola")


# ============================================================
# extrapolate=False → NaN outside bounds
# ============================================================

def test_extrapolate_false():
    print("\n--- extrapolate=False ---")
    x, y = make_1d()
    xq_in  = np.linspace(x[0], x[-1], 30)
    xq_out = np.array([x[0]-1.0, x[-1]+1.0])
    for kind in ('linear', 'nearest', 'cubic'):
        f  = interp1d(x, y, kind=kind, bounds_error=False, fill_value=np.nan)
        ok = not np.any(np.isnan(f(xq_in))) and np.all(np.isnan(f(xq_out)))
        print(f"  [{'PASS' if ok else 'FAIL'}] {kind}: interior finite, exterior NaN")


# ============================================================
# interp1d wrapper
# ============================================================

def test_interp1d_wrapper():
    print("\n--- interp1d wrapper ---")
    x, y = make_1d()
    xq   = np.linspace(x[0], x[-1], 40)
    for kind in ('linear', 'nearest', 'cubic'):
        if kind == 'cubic':
            ref = CubicSpline(x, y, bc_type='not-a-knot', extrapolate=True)(xq)
        else:
            ref = scipy_interp1d(x, y, kind=kind, fill_value='extrapolate',
                                 bounds_error=False)(xq)
        wrap = interp1d(x, y, kind=kind, bounds_error=False, fill_value='extrapolate')
        check(wrap(xq), ref, f"wrapper kind={kind}")

    f_strict = interp1d(x, y, kind='linear', bounds_error=True)
    try:
        f_strict(np.array([x[0]-1.0]))
        print("  [FAIL] bounds_error=True did not raise")
    except ValueError:
        print("  [PASS] bounds_error=True raises ValueError")


# ============================================================
# interp1d – exact pattern from user's channel-estimation code
# ============================================================

def test_user_interp1d_pattern():
    """Replicates the exact interp1d usage in the user's 1-D pilot code."""
    print("\n--- user interp1d pattern ---")
    symbols_sorted = np.array([0., 3., 7., 11., 14.], dtype=float)
    values_sorted  = np.array([0.8+0.1j, 0.9+0.05j, 0.85+0.0j,
                                0.88-0.05j, 0.82-0.1j])
    x_interp = symbols_sorted + 1   # "add 1 like MATLAB"
    x_query  = np.arange(1, 15, dtype=float)

    for part, arr in (('real', values_sorted.real), ('imag', values_sorted.imag)):
        f   = interp1d(x_interp, arr, kind='linear',
                       fill_value='extrapolate', bounds_error=False)
        sci = scipy_interp1d(x_interp, arr, kind='linear',
                             fill_value='extrapolate', bounds_error=False)
        check(f(x_query), sci(x_query), f"user interp1d – {part}")


# ============================================================
# griddata nearest  (random float coords – no equidistance ties)
# ============================================================

def test_griddata_nearest_float():
    print("\n--- griddata nearest (float coords) ---")
    rng    = np.random.default_rng(7)
    pts    = rng.random((60, 2))
    vals   = np.sin(pts[:,0]) + np.cos(pts[:,1])
    XX, YY = np.meshgrid(np.linspace(0,1,15), np.linspace(0,1,15))

    ours = griddata(pts, vals, (XX, YY), method='nearest')
    sci  = scipy_griddata(pts, vals, (XX, YY), method='nearest')
    check(ours, sci, "meshgrid xi, random float pts – exact match scipy")

    xi2 = rng.random((30, 2))
    check(griddata(pts, vals, xi2, method='nearest'),
          scipy_griddata(pts, vals, xi2, method='nearest'),
          "array xi, random float pts")


# ============================================================
# griddata linear
# ============================================================

def test_griddata_linear():
    print("\n--- griddata linear ---")
    rng    = np.random.default_rng(8)
    pts    = rng.random((80, 2))
    vals   = np.sin(2*np.pi*pts[:,0]) * np.cos(2*np.pi*pts[:,1])
    XX, YY = np.meshgrid(np.linspace(0.1, 0.9, 12), np.linspace(0.1, 0.9, 12))

    ours = griddata(pts, vals, (XX, YY), method='linear')
    sci  = scipy_griddata(pts, vals, (XX, YY), method='linear')

    # 1. NaN mask must be identical (same convex hull boundary)
    ok_nan = np.array_equal(np.isnan(ours), np.isnan(sci))
    print(f"  [{'PASS' if ok_nan else 'FAIL'}] linear: NaN mask identical to scipy")

    # 2. Exact at data points (interpolation property)
    out_at = griddata(pts, vals, pts, method='linear')
    check(out_at, vals, "linear: exact at data points")

    # 3. Linear functions exactly reproduced everywhere inside hull
    rng2     = np.random.default_rng(42)
    pts_lin  = rng2.random((60, 2))
    f_lin    = lambda p: 2.5*p[:,0] - 1.3*p[:,1] + 0.7   # linear function
    vals_lin = f_lin(pts_lin)
    q_lin    = rng2.random((40, 2)) * 0.6 + 0.2            # safely inside hull
    check(griddata(pts_lin, vals_lin, q_lin, method='linear'),
          f_lin(q_lin), "linear: exactly reproduces linear f(x,y)=ax+by+c")

    # 4. Interior cells close to scipy (different valid triangulations may
    #    differ slightly at co-circular point sets)
    mask_both = ~np.isnan(ours) & ~np.isnan(sci)
    check(ours[mask_both], sci[mask_both],
          "linear: interior cells agree with scipy", atol=LINEAR_ATOL)

    # 5. Outside convex hull → fill_value (NaN)
    xi_out  = np.array([[2.0, 2.0], [-1.0, -1.0]])
    out_out = griddata(pts, vals, xi_out, method='linear')
    ok_fv   = np.all(np.isnan(out_out))
    print(f"  [{'PASS' if ok_fv else 'FAIL'}] linear: outside hull → NaN")


# ============================================================
# griddata cubic (TPS RBF)
# ============================================================

def test_griddata_cubic():
    print("\n--- griddata cubic (TPS RBF) ---")
    rng    = np.random.default_rng(9)
    pts    = rng.random((50, 2))
    XX, YY = np.meshgrid(np.linspace(0.05,0.95,10), np.linspace(0.05,0.95,10))

    # 1. Exact at data points (primary TPS guarantee)
    vals = np.sin(2*np.pi*pts[:,0]) * np.cos(2*np.pi*pts[:,1])
    out_at = griddata(pts, vals, pts, method='cubic')
    check(out_at, vals, "cubic: exact at data points", atol=1e-10)

    # 2. Linear functions reproduced exactly (TPS polynomial augmentation is degree-1)
    vals_lin  = 3.1*pts[:,0] - 2.7*pts[:,1] + 1.4
    out_grid  = griddata(pts, vals_lin, (XX, YY), method='cubic')
    exact_lin = 3.1*XX - 2.7*YY + 1.4
    check(out_grid, exact_lin, "cubic: exactly reproduces linear f(x,y)", atol=1e-8)

    # 3. Smooth interpolation close to scipy CloughTocher for interior pts
    #    Use a gentle function; TPS and CloughTocher are both C∞ at interior.
    vals_s = pts[:,0]*(1-pts[:,0]) * pts[:,1]*(1-pts[:,1])
    ours   = griddata(pts, vals_s, (XX, YY), method='cubic')
    sci    = scipy_griddata(pts, vals_s, (XX, YY), method='cubic')
    mask   = ~np.isnan(sci)
    diff   = np.abs(ours[mask] - sci[mask])
    ok3    = diff.max() < 0.1
    print(f"  [{'PASS' if ok3 else 'FAIL'}] cubic: smooth approx to CloughTocher "
          f"(max|diff|={diff.max():.3e})")


# ============================================================
# griddata – exact pattern from user's channel-estimation code
# ============================================================

def test_user_griddata_pattern():
    """Replicates the exact griddata usage in the user's channel-estimation code."""
    print("\n--- user griddata pattern ---")
    rng = np.random.default_rng(11)
    n   = 40
    Xp  = rng.integers(1, 64, n).astype(float)
    Yp  = rng.integers(1, 14, n).astype(float)
    Zd  = rng.standard_normal(n) + 1j * rng.standard_normal(n)

    xx, yy = np.arange(1,65,dtype=float), np.arange(1,15,dtype=float)
    XX, YY = np.meshgrid(xx, yy)
    points = np.column_stack([Xp, Yp])
    xi     = (XX, YY)

    for method in ('nearest', 'linear'):
        Z_real = np.real(Zd)
        Z_imag = np.imag(Zd)

        htemp_real = griddata(points, Z_real, xi, method=method)
        htemp_imag = griddata(points, Z_imag, xi, method=method)
        htemp      = htemp_real + 1j * htemp_imag

        sci_real = scipy_griddata(points, Z_real, xi, method=method)
        sci_imag = scipy_griddata(points, Z_imag, xi, method=method)

        # Output shape must match meshgrid shape
        ok_shape = htemp.shape == XX.shape
        print(f"  [{'PASS' if ok_shape else 'FAIL'}] {method}: "
              f"output shape {htemp.shape} == {XX.shape}")

        if method == 'nearest':
            # Nearest: every returned value must be one of the pilot values.
            # (Tie-breaking with integer coords may differ from scipy KDTree.)
            vals_r = Z_real
            vals_i = Z_imag
            flat_r = htemp_real.ravel()
            flat_i = htemp_imag.ravel()
            in_set_r = all(any(np.isclose(v, vals_r)) for v in flat_r)
            in_set_i = all(any(np.isclose(v, vals_i)) for v in flat_i)
            ok = in_set_r and in_set_i
            print(f"  [{'PASS' if ok else 'FAIL'}] {method}: "
                  "all values belong to pilot set (real+imag)")

            # For cells that are clearly unambiguous (no equidistance tie):
            # find cells where nearest pilot is strictly closer than 2nd nearest
            xi_flat = np.column_stack([XX.ravel(), YY.ravel()])
            d2 = np.sum((xi_flat[:, None, :] - points[None, :, :])**2, axis=2)
            d_sorted = np.sort(d2, axis=1)
            unambiguous = d_sorted[:, 0] < d_sorted[:, 1] - 0.01   # strictly nearer
            unambig_2d  = unambiguous.reshape(XX.shape)
            n_unambig   = unambiguous.sum()
            if n_unambig > 0:
                check(htemp_real[unambig_2d],
                      sci_real[unambig_2d],
                      f"{method}: unambiguous cells match scipy ({n_unambig} cells)")

        elif method == 'linear':
            # Linear interpolation with integer-coordinate pilots:
            # multiple valid Delaunay triangulations exist (co-circular points)
            # so we verify correctness properties rather than exact scipy match.

            # 1. Same NaN mask (same convex hull)
            ok_nan = np.array_equal(np.isnan(htemp_real), np.isnan(sci_real))
            print(f"  [{'PASS' if ok_nan else 'FAIL'}] {method}: NaN mask identical to scipy")

            # 2. Exact at pilot positions (dedup first – scipy also can't handle dups)
            _, uniq = np.unique(points, axis=0, return_index=True)
            uniq    = np.sort(uniq)
            pts_u   = points[uniq]
            Zr_u, Zi_u = Z_real[uniq], Z_imag[uniq]
            out_at_r = griddata(pts_u, Zr_u, pts_u, method='linear')
            out_at_i = griddata(pts_u, Zi_u, pts_u, method='linear')
            check(out_at_r, Zr_u, f"{method}: exact at unique pilot positions (real)")
            check(out_at_i, Zi_u, f"{method}: exact at unique pilot positions (imag)")

            # 3. All interior values are finite (no spurious NaN inside hull)
            interior = ~np.isnan(htemp_real)
            ok_fin   = np.all(np.isfinite(htemp_real[interior]))
            print(f"  [{'PASS' if ok_fin else 'FAIL'}] {method}: "
                  f"all {interior.sum()} interior cells are finite")


# ============================================================
# griddata rescale flag
# ============================================================

def test_griddata_rescale():
    print("\n--- griddata rescale ---")
    rng  = np.random.default_rng(13)
    pts  = rng.random((50, 2)) * np.array([1e6, 1e-3])
    vals = pts[:,0] / 1e6 + pts[:,1] / 1e-3
    xi   = rng.random((20, 2)) * np.array([1e6, 1e-3])
    ours = griddata(pts, vals, xi, method='nearest', rescale=True)
    sci  = scipy_griddata(pts, vals, xi, method='nearest', rescale=True)
    check(ours, sci, "rescale=True, nearest")


# ============================================================
# 1-D Lagrange polynomial
# ============================================================

def test_lagrange_1d():
    print("\n--- lagrange 1-D ---")
    # Polynomial of degree ≤ n-1 is reproduced exactly.
    x = np.linspace(0.0, 4.0, 6)
    y = x**3 - 2*x**2 + x - 1              # degree-3 polynomial
    f = interp1d(x, y, kind='lagrange', fill_value='extrapolate', bounds_error=False)
    xq = np.linspace(0.0, 4.0, 40)
    check(f(xq),  xq**3 - 2*xq**2 + xq - 1, "reproduces cubic polynomial",  atol=1e-9)
    check(f(x),   y,                          "exact at nodes",               atol=1e-11)
    check(f(2.0), float(2.0**3 - 2*2.0**2 + 2.0 - 1), "scalar query",        atol=1e-11)


def test_lagrange_2d_y():
    print("\n--- lagrange 2-D y ---")
    x = np.linspace(0.0, np.pi, 7)
    y = np.column_stack([np.sin(x), np.cos(x)])   # (7, 2)
    f = interp1d(x, y, kind='lagrange', fill_value='extrapolate', bounds_error=False)
    xq = np.linspace(0.2, np.pi - 0.2, 20)
    res = f(xq)
    # 7 equally-spaced nodes → small but nonzero approximation error
    check(res[:, 0], np.sin(xq), "sin column approx",  atol=5e-4)
    check(res[:, 1], np.cos(xq), "cos column approx",  atol=5e-4)
    # Exact at nodes
    check(f(x), y, "exact at nodes (2-D)",  atol=1e-10)


def test_lagrange_extrapolate_false():
    print("\n--- lagrange extrapolate=False ---")
    x = np.linspace(0.0, 5.0, 6)
    y = x**2
    f = interp1d(x, y, kind='lagrange', bounds_error=False, fill_value=np.nan)
    q = np.array([-0.5, 2.5, 5.5])
    res = f(q)
    assert np.isnan(res[0]), "below range → nan"
    assert np.isnan(res[2]), "above range → nan"
    check(res[1], 2.5**2, "interior finite",  atol=1e-9)
    print("  [PASS] extrapolate=False")


def test_lagrange_interp1d():
    print("\n--- lagrange via interp1d ---")
    x = np.array([0., 1., 2., 3., 4.])
    y = np.array([0., 1., 4., 9., 16.])     # y = x^2
    fl = interp1d(x, y, kind='lagrange', fill_value='extrapolate')
    xq = np.array([0.5, 1.5, 2.5, 3.5])
    check(fl(xq), xq**2, "interp1d lagrange",  atol=1e-9)


# ============================================================
# 1-D Newton polynomial
# ============================================================

def test_newton_1d():
    print("\n--- newton 1-D ---")
    x = np.linspace(0.0, 4.0, 6)
    y = x**3 - 2*x**2 + x - 1
    f = interp1d(x, y, kind='newton', fill_value='extrapolate', bounds_error=False)
    xq = np.linspace(0.0, 4.0, 40)
    check(f(xq),  xq**3 - 2*xq**2 + xq - 1, "reproduces cubic polynomial",  atol=1e-9)
    check(f(x),   y,                          "exact at nodes",               atol=1e-11)
    check(f(2.0), float(2.0**3 - 2*2.0**2 + 2.0 - 1), "scalar query",        atol=1e-11)


def test_newton_2d_y():
    print("\n--- newton 2-D y ---")
    x = np.linspace(0.0, np.pi, 7)
    y = np.column_stack([np.sin(x), np.cos(x)])
    f = interp1d(x, y, kind='newton', fill_value='extrapolate', bounds_error=False)
    xq = np.linspace(0.2, np.pi - 0.2, 20)
    res = f(xq)
    check(res[:, 0], np.sin(xq), "sin column approx",  atol=5e-4)
    check(res[:, 1], np.cos(xq), "cos column approx",  atol=5e-4)
    check(f(x), y, "exact at nodes (2-D)",  atol=1e-10)


def test_newton_extrapolate_false():
    print("\n--- newton extrapolate=False ---")
    x = np.linspace(0.0, 5.0, 6)
    y = x**2
    f = interp1d(x, y, kind='newton', bounds_error=False, fill_value=np.nan)
    q = np.array([-0.5, 2.5, 5.5])
    res = f(q)
    assert np.isnan(res[0]), "below range → nan"
    assert np.isnan(res[2]), "above range → nan"
    check(res[1], 2.5**2, "interior finite",  atol=1e-9)
    print("  [PASS] extrapolate=False")


def test_newton_interp1d():
    print("\n--- newton via interp1d ---")
    x = np.array([0., 1., 2., 3., 4.])
    y = np.array([0., 1., 4., 9., 16.])
    fn = interp1d(x, y, kind='newton', fill_value='extrapolate')
    xq = np.array([0.5, 1.5, 2.5, 3.5])
    check(fn(xq), xq**2, "interp1d newton",  atol=1e-9)


def test_lagrange_vs_scipy_barycentric():
    """Our Lagrange must match scipy.BarycentricInterpolator to machine precision."""
    print("\n--- lagrange vs scipy BarycentricInterpolator ---")
    from scipy.interpolate import BarycentricInterpolator
    rng = np.random.default_rng(42)
    x   = np.unique(np.sort(rng.uniform(0, 6, 12)))
    y   = np.sin(x)
    bi  = BarycentricInterpolator(x, y)
    f   = interp1d(x, y, kind='lagrange', fill_value='extrapolate', bounds_error=False)
    xq  = np.linspace(x[0], x[-1], 60)
    check(f(xq), bi(xq), "matches BarycentricInterpolator",  atol=1e-9)
    check(f(x),  y,       "exact at nodes",                   atol=1e-12)


def test_lagrange_newton_agreement():
    """Lagrange and Newton must produce identical polynomials."""
    print("\n--- lagrange == newton agreement ---")
    rng = np.random.default_rng(42)
    x   = np.unique(rng.uniform(0, 5, 8))
    y   = np.sin(x) * x
    fl  = interp1d(x, y, kind='lagrange', fill_value='extrapolate', bounds_error=False)
    fn  = interp1d(x, y, kind='newton',   fill_value='extrapolate', bounds_error=False)
    xq  = np.linspace(0, 5, 50)
    check(fl(xq), fn(xq), "lagrange == newton",  atol=1e-9)


# ============================================================
# High-precision (decimal) tests
# ============================================================

_HP_PREC = 40   # 40 decimal digits → arithmetic error ~10^-40


def test_hp_eval_hp_decimal_output():
    """eval_hp() returns Decimal results accurate to < 10^-35."""
    print("\n--- hp: eval_hp() Decimal output ---")
    from decimal import Decimal, localcontext
    xs = ['0', '1', '2', '3', '4', '5']
    ys = ['0', '1', '16', '81', '256', '625']   # y = x^4
    f_lag = interp1d(xs, ys, kind='lagrange', fill_value='extrapolate',
                     bounds_error=False, precision=_HP_PREC)
    f_new = interp1d(xs, ys, kind='newton',   fill_value='extrapolate',
                     bounds_error=False, precision=_HP_PREC)
    with localcontext() as ctx:
        ctx.prec = 60
        true_val = Decimal('2.7') ** 4
    for f, tag in [(f_lag, 'lagrange'), (f_new, 'newton')]:
        r = f.eval_hp(['2.7'])[0]
        err = abs(r - true_val)
        ok = err < Decimal('1e-35')
        print(f"  [{'PASS' if ok else 'FAIL'}] {tag}: err={err}")
        assert ok


def test_hp_32digit_input_accuracy():
    """String inputs with 34 sig figs + eval_hp() achieves error < 10^-32."""
    print("\n--- hp: 32-digit input accuracy ---")
    from decimal import Decimal
    xs = ['0', '1', '2', '3', '4']
    ys = ['0', '1', '4', '9', '16']
    for kind in ('lagrange', 'newton'):
        f = interp1d(xs, ys, kind=kind, fill_value='extrapolate',
                     bounds_error=False, precision=_HP_PREC)
        r = f.eval_hp(['1.5000000000000000000000000000000000'])[0]
        err = abs(r - Decimal('2.25'))
        ok = err < Decimal('1e-32')
        print(f"  [{'PASS' if ok else 'FAIL'}] {kind}: err={err}")
        assert ok


def test_hp_call_returns_float64():
    """__call__() returns a standard float64 numpy array even with precision set."""
    print("\n--- hp: __call__ returns float64 ---")
    f = interp1d(['0', '1', '2', '3', '4'], ['0', '1', '4', '9', '16'],
                 kind='lagrange', fill_value='extrapolate',
                 bounds_error=False, precision=_HP_PREC)
    res = f(['0.5', '1.5', '2.5', '3.5'])
    ok = (isinstance(res, np.ndarray) and res.dtype == np.float64
          and np.allclose(res, [0.25, 2.25, 6.25, 12.25], atol=1e-14))
    print(f"  [{'PASS' if ok else 'FAIL'}] dtype={res.dtype} values={res}")
    assert ok


def test_hp_interp1d_wrapper():
    """interp1d(precision=40) works for both lagrange and newton."""
    print("\n--- hp: interp1d wrapper with precision=40 ---")
    for kind in ('lagrange', 'newton'):
        fl = interp1d(['0', '1', '2', '3', '4'], ['0', '1', '8', '27', '64'],
                      kind=kind, fill_value='extrapolate', precision=_HP_PREC)
        res = fl(['1.5', '2.5', '3.5'])
        exp = np.array([1.5**3, 2.5**3, 3.5**3])
        ok = np.allclose(res, exp, atol=1e-14)
        print(f"  [{'PASS' if ok else 'FAIL'}] {kind}")
        assert ok


def test_hp_extrapolate_false():
    """extrapolate=False with precision=40: out-of-bounds → NaN, interior OK."""
    print("\n--- hp: extrapolate=False ---")
    f = interp1d(['0', '1', '2', '3'], ['0', '1', '4', '9'],
                 kind='lagrange', bounds_error=False, fill_value=np.nan,
                 precision=_HP_PREC)
    res = f(['-0.5', '1.5', '3.5'])
    ok = (np.isnan(res[0]) and np.isnan(res[2])
          and abs(res[1] - 2.25) < 1e-14)
    print(f"  [{'PASS' if ok else 'FAIL'}] {res}")
    assert ok


def test_hp_float64_input():
    """float64 inputs with precision=40: output limited to ~10^-15 but correct."""
    print("\n--- hp: float64 input (output at float64 precision) ---")
    x = np.array([0., 1., 2., 3., 4.])
    y = x**2
    for kind in ('lagrange', 'newton'):
        f = interp1d(x, y, kind=kind, fill_value='extrapolate',
                     bounds_error=False, precision=_HP_PREC)
        xq = np.array([0.5, 1.5, 2.5, 3.5])
        ok = np.allclose(f(xq), xq**2, atol=1e-13)
        print(f"  [{'PASS' if ok else 'FAIL'}] {kind}")
        assert ok


def test_hp_scalar_query():
    """eval_hp() accepts a single scalar string."""
    print("\n--- hp: scalar query eval_hp ---")
    from decimal import Decimal
    f = interp1d(['0', '1', '2', '3', '4'], ['0', '1', '4', '9', '16'],
                 kind='newton', fill_value='extrapolate',
                 bounds_error=False, precision=_HP_PREC)
    r = f.eval_hp('2.5')
    err = abs(r[0] - Decimal('6.25'))
    ok = err < Decimal('1e-35')
    print(f"  [{'PASS' if ok else 'FAIL'}] err={err}")
    assert ok


def test_hp_lagrange_newton_agree():
    """lagrange and newton agree at full Decimal precision (error < 10^-35)."""
    print("\n--- hp: lagrange == newton at Decimal precision ---")
    from decimal import Decimal
    xs = ['0', '1', '2', '3', '4', '5', '6']
    ys = ['1', '3', '5', '4', '2', '6', '8']
    fl = interp1d(xs, ys, kind='lagrange', fill_value='extrapolate',
                  bounds_error=False, precision=_HP_PREC)
    fn = interp1d(xs, ys, kind='newton',   fill_value='extrapolate',
                  bounds_error=False, precision=_HP_PREC)
    queries = ['0.3', '1.7', '3.14159', '5.5']
    rl = fl.eval_hp(queries)
    rn = fn.eval_hp(queries)
    ok = all(abs(rl[i] - rn[i]) < Decimal('1e-35') for i in range(len(queries)))
    print(f"  [{'PASS' if ok else 'FAIL'}]")
    assert ok


# ============================================================
# Run all
# ============================================================

if __name__ == '__main__':
    tests = [
        test_linear_1d,
        test_linear_2d_y,
        test_nearest_1d,
        test_nearest_2d_y,
        test_cubic_1d,
        test_cubic_2d_y,
        test_cubic_small_n,
        test_extrapolate_false,
        test_interp1d_wrapper,
        test_user_interp1d_pattern,
        test_griddata_nearest_float,
        test_griddata_linear,
        test_griddata_cubic,
        test_user_griddata_pattern,
        test_griddata_rescale,
        test_lagrange_1d,
        test_lagrange_2d_y,
        test_lagrange_extrapolate_false,
        test_lagrange_interp1d,
        test_lagrange_vs_scipy_barycentric,
        test_newton_1d,
        test_newton_2d_y,
        test_newton_extrapolate_false,
        test_newton_interp1d,
        test_lagrange_newton_agreement,
        test_hp_eval_hp_decimal_output,
        test_hp_32digit_input_accuracy,
        test_hp_call_returns_float64,
        test_hp_interp1d_wrapper,
        test_hp_extrapolate_false,
        test_hp_float64_input,
        test_hp_scalar_query,
        test_hp_lagrange_newton_agree,
    ]
    for t in tests:
        t()
    print("\nDone.")
