"""
Numpy-only drop-in replacements for scipy.interpolate.interp1d and griddata.

Public API (matching scipy signatures):
  NumpyInterpolator(x, y, kind, extrapolate)   – 1-D interpolant object
  interp1d(x, y, kind, fill_value, bounds_error) – scipy.interpolate.interp1d clone
  griddata(points, values, xi, method, fill_value, rescale) – scipy.interpolate.griddata clone

Supported methods
  interp1d / NumpyInterpolator : 'linear', 'nearest', 'cubic', 'lagrange', 'newton'
  griddata                     : 'nearest', 'linear', 'cubic'

Algorithms
  1-D linear   : searchsorted + de Boor weighting  (exact match with scipy)
  1-D nearest  : midpoint-boundary searchsorted    (exact match with scipy)
  1-D cubic    : not-a-knot C² spline via Thomas-algorithm tridiagonal solve
                 (exact match with scipy.interpolate.CubicSpline default)
  1-D lagrange : barycentric second form (numerically stable; O(n²) setup,
                 O(n) evaluation; exact at nodes; subject to Runge's phenomenon
                 on equally-spaced grids)
  1-D newton   : divided-difference table + Horner evaluation (O(n²) setup,
                 O(n) evaluation; equivalent polynomial to Lagrange)
  2-D nearest  : vectorised Euclidean distance (exact match)
  2-D linear   : Bowyer-Watson Delaunay + barycentric interpolation
                 (matches scipy LinearNDInterpolator; fill_value outside hull)
  2-D cubic    : Thin-plate spline RBF (C∞, exact at data points;
                 smooth approximation to scipy CloughTocher2DInterpolator)
"""

import numpy as np


# ============================================================
# SECTION 1 – 1-D helpers
# ============================================================

def _thomas_solve(lower, main, upper, rhs):
    """Tridiagonal system solver (Thomas / LU algorithm).

    lower : (n-1,)  sub-diagonal
    main  : (n,)    main diagonal
    upper : (n-1,)  super-diagonal
    rhs   : (n, m)  right-hand side
    Returns solution of shape (n, m).
    """
    n = len(main)
    b = main.astype(float, copy=True)
    d = rhs.astype(float, copy=True)

    for i in range(1, n):
        w = lower[i - 1] / b[i - 1]
        b[i] -= w * upper[i - 1]
        d[i] -= w * d[i - 1]

    x = np.empty_like(d)
    x[n - 1] = d[n - 1] / b[n - 1]
    for i in range(n - 2, -1, -1):
        x[i] = (d[i] - upper[i] * x[i + 1]) / b[i]
    return x


def _cubic_spline_coeffs(x, y):
    """Not-a-knot cubic spline coefficients matching scipy.CubicSpline.

    On interval i the polynomial is:
        S(t) = c[0,i]*(t-x[i])^3 + c[1,i]*(t-x[i])^2
             + c[2,i]*(t-x[i])   + c[3,i]

    Returns c of shape (4, n-1) for 1-D y, or (4, n-1, m) for 2-D y.
    """
    n = len(x)
    y1d = y.ndim == 1
    Y = y[:, np.newaxis] if y1d else y          # (n, m)
    m = Y.shape[1]

    h     = np.diff(x)                              # (n-1,)
    slope = np.diff(Y, axis=0) / h[:, np.newaxis]  # (n-1, m)

    if n == 2:
        # Single interval → linear polynomial (both derivatives = slope)
        s = np.vstack([slope, slope])               # (2, m)

    elif n == 3:
        # not-a-knot with 3 points: solve explicit 3×3 system
        A = np.array([
            [1.0,  1.0,           0.0],
            [h[1], 2*(h[0]+h[1]), h[0]],
            [0.0,  1.0,           1.0],
        ])
        b3 = np.zeros((3, m))
        b3[0] = 2 * slope[0]
        b3[1] = 3 * (h[1] * slope[0] + h[0] * slope[1])
        b3[2] = 2 * slope[1]
        s = np.linalg.solve(A, b3)                  # (3, m)

    else:
        # General case: banded tridiagonal + not-a-knot boundary conditions.
        # Mirrors scipy's CubicSpline solve_banded path exactly.
        d_s = x[2]  - x[0]    # h[0] + h[1]
        d_e = x[-1] - x[-3]   # h[-2] + h[-1]

        main_d = np.zeros(n)
        upper  = np.zeros(n - 1)
        lower  = np.zeros(n - 1)
        rhs    = np.zeros((n, m))

        # Interior rows 1 .. n-2
        main_d[1:-1] = 2 * (h[:-1] + h[1:])
        upper[1:]    = h[:-1]
        lower[:-1]   = h[1:]
        rhs[1:-1]    = 3 * (h[1:, None] * slope[:-1] + h[:-1, None] * slope[1:])

        # not-a-knot start (row 0)
        main_d[0] = h[1]
        upper[0]  = d_s
        rhs[0]    = ((h[0] + 2*d_s)*h[1]*slope[0] + h[0]**2*slope[1]) / d_s

        # not-a-knot end (row n-1)
        main_d[-1] = h[-2]
        lower[-1]  = d_e
        rhs[-1]    = (h[-1]**2*slope[-2]
                      + (2*d_e + h[-1])*h[-2]*slope[-1]) / d_e

        s = _thomas_solve(lower, main_d, upper, rhs)   # (n, m)

    # Convert Hermite data (Y, s) → PPoly coefficients.
    # Mirrors scipy CubicHermiteSpline.__init__ exactly.
    hr = h[:, np.newaxis]                           # (n-1, 1)
    t  = (s[:-1] + s[1:] - 2 * slope) / hr         # (n-1, m)

    c = np.stack([                                  # (4, n-1, m)
        t / hr,
        (slope - s[:-1]) / hr - t,
        s[:-1],
        Y[:-1],
    ])
    return c[:, :, 0] if y1d else c


# ------------------------------------------------------------------
# 1b  Lagrange polynomial helpers  (barycentric second form)
# ------------------------------------------------------------------

def _barycentric_weights(x):
    """Barycentric weights for Lagrange interpolation.

    w_j = 1 / prod_{m != j} (x_j - x_m), scaled by max(|w|) for stability.
    """
    n    = len(x)
    diff = x[:, None] - x[None, :]    # (n, n)  diff[j, m] = x_j - x_m
    np.fill_diagonal(diff, 1.0)        # avoid 0 on diagonal
    w    = 1.0 / diff.prod(axis=1)     # (n,)
    w   /= np.abs(w).max()
    return w


def _lagrange_eval(x_new, x, y, w):
    """Evaluate Lagrange polynomial via the barycentric second form.

    L(x) = [Σ_j w_j/(x-x_j) * y_j] / [Σ_j w_j/(x-x_j)]

    Exactly returns y_j when x equals a node x_j.
    """
    diff    = x_new[:, None] - x[None, :]          # (Q, n)
    atol    = np.finfo(float).eps * (np.abs(x).max() + 1.0) * 16.0
    on_node = np.abs(diff) < atol                   # (Q, n)

    with np.errstate(divide='ignore', invalid='ignore'):
        terms = w[None, :] / diff                   # (Q, n)
    terms[on_node] = 0.0                            # prevent 0*inf → nan

    if y.ndim == 1:
        numer = (terms * y[None, :]).sum(axis=1)    # (Q,)
        denom = terms.sum(axis=1)                   # (Q,)
        denom_safe = np.where(np.abs(denom) < 1e-300, 1.0, denom)
        out = numer / denom_safe
        for qi in np.where(on_node.any(axis=1))[0]:
            out[qi] = y[np.argmax(on_node[qi])]
    else:
        numer = (terms[:, :, None] * y[None, :, :]).sum(axis=1)  # (Q, m)
        denom = terms.sum(axis=1)                                  # (Q,)
        denom_safe = np.where(np.abs(denom) < 1e-300, 1.0, denom)
        out = numer / denom_safe[:, None]
        for qi in np.where(on_node.any(axis=1))[0]:
            out[qi] = y[np.argmax(on_node[qi])]

    return out


# ------------------------------------------------------------------
# 1c  Newton polynomial helpers  (divided differences + Horner)
# ------------------------------------------------------------------

def _divided_differences(x, y):
    """Compute Newton divided-difference coefficients.

    Returns a where a[j] is the j-th divided difference [y_0, ..., y_j].
    Shape: (n,) for 1-D y, or (n, m) for 2-D y.
    """
    n   = len(x)
    y1d = y.ndim == 1
    table = np.asarray(y, dtype=float).copy()
    if y1d:
        table = table[:, np.newaxis]               # (n, 1) for uniform code

    a = np.zeros_like(table)
    a[0] = table[0]

    for j in range(1, n):
        table[:n - j] = ((table[1:n - j + 1] - table[:n - j])
                         / (x[j:] - x[:n - j])[:, np.newaxis])
        a[j] = table[0]

    return a[:, 0] if y1d else a


def _newton_eval(x_new, x, a):
    """Evaluate Newton polynomial at x_new via Horner's method.

    p(x) = a[0] + (x-x[0])*(a[1] + (x-x[1])*(... + (x-x[n-2])*a[n-1]...))
    """
    n   = len(x)
    y1d = a.ndim == 1
    A   = a[:, np.newaxis] if y1d else a            # (n, m)

    result = np.tile(A[n - 1], (len(x_new), 1)).astype(float)  # (Q, m)
    for j in range(n - 2, -1, -1):
        result = result * (x_new - x[j])[:, np.newaxis] + A[j][np.newaxis, :]

    return result[:, 0] if y1d else result


# ============================================================
# SECTION 2 – NumpyInterpolator  (1-D)
# ============================================================

class NumpyInterpolator:
    """1-D interpolation using NumPy only, matching SciPy behavior.

    Parameters
    ----------
    x : array_like (n,)
        Strictly increasing sample points.
    y : array_like (n,) or (n, m)
        Sample values.
    kind : {'linear', 'nearest', 'cubic', 'lagrange', 'newton'}
    extrapolate : bool, default True
        Extrapolate outside [x[0], x[-1]] using the first/last interval.
        If False, out-of-bounds queries return NaN.

    Call
    ----
    f(x_new)  →  interpolated values
    """

    def __init__(self, x, y, kind='linear', extrapolate=True):
        x = np.asarray(x, dtype=float)
        y = np.asarray(y, dtype=float)

        if x.ndim != 1:
            raise ValueError("x must be 1-D")
        if y.shape[0] != len(x):
            raise ValueError("y.shape[0] must equal len(x)")
        if len(x) < 2:
            raise ValueError("Need at least 2 data points")
        if kind not in ('linear', 'nearest', 'cubic', 'lagrange', 'newton'):
            raise ValueError(
                "kind must be 'linear', 'nearest', 'cubic', 'lagrange', or 'newton'"
            )

        self.x           = x
        self.y           = y
        self.kind        = kind
        self.extrapolate = extrapolate
        self._y1d        = (y.ndim == 1)

        if kind == 'nearest':
            # Midpoint boundaries – divide before add to avoid overflow.
            # Exactly matches scipy interp1d kind='nearest'.
            self._x_bds = x[:-1] / 2.0 + x[1:] / 2.0

        elif kind == 'cubic':
            self._c = _cubic_spline_coeffs(x, y)

        elif kind == 'lagrange':
            # Precompute barycentric weights O(n²); evaluation is O(n).
            self._bary_w = _barycentric_weights(x)

        elif kind == 'newton':
            # Precompute divided-difference table O(n²); Horner evaluation O(n).
            self._dd_a = _divided_differences(x, y)

    def __call__(self, x_new):
        """Evaluate the interpolant at x_new."""
        x_new  = np.asarray(x_new, dtype=float)
        scalar = x_new.ndim == 0
        x_new  = np.atleast_1d(x_new)

        if self.kind == 'linear':
            out = self._linear(x_new)
        elif self.kind == 'nearest':
            out = self._nearest(x_new)
        elif self.kind == 'cubic':
            out = self._cubic(x_new)
        elif self.kind == 'lagrange':
            out = self._lagrange(x_new)
        else:
            out = self._newton(x_new)

        if not self.extrapolate:
            oob = (x_new < self.x[0]) | (x_new > self.x[-1])
            if oob.any():
                out = out.copy()
                out[oob] = np.nan

        if scalar:
            out = out.squeeze()
            return float(out) if out.ndim == 0 else out
        return out

    def _linear(self, x_new):
        # side='left' + clip to [1, n-1] matches scipy interp1d._call_linear.
        idx    = np.searchsorted(self.x, x_new).clip(1, len(self.x) - 1)
        lo, hi = idx - 1, idx
        w = (x_new - self.x[lo]) / (self.x[hi] - self.x[lo])
        if self._y1d:
            return self.y[lo] * (1.0 - w) + self.y[hi] * w
        return self.y[lo] + w[:, None] * (self.y[hi] - self.y[lo])

    def _nearest(self, x_new):
        # Midpoint searchsorted matches scipy interp1d kind='nearest' exactly.
        idx = np.searchsorted(self._x_bds, x_new, side='left')
        idx = idx.clip(0, len(self.x) - 1)
        return self.y[idx]

    def _cubic(self, x_new):
        # side='right' - 1 places x exactly on a knot in the left interval.
        idx = (np.searchsorted(self.x, x_new, side='right') - 1).clip(
            0, len(self.x) - 2
        )
        dx = x_new - self.x[idx]
        c  = self._c
        if self._y1d:
            return (((c[0, idx]*dx + c[1, idx])*dx + c[2, idx])*dx + c[3, idx])
        dx = dx[:, None]
        return (((c[0][idx]*dx + c[1][idx])*dx + c[2][idx])*dx + c[3][idx])

    def _lagrange(self, x_new):
        return _lagrange_eval(x_new, self.x, self.y, self._bary_w)

    def _newton(self, x_new):
        return _newton_eval(x_new, self.x, self._dd_a)


# ============================================================
# SECTION 3 – interp1d  (scipy.interpolate.interp1d drop-in)
# ============================================================

def interp1d(x, y, kind='linear', axis=0, bounds_error=True,
             fill_value=np.nan, assume_sorted=False):
    """1-D interpolation function matching scipy.interpolate.interp1d.

    Parameters
    ----------
    x : (n,) array_like          Strictly monotonic sample x-coordinates.
    y : array_like               Sample values; interpolation axis is `axis`.
    kind : str                   'linear', 'nearest', 'cubic', 'lagrange', or 'newton'.
    axis : int                   Axis of y along which to interpolate (default 0).
    bounds_error : bool          Raise ValueError for out-of-bounds if True.
    fill_value : float or 'extrapolate'
                                 Fill for out-of-bounds when bounds_error=False.
                                 Use 'extrapolate' to extrapolate.
    assume_sorted : bool         Skip sorting if True.

    Returns
    -------
    Callable f such that f(x_new) returns interpolated values.

    Examples
    --------
    >>> import numpy as np
    >>> xi = np.array([0., 1., 2., 3.])
    >>> yi = np.array([0., 1., 4., 9.])
    >>> f = interp1d(xi, yi, kind='cubic', fill_value='extrapolate')
    >>> f(1.5)          # ≈ 2.25
    """
    x  = np.asarray(x, dtype=float)
    y  = np.asarray(y, dtype=float)

    # Move interpolation axis to position 0
    y = np.moveaxis(y, axis, 0)

    if not assume_sorted:
        sort_idx = np.argsort(x)
        x = x[sort_idx]
        y = y[sort_idx]

    extrapolate = (isinstance(fill_value, str) and
                   fill_value.lower() == 'extrapolate')
    interp = NumpyInterpolator(x, y, kind=kind, extrapolate=extrapolate)

    if bounds_error and not extrapolate:
        # Python resolves __call__ on the class, not the instance,
        # so return a closure instead of patching the instance attribute.
        def _checked(x_new):
            x_new = np.asarray(x_new, dtype=float)
            if np.any(x_new < x[0]) or np.any(x_new > x[-1]):
                raise ValueError(
                    "A value in x_new is below/above the interpolation range."
                )
            return interp(x_new)
        return _checked

    return interp


# ============================================================
# SECTION 4 – 2-D scattered interpolation helpers
# ============================================================

# ------------------------------------------------------------------
# 4a  Nearest-neighbor
# ------------------------------------------------------------------

def _nearest_2d(pts, vals, xi):
    """Return value of nearest data point for each query in xi.

    pts  : (N, 2)  data coordinates
    vals : (N,) or (N, m)
    xi   : (Q, 2)  query coordinates
    """
    # Chunked to limit peak memory (each chunk ≈ CHUNK*N*2*8 bytes)
    CHUNK = 2048
    Q = len(xi)
    out = np.empty((Q,) + vals.shape[1:], dtype=vals.dtype)

    for start in range(0, Q, CHUNK):
        end  = min(start + CHUNK, Q)
        diff = xi[start:end, None, :] - pts[None, :, :]      # (Qc, N, 2)
        d2   = (diff * diff).sum(axis=2)                      # (Qc, N)
        idx  = d2.argmin(axis=1)                              # (Qc,)
        out[start:end] = vals[idx]

    return out


# ------------------------------------------------------------------
# 4b  Delaunay triangulation (Bowyer-Watson)
# ------------------------------------------------------------------

def _delaunay_2d(pts):
    """2-D Delaunay triangulation using the Bowyer-Watson algorithm.

    pts : (N, 2) array of data point coordinates.
    Returns triangles : (T, 3) int array of vertex indices into pts.
    """
    pts = np.asarray(pts, dtype=float)
    n   = len(pts)

    if n < 3:
        raise ValueError("Delaunay triangulation requires at least 3 points")

    # Normalise to [0,1]² for numerical stability
    lo    = pts.min(0)
    hi    = pts.max(0)
    span  = hi - lo
    span[span < 1e-12] = 1.0
    p = (pts - lo) / span                       # normalised coords

    # Super-triangle vertices at indices n, n+1, n+2
    big = 100.0
    sp  = np.array([[-big, -big], [2*big, -big], [-big, 2*big]], dtype=float)
    all_p = np.vstack([p, sp])

    def _ccw(t):
        """Return [i, j, k] in counter-clockwise order."""
        a, b, c = all_p[t[0]], all_p[t[1]], all_p[t[2]]
        if (b[0]-a[0])*(c[1]-a[1]) - (b[1]-a[1])*(c[0]-a[0]) < 0:
            return [t[0], t[2], t[1]]
        return list(t)

    def _in_circumcircle(t, d):
        """True iff point d is strictly inside the circumcircle of CCW triangle t."""
        # Lifted-paraboloid determinant test (robust, O(1)):
        #   | ax-dx  ay-dy  (ax-dx)²+(ay-dy)² |
        #   | bx-dx  by-dy  (bx-dx)²+(by-dy)² | > 0
        #   | cx-dx  cy-dy  (cx-dx)²+(cy-dy)² |
        a = all_p[t[0]] - d
        b = all_p[t[1]] - d
        c = all_p[t[2]] - d
        ar2 = a[0]*a[0] + a[1]*a[1]
        br2 = b[0]*b[0] + b[1]*b[1]
        cr2 = c[0]*c[0] + c[1]*c[1]
        det = (a[0]*(b[1]*cr2 - c[1]*br2)
               - a[1]*(b[0]*cr2 - c[0]*br2)
               + ar2 *(b[0]*c[1] - c[0]*b[1]))
        return det > 0.0

    # Start with super-triangle
    tris = [[n, n + 1, n + 2]]

    for i in range(n):
        pt  = p[i]

        # Find triangles whose circumcircle contains pt
        bad = [j for j, t in enumerate(tris) if _in_circumcircle(t, pt)]

        # Collect boundary edges (shared by exactly one bad triangle)
        edge_cnt = {}
        for j in bad:
            t = tris[j]
            for e in ((t[0], t[1]), (t[1], t[2]), (t[2], t[0])):
                key = (min(e), max(e))
                edge_cnt[key] = edge_cnt.get(key, 0) + 1
        boundary = [e for e, cnt in edge_cnt.items() if cnt == 1]

        # Remove bad triangles (high index first to preserve positions)
        for j in sorted(bad, reverse=True):
            tris.pop(j)

        # Re-triangulate the hole
        for e in boundary:
            tris.append(_ccw([e[0], e[1], i]))

    # Discard triangles that share a vertex with the super-triangle
    result = [t for t in tris if all(v < n for v in t)]
    if len(result) == 0:
        raise ValueError("Delaunay triangulation produced no triangles "
                         "(all points may be collinear)")
    return np.array(result, dtype=int)


# ------------------------------------------------------------------
# 4c  Barycentric linear interpolation on a triangulation
# ------------------------------------------------------------------

def _linear_2d(pts, vals, xi, fill_value):
    """Linear barycentric interpolation matching scipy.LinearNDInterpolator.

    pts        : (N, 2)  data point coordinates
    vals       : (N,) or (N, m)  data values
    xi         : (Q, 2)  query coordinates
    fill_value : scalar  used for points outside the convex hull
    """
    # Deduplicate: keep first occurrence of each unique (x, y) pair.
    # Duplicate points with different values are inconsistent for linear interp;
    # keeping the first occurrence matches scipy's behaviour.
    _, uniq_idx = np.unique(pts, axis=0, return_index=True)
    if len(uniq_idx) < len(pts):
        uniq_idx = np.sort(uniq_idx)   # preserve original ordering
        pts  = pts[uniq_idx]
        vals = vals[uniq_idx]

    tris = _delaunay_2d(pts)            # (T, 3)
    T    = len(tris)
    Q    = len(xi)

    P0 = pts[tris[:, 0]]               # (T, 2)
    P1 = pts[tris[:, 1]]
    P2 = pts[tris[:, 2]]

    # Precompute denominator for barycentric coordinates (same for all queries)
    denom = ((P1[:, 1] - P2[:, 1]) * (P0[:, 0] - P2[:, 0]) +
             (P2[:, 0] - P1[:, 0]) * (P0[:, 1] - P2[:, 1]))   # (T,)

    val_dtype = np.result_type(vals, 1.0)
    out = np.full((Q,) + vals.shape[1:], fill_value, dtype=val_dtype)

    # Process queries in chunks to bound memory use
    CHUNK = 2048
    for start in range(0, Q, CHUNK):
        end = min(start + CHUNK, Q)
        Qc  = end - start

        qx = xi[start:end, 0]          # (Qc,)
        qy = xi[start:end, 1]

        dx = qx[:, None] - P2[None, :, 0]   # (Qc, T)
        dy = qy[:, None] - P2[None, :, 1]

        n01 = (P1[:, 1] - P2[:, 1])[None, :]   # (1, T)
        n11 = (P2[:, 0] - P1[:, 0])[None, :]
        n02 = (P2[:, 1] - P0[:, 1])[None, :]
        n12 = (P0[:, 0] - P2[:, 0])[None, :]

        l0 = (n01*dx + n11*dy) / denom[None, :]   # (Qc, T)
        l1 = (n02*dx + n12*dy) / denom[None, :]
        l2 = 1.0 - l0 - l1

        eps    = -1e-10
        inside = (l0 >= eps) & (l1 >= eps) & (l2 >= eps)   # (Qc, T)

        valid   = inside.any(axis=1)                        # (Qc,)
        tri_idx = np.argmax(inside, axis=1)                 # (Qc,) first match

        qv = np.where(valid)[0]
        tv = tri_idx[qv]

        lv0 = l0[qv, tv]
        lv1 = l1[qv, tv]
        lv2 = l2[qv, tv]

        v0 = vals[tris[tv, 0]]
        v1 = vals[tris[tv, 1]]
        v2 = vals[tris[tv, 2]]

        if vals.ndim == 1:
            out[start + qv] = lv0*v0 + lv1*v1 + lv2*v2
        else:
            out[start + qv] = (lv0[:, None]*v0 + lv1[:, None]*v1
                                + lv2[:, None]*v2)

    return out


# ------------------------------------------------------------------
# 4d  Thin-plate spline RBF (cubic-equivalent for 2-D scattered data)
# ------------------------------------------------------------------

def _tps_phi(r2):
    """Thin-plate spline kernel: φ(r) = r² ln(r), with φ(0)=0."""
    with np.errstate(divide='ignore', invalid='ignore'):
        out = 0.5 * r2 * np.log(r2)
    out[r2 < 1e-300] = 0.0
    return out


def _cubic_2d(pts, vals, xi):
    """Thin-plate spline RBF interpolation (2-D cubic).

    Gives exact interpolation at data points and a C∞ smooth surface.
    Approximates scipy CloughTocher2DInterpolator for interior queries;
    extrapolates globally (no fill_value region).

    pts  : (N, 2)
    vals : (N,) or (N, m)
    xi   : (Q, 2)
    """
    N = len(pts)

    # --- Build (N+3) × (N+3) TPS system ---
    diff   = pts[:, None, :] - pts[None, :, :]           # (N, N, 2)
    r2_mat = (diff * diff).sum(axis=2)                    # (N, N)
    Phi    = _tps_phi(r2_mat.copy())                      # (N, N)

    P = np.column_stack([np.ones(N), pts[:, 0], pts[:, 1]])  # (N, 3)

    M          = np.zeros((N + 3, N + 3))
    M[:N, :N]  = Phi
    M[:N, N:]  = P
    M[N:, :N]  = P.T
    # M[N:, N:] = 0   (already zero)

    rhs      = np.zeros((N + 3,) + vals.shape[1:], dtype=float)
    rhs[:N]  = vals.real if np.iscomplexobj(vals) else vals

    try:
        coeffs_r = np.linalg.solve(M, rhs.reshape(N + 3, -1))
    except np.linalg.LinAlgError:
        coeffs_r, *_ = np.linalg.lstsq(M, rhs.reshape(N + 3, -1), rcond=None)

    if np.iscomplexobj(vals):
        rhs[:N] = vals.imag
        try:
            coeffs_i = np.linalg.solve(M, rhs.reshape(N + 3, -1))
        except np.linalg.LinAlgError:
            coeffs_i, *_ = np.linalg.lstsq(M, rhs.reshape(N + 3, -1), rcond=None)

    w_r = coeffs_r[:N]     # (N, ...)  RBF weights
    a_r = coeffs_r[N:]     # (3, ...)  polynomial coefficients

    # --- Evaluate at query points in chunks ---
    Q     = len(xi)
    CHUNK = 1024
    out_r = np.empty((Q,) + vals.shape[1:], dtype=float)

    for start in range(0, Q, CHUNK):
        end  = min(start + CHUNK, Q)
        diff_q = xi[start:end, None, :] - pts[None, :, :]  # (Qc, N, 2)
        r2_q   = (diff_q * diff_q).sum(axis=2)              # (Qc, N)
        Phi_q  = _tps_phi(r2_q.copy())                      # (Qc, N)

        P_q = np.column_stack([
            np.ones(end - start),
            xi[start:end, 0],
            xi[start:end, 1],
        ])                                                    # (Qc, 3)

        out_r[start:end] = (Phi_q @ w_r + P_q @ a_r).reshape(
            (end - start,) + vals.shape[1:]
        )

    if not np.iscomplexobj(vals):
        return out_r

    # Imaginary part
    out_i = np.empty_like(out_r)
    w_i   = coeffs_i[:N]
    a_i   = coeffs_i[N:]
    for start in range(0, Q, CHUNK):
        end    = min(start + CHUNK, Q)
        diff_q = xi[start:end, None, :] - pts[None, :, :]
        r2_q   = (diff_q * diff_q).sum(axis=2)
        Phi_q  = _tps_phi(r2_q.copy())
        P_q = np.column_stack([
            np.ones(end - start),
            xi[start:end, 0],
            xi[start:end, 1],
        ])
        out_i[start:end] = (Phi_q @ w_i + P_q @ a_i).reshape(
            (end - start,) + vals.shape[1:]
        )

    return out_r + 1j * out_i


# ============================================================
# SECTION 5 – griddata  (scipy.interpolate.griddata drop-in)
# ============================================================

def griddata(points, values, xi, method='linear', fill_value=np.nan,
             rescale=False):
    """Interpolate unstructured 2-D data.

    Matches scipy.interpolate.griddata for 2-D scattered data.

    Parameters
    ----------
    points : (N, 2) array_like or length-2 tuple of (N,) arrays
        Known data point coordinates.
    values : (N,) or (N, m) array_like
        Data values at the known points.  Complex arrays are supported.
    xi : (M, 2) array_like  **or**  length-2 tuple of arrays (e.g. meshgrid)
        Query coordinates.  A tuple ``(XX, YY)`` from ``np.meshgrid`` is
        flattened internally and the result is reshaped to match ``XX.shape``.
    method : {'linear', 'nearest', 'cubic'}
        * 'nearest' – value at closest data point (no fill_value region).
        * 'linear'  – Delaunay triangulation + barycentric interpolation;
                      fill_value used outside the convex hull.
        * 'cubic'   – thin-plate spline RBF (C∞, exact at data points;
                      smooth approximation to scipy CloughTocher).
    fill_value : float, default np.nan
        Value for queries outside the convex hull ('linear' only).
    rescale : bool, default False
        Rescale coordinates to unit cube for numerical stability.

    Returns
    -------
    ndarray
        Interpolated values.  Shape matches xi (with trailing value dims).

    Examples
    --------
    >>> import numpy as np
    >>> rng = np.random.default_rng(0)
    >>> pts = rng.random((50, 2))
    >>> vals = np.sin(pts[:, 0]) * np.cos(pts[:, 1])
    >>> XX, YY = np.meshgrid(np.linspace(0,1,20), np.linspace(0,1,20))
    >>> Z = griddata(pts, vals, (XX, YY), method='linear')
    >>> Z.shape
    (20, 20)
    """
    # ---------- normalise points ----------
    if isinstance(points, tuple):
        points = np.column_stack([np.asarray(a, dtype=float).ravel()
                                  for a in points])
    else:
        points = np.asarray(points, dtype=float)
        if points.ndim == 1:
            points = points[:, np.newaxis]

    values = np.asarray(values)
    if values.ndim == 0 or len(values) != len(points):
        raise ValueError("values must have the same length as points")

    # ---------- normalise xi ----------
    if isinstance(xi, tuple):
        xi_arrs   = [np.asarray(a, dtype=float) for a in xi]
        out_shape = xi_arrs[0].shape
        xi_flat   = np.column_stack([a.ravel() for a in xi_arrs])
    else:
        xi_flat   = np.asarray(xi, dtype=float)
        if xi_flat.ndim == 1:
            xi_flat = xi_flat[:, np.newaxis]
        out_shape = xi_flat.shape[:-1]
        xi_flat   = xi_flat.reshape(-1, xi_flat.shape[-1])

    if points.shape[1] != xi_flat.shape[1]:
        raise ValueError("points and xi must have the same number of dimensions")

    # ---------- optional rescale ----------
    if rescale:
        pmin = points.min(0)
        pmax = points.max(0)
        span = pmax - pmin
        span[span < 1e-12] = 1.0
        points  = (points  - pmin) / span
        xi_flat = (xi_flat - pmin) / span

    # ---------- dispatch ----------
    method = method.lower()

    if method == 'nearest':
        out_flat = _nearest_2d(points, values, xi_flat)

    elif method == 'linear':
        out_flat = _linear_2d(points, values, xi_flat, fill_value)

    elif method == 'cubic':
        ndim = points.shape[1]
        if ndim != 2:
            raise ValueError("cubic griddata is only supported for 2-D data")
        out_flat = _cubic_2d(points, values, xi_flat)

    else:
        raise ValueError(
            f"method must be 'nearest', 'linear', or 'cubic', got {method!r}"
        )

    # ---------- reshape output ----------
    return out_flat.reshape(out_shape + values.shape[1:])
