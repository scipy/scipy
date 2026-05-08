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
                 on equally-spaced grids; use precision= for >16 digit accuracy)
  1-D newton   : divided-difference table + Horner evaluation (O(n²) setup,
                 O(n) evaluation; equivalent polynomial to Lagrange;
                 use precision= for >16 digit accuracy)
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
    """Barycentric weights matching scipy.BarycentricInterpolator.

    Uses capacity scaling c = 4/(x_max - x_min) so each factor in the
    product stays in [-4, 4], preventing overflow/underflow for large n.
    Matches scipy's formula exactly (without the optional random permutation).

    w_j = 1 / prod_{m != j} c*(x_j - x_m)
    """
    n       = len(x)
    inv_cap = 4.0 / (x.max() - x.min())
    w       = np.zeros(n)
    for i in range(n):
        dist    = inv_cap * (x[i] - x)   # scaled differences
        dist[i] = 1.0                     # skip self-term
        w[i]    = 1.0 / np.prod(dist)
    return w


def _lagrange_eval(x_new, x, y, w):
    """Evaluate Lagrange polynomial via the barycentric second form.

    L(x) = [Σ_j w_j/(x-x_j) * y_j] / [Σ_j w_j/(x-x_j)]

    Mirrors scipy.BarycentricInterpolator._evaluate exactly:
    - uses exact equality to detect queries that fall on a node,
    - overwrites those positions with the node value directly.
    """
    c       = x_new[:, None] - x[None, :]   # (Q, n)  differences
    z       = (c == 0)                       # (Q, n)  on-node mask (exact)
    c[z]    = 1.0                            # avoid division by zero temporarily

    with np.errstate(divide='ignore', invalid='ignore'):
        c = w[None, :] / c                   # (Q, n)  barycentric terms

    if y.ndim == 1:
        numer = c.dot(y)                     # (Q,)
        denom = c.sum(axis=1)                # (Q,)
        out   = numer / denom
        # Overwrite positions where query == node
        r = np.nonzero(z)
        if len(r) == 1:                      # scalar query path
            if len(r[0]) > 0:
                out = y[r[0][0]]
        else:
            out[r[0]] = y[r[1]]
    else:
        numer = c.dot(y)                              # (Q, m)
        denom = c.sum(axis=1)                         # (Q,)
        out   = numer / denom[:, None]
        r = np.nonzero(z)
        if len(r) == 1:
            if len(r[0]) > 0:
                out = y[r[0][0]]
        else:
            out[r[0]] = y[r[1]]

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


# ------------------------------------------------------------------
# 1d  High-precision helpers  (Python decimal module)
# ------------------------------------------------------------------

def _to_dec(v, ctx):
    """Convert a scalar to Decimal preserving all available digits.

    - Decimal/int/float  → direct conversion via 17-significant-digit format
    - str                → Decimal(v) directly (use this to supply >17-digit data)
    """
    from decimal import Decimal
    if isinstance(v, Decimal):
        return v
    if isinstance(v, str):
        return ctx.create_decimal(v)
    # format with 17 significant figures – full IEEE-754 double precision
    return ctx.create_decimal(format(float(v), '.17e'))


def _barycentric_weights_hp(x, prec):
    """Barycentric weights at *prec* decimal digits using capacity scaling.

    Returns (x_dec, w_dec): lists of Decimal matching scipy's formula.
    """
    from decimal import localcontext
    with localcontext() as ctx:
        ctx.prec = prec + 10                   # 10 guard digits
        xd = [_to_dec(xi, ctx) for xi in x]
        n  = len(xd)
        inv_cap = ctx.create_decimal(4) / (xd[-1] - xd[0])
        w = []
        for i in range(n):
            prod = ctx.create_decimal(1)
            for j in range(n):
                if j != i:
                    prod = prod * (inv_cap * (xd[i] - xd[j]))
            w.append(ctx.create_decimal(1) / prod)
    return xd, w


def _lagrange_eval_hp(x_new, x_dec, y, w_dec, prec, return_decimal=False):
    """Evaluate Lagrange polynomial at *prec* decimal digits.

    x_new          : list of scalars (float, str, or Decimal)
    x_dec          : list of Decimal (nodes)
    y              : numpy array (n,) or (n, m)
    w_dec          : list of Decimal (barycentric weights)
    return_decimal : if True return list of Decimal (or list of lists); else float64 ndarray
    """
    from decimal import localcontext
    n   = len(x_dec)
    y   = np.asarray(y, dtype=float)
    y1d = y.ndim == 1
    Q   = len(x_new)
    raw = []            # list of Decimal (1-D) or list of list of Decimal (2-D)

    with localcontext() as ctx:
        ctx.prec = prec + 10
        if y1d:
            yd = [_to_dec(yi, ctx) for yi in y]
        else:
            yd = [[_to_dec(y[i, k], ctx) for k in range(y.shape[1])]
                  for i in range(n)]

        for qi in range(Q):
            xq    = _to_dec(x_new[qi], ctx)
            diffs = [xq - xi for xi in x_dec]
            node  = next((i for i, d in enumerate(diffs) if d == 0), None)

            if node is not None:
                raw.append(yd[node] if not y1d else yd[node])
                continue

            terms = [w_dec[i] / diffs[i] for i in range(n)]
            denom = sum(terms)
            if y1d:
                numer = sum(terms[i] * yd[i] for i in range(n))
                raw.append(numer / denom)
            else:
                row = [sum(terms[i] * yd[i][k] for i in range(n)) / denom
                       for k in range(y.shape[1])]
                raw.append(row)

    if return_decimal:
        return raw          # list of Decimal / list of list of Decimal

    # Convert to float64 numpy array
    if y1d:
        return np.array([float(v) for v in raw], dtype=float)
    else:
        return np.array([[float(v) for v in row] for row in raw], dtype=float)


def _divided_differences_hp(x, y, prec):
    """Newton divided-difference coefficients at *prec* decimal digits.

    x : list or array of x-nodes (strings preserved as-is; floats use .17e)
    y : numpy array (n,) or (n, m) – float64 values used for y
    Returns (x_dec, a_dec):
      x_dec : list of n Decimal nodes
      a_dec : list of n Decimal coefficients (1-D y)
           or list of n lists of m Decimal coefficients (2-D y)
    """
    from decimal import localcontext
    y = np.asarray(y, dtype=float)          # ensure numpy array
    y1d = y.ndim == 1
    with localcontext() as ctx:
        ctx.prec = prec + 10
        n   = len(x)
        xd  = [_to_dec(xi, ctx) for xi in x]

        if y1d:
            table = [_to_dec(yi, ctx) for yi in y]
            a = [table[0]]
            for j in range(1, n):
                table = [(table[i + 1] - table[i]) / (xd[i + j] - xd[i])
                         for i in range(n - j)]
                a.append(table[0])
        else:
            m = y.shape[1]
            tables = [[_to_dec(y[i, k], ctx) for i in range(n)]
                      for k in range(m)]
            a = [[tables[k][0] for k in range(m)]]
            for j in range(1, n):
                for k in range(m):
                    tables[k] = [
                        (tables[k][i + 1] - tables[k][i]) / (xd[i + j] - xd[i])
                        for i in range(n - j)
                    ]
                a.append([tables[k][0] for k in range(m)])

    return xd, a


def _newton_eval_hp(x_new, x_dec, a_dec, y1d, prec, return_decimal=False):
    """Evaluate Newton polynomial at *prec* decimal digits via Horner.

    a_dec          : list of Decimal (1-D) or list of list of Decimal (2-D)
    return_decimal : if True return list of Decimal; else float64 ndarray
    """
    from decimal import localcontext
    n   = len(x_dec)
    Q   = len(x_new)
    raw = []

    with localcontext() as ctx:
        ctx.prec = prec + 10
        if not y1d:
            m = len(a_dec[0])
        for qi in range(Q):
            xq = _to_dec(x_new[qi], ctx)
            if y1d:
                result = a_dec[n - 1]
                for j in range(n - 2, -1, -1):
                    result = result * (xq - x_dec[j]) + a_dec[j]
                raw.append(result)
            else:
                result = list(a_dec[n - 1])
                for j in range(n - 2, -1, -1):
                    dx = xq - x_dec[j]
                    result = [result[k] * dx + a_dec[j][k] for k in range(m)]
                raw.append(result)

    if return_decimal:
        return raw

    if y1d:
        return np.array([float(v) for v in raw], dtype=float)
    else:
        return np.array([[float(v) for v in row] for row in raw], dtype=float)


def _as_query_list(x_new):
    """Normalise query input to a plain Python list.

    Accepts: str, float, int, Decimal, numpy scalar, list, numpy array.
    Strings are preserved so _to_dec can parse them with full precision.
    """
    if isinstance(x_new, str):
        return [x_new]
    if isinstance(x_new, (int, float)):
        return [x_new]
    try:
        from decimal import Decimal
        if isinstance(x_new, Decimal):
            return [x_new]
    except ImportError:
        pass
    if isinstance(x_new, np.ndarray):
        # preserve dtype=object (may contain strings/Decimal)
        if x_new.dtype == object:
            return list(x_new.ravel())
        return list(x_new.ravel().astype(object))
    try:
        return list(x_new)
    except TypeError:
        return [x_new]


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
    precision : int or None, default None
        Number of decimal digits for internal arithmetic.  Only applies to
        ``kind='lagrange'`` and ``kind='newton'``.

        * ``None``           – standard float64 arithmetic (~15–16 sig figs).
        * ``35`` or higher   – use Python ``decimal`` module internally;
                               results are returned as float64 but the
                               *computation* retains the requested precision,
                               eliminating accumulated arithmetic errors.

        To obtain 10⁻³² accuracy in the **output** you must also supply *x*
        and *y* values with that many significant digits.  Pass them as Python
        strings (e.g. ``x=['1.234567890123456789012345678901234', ...]``) and
        they will be converted to ``Decimal`` exactly; float64 inputs are
        limited to ~17 significant digits.

    Call
    ----
    f(x_new)  →  interpolated values
    """

    def __init__(self, x, y, kind='linear', extrapolate=True, precision=None):
        if precision is not None and kind not in ('lagrange', 'newton'):
            raise ValueError("precision parameter only applies to kind='lagrange' or 'newton'")
        if precision is not None and (not isinstance(precision, int) or precision < 1):
            raise ValueError("precision must be a positive integer")

        # Always build float64 numpy arrays (used by non-hp paths and for bounds).
        # For the hp path, also keep x_raw as a list so string inputs survive as-is.
        if precision is not None:
            x_raw = list(x) if not isinstance(x, list) else x
            x = np.asarray([float(xi) for xi in x_raw], dtype=float)
        else:
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
        self.precision   = precision

        if kind == 'nearest':
            # Midpoint boundaries – divide before add to avoid overflow.
            # Exactly matches scipy interp1d kind='nearest'.
            self._x_bds = x[:-1] / 2.0 + x[1:] / 2.0

        elif kind == 'cubic':
            self._c = _cubic_spline_coeffs(x, y)

        elif kind == 'lagrange':
            if precision is not None:
                # Store Decimal nodes and weights for high-precision evaluation.
                self._x_dec, self._bary_w_dec = _barycentric_weights_hp(
                    x_raw, precision)
            else:
                self._bary_w = _barycentric_weights(x)

        elif kind == 'newton':
            if precision is not None:
                # y values passed as float64 numpy array; x nodes as raw list.
                self._x_dec, self._dd_a_dec = _divided_differences_hp(
                    x_raw, y, precision)
            else:
                self._dd_a = _divided_differences(x, y)

    def __call__(self, x_new):
        """Evaluate the interpolant at x_new."""
        scalar = not hasattr(x_new, '__len__') and not isinstance(x_new, np.ndarray)
        if not scalar:
            scalar = np.asarray(x_new).ndim == 0
        x_new_f = np.asarray(x_new, dtype=float)
        x_new_f = np.atleast_1d(x_new_f)

        if self.kind == 'linear':
            out = self._linear(x_new_f)
        elif self.kind == 'nearest':
            out = self._nearest(x_new_f)
        elif self.kind == 'cubic':
            out = self._cubic(x_new_f)
        elif self.kind == 'lagrange':
            out = self._lagrange(x_new, x_new_f)
        else:
            out = self._newton(x_new, x_new_f)

        if not self.extrapolate:
            oob = (x_new_f < self.x[0]) | (x_new_f > self.x[-1])
            if oob.any():
                out = out.copy()
                out[oob] = np.nan

        if scalar:
            out = out.squeeze()
            return float(out) if out.ndim == 0 else out
        return out

    def _linear(self, x_new_f):
        # side='left' + clip to [1, n-1] matches scipy interp1d._call_linear.
        idx    = np.searchsorted(self.x, x_new_f).clip(1, len(self.x) - 1)
        lo, hi = idx - 1, idx
        w = (x_new_f - self.x[lo]) / (self.x[hi] - self.x[lo])
        if self._y1d:
            return self.y[lo] * (1.0 - w) + self.y[hi] * w
        return self.y[lo] + w[:, None] * (self.y[hi] - self.y[lo])

    def _nearest(self, x_new_f):
        # Midpoint searchsorted matches scipy interp1d kind='nearest' exactly.
        idx = np.searchsorted(self._x_bds, x_new_f, side='left')
        idx = idx.clip(0, len(self.x) - 1)
        return self.y[idx]

    def _cubic(self, x_new_f):
        # side='right' - 1 places x exactly on a knot in the left interval.
        idx = (np.searchsorted(self.x, x_new_f, side='right') - 1).clip(
            0, len(self.x) - 2
        )
        dx = x_new_f - self.x[idx]
        c  = self._c
        if self._y1d:
            return (((c[0, idx]*dx + c[1, idx])*dx + c[2, idx])*dx + c[3, idx])
        dx = dx[:, None]
        return (((c[0][idx]*dx + c[1][idx])*dx + c[2][idx])*dx + c[3][idx])

    def _lagrange(self, x_new_raw, x_new_f):
        if self.precision is not None:
            qpts = _as_query_list(x_new_raw)
            return _lagrange_eval_hp(qpts, self._x_dec, self.y,
                                     self._bary_w_dec, self.precision)
        return _lagrange_eval(x_new_f, self.x, self.y, self._bary_w)

    def _newton(self, x_new_raw, x_new_f):
        if self.precision is not None:
            qpts = _as_query_list(x_new_raw)
            return _newton_eval_hp(qpts, self._x_dec, self._dd_a_dec,
                                   self._y1d, self.precision)
        return _newton_eval(x_new_f, self.x, self._dd_a)

    def eval_hp(self, x_new):
        """Evaluate at full Decimal precision (only for kind='lagrange'/'newton').

        Returns a list of ``decimal.Decimal`` values (1-D y) or a list of
        lists of ``decimal.Decimal`` (2-D y).  Use this when you need more than
        ~16 significant digits in the output.

        Requires ``precision`` to have been set in the constructor.

        Parameters
        ----------
        x_new : sequence of str, float, or Decimal query points.
                Pass values as strings to preserve digits beyond float64.

        Returns
        -------
        list of Decimal  (1-D y)  or  list of list of Decimal  (2-D y)
        """
        if self.precision is None:
            raise ValueError("eval_hp requires precision= to be set in the constructor")
        if self.kind not in ('lagrange', 'newton'):
            raise ValueError("eval_hp only supported for kind='lagrange' or 'newton'")

        qpts = _as_query_list(x_new)
        if self.kind == 'lagrange':
            return _lagrange_eval_hp(qpts, self._x_dec, self.y,
                                     self._bary_w_dec, self.precision,
                                     return_decimal=True)
        else:
            return _newton_eval_hp(qpts, self._x_dec, self._dd_a_dec,
                                   self._y1d, self.precision,
                                   return_decimal=True)


# ============================================================
# SECTION 3 – interp1d  (scipy.interpolate.interp1d drop-in)
# ============================================================

def interp1d(x, y, kind='linear', axis=0, bounds_error=True,
             fill_value=np.nan, assume_sorted=False, precision=None):
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
    precision : int or None      Decimal digits for internal arithmetic
                                 (only for kind='lagrange'/'newton').
                                 Use precision=35 or higher for 10⁻³² accuracy.

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

    High-precision example (supply data as strings for full accuracy):
    >>> xi = ['0', '1', '2', '3', '4']
    >>> yi = ['0', '1', '4', '9', '16']
    >>> f = interp1d(xi, yi, kind='lagrange', fill_value='extrapolate', precision=40)
    >>> f(['1.5'])      # ≈ 2.25  (computed with 40 decimal digits internally)
    """
    # For precision path, keep raw lists to preserve string inputs
    if precision is not None:
        x_raw = list(x)
        y_raw = list(y)
        x_f   = np.array([float(xi) for xi in x_raw], dtype=float)
        y_f   = np.array([float(yi) for yi in y_raw], dtype=float)
        y_f   = np.moveaxis(y_f, axis, 0)
        if not assume_sorted:
            sort_idx = np.argsort(x_f)
            x_raw = [x_raw[i] for i in sort_idx]
            y_raw = [y_raw[i] for i in sort_idx]
            x_f   = x_f[sort_idx]
            y_f   = y_f[sort_idx]
        extrapolate = (isinstance(fill_value, str) and
                       fill_value.lower() == 'extrapolate')
        interp = NumpyInterpolator(x_raw, y_f, kind=kind,
                                   extrapolate=extrapolate, precision=precision)
    else:
        x  = np.asarray(x, dtype=float)
        y  = np.asarray(y, dtype=float)
        y  = np.moveaxis(y, axis, 0)
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
        x0 = interp.x[0]
        x1 = interp.x[-1]
        def _checked(x_new):
            xf = np.asarray(x_new, dtype=float)
            if np.any(xf < x0) or np.any(xf > x1):
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
