"""
Arc-length parameterization utilities for splines.
"""
import numpy as np
from scipy.integrate import quad, cumulative_trapezoid
from scipy.interpolate import splev, interp1d
from scipy.optimize import brentq


__all__ = ['splev_arclength', 'spline_arclength']


def _make_arc_derivative(tck):
    """Create arc-length derivative function for a spline.
    
    Parameters
    ----------
    tck : tuple
        A tuple (t, c, k) containing the vector of knots, the B-spline
        coefficients, and the degree of the spline.
    
    Returns
    -------
    arc_deriv : callable
        Function that computes ||dC/du|| at parameter u.
    parametric : bool
        Whether the curve is parametric (multi-dimensional).
    ndim : int
        Number of dimensions.
    u_min : float
        Minimum parameter value.
    u_max : float
        Maximum parameter value.
    """
    t, c, k = tck
    
    # Determine if parametric and extract dimensions
    try:
        c[0][0]
        parametric = True
        ndim = len(c)
    except (TypeError, IndexError):
        parametric = False
        ndim = 1
    
    # Get parameter range from knot vector
    u_min = t[k]
    u_max = t[-k-1]
    
    # Create appropriate derivative function
    if parametric:
        def arc_deriv(u):
            """Compute ||dC/du|| for parametric curve."""
            derivs = [splev(u, [t, c[i], k], der=1) for i in range(ndim)]
            return np.sqrt(sum(d**2 for d in derivs))
    else:
        def arc_deriv(u):
            """Compute |dy/du| for 1D curve."""
            return np.abs(splev(u, tck, der=1))
    
    return arc_deriv, parametric, ndim, u_min, u_max


def splev_arclength(tck, n_samples, u_range=None):
    """
    Evaluate a B-spline at points equally spaced along the arc length.

    Given a parametric spline curve, this function samples `n_samples` points
    that are equally spaced along the curve's arc length, rather than equally
    spaced in the parameter domain.

    Parameters
    ----------
    tck : tuple
        A tuple (t, c, k) containing the vector of knots, the B-spline
        coefficients, and the degree of the spline. For parametric curves,
        c should be a list of coefficient arrays for each dimension.
    n_samples : int
        Number of points to sample along the curve. Must be at least 1.
    u_range : tuple of float, optional
        (u_min, u_max) parameter range to consider. If None, uses the
        full knot range (t[k], t[-k-1]).

    Returns
    -------
    points : ndarray
        Array of shape (n_samples, ndim) for parametric curves, or
        (n_samples,) for 1D curves, containing the sampled points.
    u_values : ndarray
        Array of shape (n_samples,) containing the parameter values
        corresponding to each sampled point.

    Notes
    -----
    This function computes arc length by numerical integration of the
    curve's derivative. For parametric curves in n dimensions:

    .. math::
        s(u) = \\int_{u_0}^{u} \\sqrt{\\sum_{i=1}^{n} (dx_i/du)^2} du

    The algorithm:
    1. Pre-samples the curve densely to build a cumulative arc-length table
    2. Inverts this table via interpolation to map target arc lengths to
       parameter values
    3. Evaluates the spline at these parameter values

    This approach is O(n) in the number of samples, compared to O(n²) for
    naive root-finding methods.

    Examples
    --------
    Sample a 2D parametric curve (e.g., a circle) with equal arc-length spacing:

    >>> import numpy as np
    >>> from scipy.interpolate import splprep, splev_arclength
    >>> # Create a circle
    >>> theta = np.linspace(0, 2*np.pi, 20)
    >>> x = np.cos(theta)
    >>> y = np.sin(theta)
    >>> # Fit a periodic spline
    >>> tck, u = splprep([x, y], s=0, per=True)
    >>> # Sample 100 equally-spaced points along arc length
    >>> points, u_values = splev_arclength(tck, 100)
    >>> # points[i] and points[i+1] are approximately the same distance apart

    See Also
    --------
    splprep : Find the B-spline representation of an N-D curve.
    splev : Evaluate a B-spline or its derivatives.
    spline_arclength : Compute total arc length of a spline.

    References
    ----------
    .. [1] de Boor, C. (1978). A Practical Guide to Splines.
           Springer-Verlag.
    """
    # Input validation
    if not isinstance(n_samples, (int, np.integer)) or n_samples < 1:
        raise ValueError("n_samples must be a positive integer")
    
    # Extract arc-length derivative and metadata
    arc_deriv, parametric, ndim, u_min_default, u_max_default = _make_arc_derivative(tck)
    t, c, k = tck
    
    # Use provided range or default
    if u_range is None:
        u_min, u_max = u_min_default, u_max_default
    else:
        u_min, u_max = u_range
    
    # Special case: single sample
    if n_samples == 1:
        u_values = np.array([u_min])
        if parametric:
            points = np.array([splev(u_min, [t, c[i], k]) for i in range(ndim)]).T
        else:
            points = np.array([splev(u_min, tck)])
        return points, u_values
    
    # Pre-sample the curve to build cumulative arc-length table
    # Use adaptive density: more samples for longer curves
    n_dense = max(500, n_samples * 5)
    u_dense = np.linspace(u_min, u_max, n_dense)
    
    # Compute speeds at sample points
    speeds = np.array([arc_deriv(u) for u in u_dense])
    
    # Compute cumulative arc length using trapezoidal rule
    s_dense = cumulative_trapezoid(speeds, u_dense, initial=0)
    total_length = s_dense[-1]
    
    # Handle degenerate case
    if total_length == 0:
        u_values = np.linspace(u_min, u_max, n_samples)
    else:
        # Target arc lengths for each sample
        target_lengths = np.linspace(0, total_length, n_samples)
        
        # Invert cumulative arc length via interpolation
        # Use cubic interpolation for smooth parameter distribution
        s_to_u = interp1d(s_dense, u_dense, kind='cubic', 
                         bounds_error=False, fill_value=(u_min, u_max))
        u_values = s_to_u(target_lengths)
    
    # Evaluate the spline at the computed parameter values
    if parametric:
        points = np.column_stack([splev(u_values, [t, c[i], k]) 
                                 for i in range(ndim)])
    else:
        points = splev(u_values, tck)
    
    return points, u_values


def spline_arclength(tck, u_range=None):
    """
    Compute the arc length of a B-spline curve.

    Parameters
    ----------
    tck : tuple
        A tuple (t, c, k) containing the vector of knots, the B-spline
        coefficients, and the degree of the spline.
    u_range : tuple of float, optional
        (u_min, u_max) parameter range. If None, uses the full knot range.

    Returns
    -------
    length : float
        The arc length of the curve.

    Notes
    -----
    The arc length is computed by numerical integration of the curve's
    speed (magnitude of derivative) over the parameter domain.

    Examples
    --------
    >>> import numpy as np
    >>> from scipy.interpolate import splprep, spline_arclength
    >>> # Create a quarter circle
    >>> theta = np.linspace(0, np.pi/2, 10)
    >>> x = np.cos(theta)
    >>> y = np.sin(theta)
    >>> tck, u = splprep([x, y], s=0)
    >>> length = spline_arclength(tck)
    >>> np.abs(length - np.pi/2) < 0.01
    True

    See Also
    --------
    splev_arclength : Sample spline at equal arc-length intervals.
    splint : Integrate a B-spline between two points.
    """
    # Extract arc-length derivative and parameter range
    arc_deriv, _, _, u_min_default, u_max_default = _make_arc_derivative(tck)
    
    # Use provided range or default
    if u_range is None:
        u_min, u_max = u_min_default, u_max_default
    else:
        u_min, u_max = u_range
    
    # Integrate speed to get total length
    length, _ = quad(arc_deriv, u_min, u_max, limit=100)
    return length
