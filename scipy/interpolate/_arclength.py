"""
Arc-length parameterization utilities for splines.
"""
import numpy as np
from scipy.integrate import quad
from scipy.interpolate import splev
from scipy.optimize import brentq


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
        Number of points to sample along the curve.
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
    1. Computes the total arc length of the curve
    2. Divides it into n_samples equal segments
    3. For each segment, finds the parameter value u that corresponds
       to that arc length using root-finding

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

    References
    ----------
    .. [1] de Boor, C. (1978). A Practical Guide to Splines.
           Springer-Verlag.
    """
    t, c, k = tck
    
    # Determine if this is a parametric curve
    try:
        c[0][0]
        parametric = True
        ndim = len(c)
    except (TypeError, IndexError):
        parametric = False
        ndim = 1
    
    # Determine parameter range
    if u_range is None:
        u_min = t[k]
        u_max = t[-k-1]
    else:
        u_min, u_max = u_range
    
    # Function to compute the arc length derivative (speed) at parameter u
    def arc_length_derivative(u):
        """Compute ||dC/du|| at parameter u."""
        if parametric:
            # Evaluate first derivative for each dimension
            derivs = [splev(u, [t, c[i], k], der=1) for i in range(ndim)]
            # Compute Euclidean norm
            return np.sqrt(sum(d**2 for d in derivs))
        else:
            # For 1D curves, just absolute value of derivative
            deriv = splev(u, tck, der=1)
            return np.abs(deriv)
    
    # Compute cumulative arc length from u_min to a given u
    def cumulative_arc_length(u):
        """Integrate arc length from u_min to u."""
        if u <= u_min:
            return 0.0
        result, _ = quad(arc_length_derivative, u_min, u, limit=100)
        return result
    
    # Compute total arc length
    total_length = cumulative_arc_length(u_max)
    
    if total_length == 0:
        # Degenerate case: curve has no length
        u_values = np.linspace(u_min, u_max, n_samples)
    else:
        # Target arc lengths for each sample
        target_lengths = np.linspace(0, total_length, n_samples)
        
        # Find parameter values corresponding to each target arc length
        u_values = np.zeros(n_samples)
        u_values[0] = u_min
        u_values[-1] = u_max
        
        for i in range(1, n_samples - 1):
            target = target_lengths[i]
            
            # Find u such that cumulative_arc_length(u) = target
            # Use previous u as starting point for search
            u_search_min = u_values[i-1]
            u_search_max = u_max
            
            try:
                # Root finding: find u where cumulative_arc_length(u) - target = 0
                u_values[i] = brentq(
                    lambda u: cumulative_arc_length(u) - target,
                    u_search_min, u_search_max,
                    xtol=1e-8
                )
            except ValueError:
                # If root finding fails, fall back to linear interpolation
                u_values[i] = u_min + (u_max - u_min) * target / total_length
    
    # Evaluate the spline at the computed parameter values
    if parametric:
        points = np.column_stack([splev(u_values, [t, c[i], k]) for i in range(ndim)])
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

    Examples
    --------
    >>> import numpy as np
    >>> from scipy.interpolate import splprep
    >>> from scipy.interpolate import spline_arclength
    >>> # Create a quarter circle
    >>> theta = np.linspace(0, np.pi/2, 10)
    >>> x = np.cos(theta)
    >>> y = np.sin(theta)
    >>> tck, u = splprep([x, y], s=0)
    >>> length = spline_arclength(tck)
    >>> print(f"Arc length: {length:.4f}")  # Should be close to pi/2
    """
    t, c, k = tck
    
    # Determine if parametric
    try:
        c[0][0]
        parametric = True
        ndim = len(c)
    except (TypeError, IndexError):
        parametric = False
        ndim = 1
    
    # Determine parameter range
    if u_range is None:
        u_min = t[k]
        u_max = t[-k-1]
    else:
        u_min, u_max = u_range
    
    # Arc length derivative
    def arc_length_derivative(u):
        if parametric:
            derivs = [splev(u, [t, c[i], k], der=1) for i in range(ndim)]
            return np.sqrt(sum(d**2 for d in derivs))
        else:
            deriv = splev(u, tck, der=1)
            return np.abs(deriv)
    
    # Integrate to get total length
    length, _ = quad(arc_length_derivative, u_min, u_max, limit=100)
    return length
