import math
import numpy as np

#pythran export _compute_prob_inside_method(int, int, int, int)
def _compute_prob_inside_method(m, n, g, h):
    """
    Count the proportion of paths that stay strictly inside two diagonal lines.

    Parameters
    ----------
    m : integer
        m > 0
    n : integer
        n > 0
    g : integer
        g is greatest common divisor of m and n
    h : integer
        0 <= h <= lcm(m,n)

    Returns
    -------
    p : float
        The proportion of paths that stay inside the two lines.

    Count the integer lattice paths from (0, 0) to (m, n) which satisfy
    |x/m - y/n| < h / lcm(m, n).
    The paths make steps of size +1 in either positive x or positive y
    directions.

    We generally follow Hodges' treatment of Drion/Gnedenko/Korolyuk.
    Hodges, J.L. Jr.,
    "The Significance Probability of the Smirnov Two-Sample Test,"
    Arkiv fiur Matematik, 3, No. 43 (1958), 469-86.

    """
    # Probability is symmetrical in m, n.  Computation below uses m >= n.
    if m < n:
        m, n = n, m
    mg = m // g
    ng = n // g

    # Count the integer lattice paths from (0, 0) to (m, n) which satisfy
    # |nx/g - my/g| < h.
    # Compute matrix A such that:
    #  A(x, 0) = A(0, y) = 1
    #  A(x, y) = A(x, y-1) + A(x-1, y), for x,y>=1, except that
    #  A(x, y) = 0 if |x/m - y/n|>= h
    # Probability is A(m, n)/binom(m+n, n)
    # Optimizations exist for m==n, m==n*p.
    # Only need to preserve a single column of A, and only a
    # sliding window of it.
    # minj keeps track of the slide.
    minj, maxj = 0, min(int(np.ceil(h / mg)), n + 1)
    curlen = maxj - minj
    # Make a vector long enough to hold maximum window needed.
    lenA = min(2 * maxj + 2, n + 1)
    # This is an integer calculation, but the entries are essentially
    # binomial coefficients, hence grow quickly.
    # Scaling after each column is computed avoids dividing by a
    # large binomial coefficent at the end, but is not sufficient to avoid
    # the large dyanamic range which appears during the calculation.
    # Instead we rescale based on the magnitude of the right most term in
    # the column and keep track of an exponent separately and apply
    # it at the end of the calculation.  Similarly when multiplying by
    # the binomial coefficint
    dtype = np.float64
    A = np.zeros(lenA, dtype=dtype)
    # Initialize the first column
    A[minj:maxj] = 1
    expnt = 0
    for i in range(1, m + 1):
        # Generate the next column.
        # First calculate the sliding window
        lastminj, lastlen = minj, curlen
        minj = max(int(np.floor((ng * i - h) / mg)) + 1, 0)
        minj = min(minj, n)
        maxj = min(int(np.ceil((ng * i + h) / mg)), n + 1)
        if maxj <= minj:
            return 0
        # Now fill in the values
        A[0:maxj - minj] = np.cumsum(A[minj - lastminj:maxj - lastminj])
        curlen = maxj - minj
        if lastlen > curlen:
            # Set some carried-over elements to 0
            A[maxj - minj:maxj - minj + (lastlen - curlen)] = 0
        # Rescale if the right most value is over 2**900
        val = A[maxj - minj - 1]
        _, valexpt = math.frexp(val)
        if valexpt > 900:
            # Scaling to bring down to about 2**800 appears
            # sufficient for sizes under 10000.
            valexpt -= 800
            A = np.ldexp(A, -valexpt)
            expnt += valexpt

    val = A[maxj - minj - 1]
    # Now divide by the binomial (m+n)!/m!/n!
    for i in range(1, n + 1):
        val = (val * i) / (m + i)
        _, valexpt = math.frexp(val)
        if valexpt < -128:
            val = np.ldexp(val, -valexpt)
            expnt += valexpt
    # Finally scale if needed.
    return np.ldexp(val, expnt)

# pythran export _compute_prob_outside_square(int, int)
def _compute_prob_outside_square(n, h):
    """
    Compute the proportion of paths that pass outside the two diagonal lines.

    Parameters
    ----------
    n : integer
        n > 0
    h : integer
        0 <= h <= n

    Returns
    -------
    p : float
        The proportion of paths that pass outside the lines x-y = +/-h.

    """
    # Compute Pr(D_{n,n} >= h/n)
    # Prob = 2 * ( binom(2n, n-h) - binom(2n, n-2a) + binom(2n, n-3a) - ... )
    # / binom(2n, n)
    # This formulation exhibits subtractive cancellation.
    # Instead divide each term by binom(2n, n), then factor common terms
    # and use a Horner-like algorithm
    # P = 2 * A0 * (1 - A1*(1 - A2*(1 - A3*(1 - A4*(...)))))
    P = 0.0
    k = int(np.floor(n / h))
    while k >= 0:
        p1 = 1.0
        # Each of the Ai terms has numerator and denominator with
        # h simple terms.
        for j in range(h):
            p1 = (n - k * h - j) * p1 / (n + k * h + j + 1)
        P = p1 * (1.0 - P)
        k -= 1
    return 2 * P

