"""
Binary search helper for PPoly, BPoly & BSpline.

"""


@cython.wraparound(False)
@cython.boundscheck(False)
@cython.cdivision(True)
cdef int find_interval(double[::1] x,
                       double xval,
                       int prev_interval=0,
                       int extrapolate=1) nogil:
    """
    Find an interval such that x[interval] <= xval < x[interval+1]
    or interval == 0 and xval < x[0]
    or interval == n-2 and xval > x[n-1]

    Parameters
    ----------
    x : array of double, shape (m,)
        Piecewise polynomial breakpoints
    xval : double
        Point to find
    prev_interval : int, optional
        Interval where a previous point was found
    extrapolate : int, optional
        Whether to return the last of the first interval if the
        point is ouf-of-bounds. 

    Returns
    -------
    interval : int
        Suitable interval or -1 if nan.

    """
    cdef int interval, high, low, mid
    cdef double a, b

    a = x[0]
    b = x[x.shape[0]-1]

    interval = prev_interval
    if interval < 0 or interval >= x.shape[0]:
        interval = 0

    if not (a <= xval <= b):
        # Out-of-bounds (or nan)
        if xval < a and extrapolate:
            # below
            interval = 0
        elif xval > b and extrapolate:
            # above
            interval = x.shape[0] - 2
        else:
            # nan or no extrapolation
            interval = -1
    elif xval == b:
        # Make the interval closed from the right
        interval = x.shape[0] - 2
    else:
        # Find the interval the coordinate is in
        # (binary search with locality)
        if xval >= x[interval]:
            low = interval
            high = x.shape[0]-2
        else:
            low = 0
            high = interval

        if xval < x[low+1]:
            high = low

        while low < high:
            mid = (high + low)//2
            if xval < x[mid]:
                # mid < high
                high = mid
            elif xval >= x[mid + 1]:
                low = mid + 1
            else:
                # x[mid] <= xval < x[mid+1]
                low = mid
                break

        interval = low

    return interval
