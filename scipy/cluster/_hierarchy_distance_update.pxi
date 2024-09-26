"""
A `linkage_distance_update` function calculates the distance from cluster i
to the new cluster xy after merging cluster x and cluster y

Parameters
----------
d_xi : double
    Distance from cluster x to cluster i
d_yi : double
    Distance from cluster y to cluster i
d_xy : double
    Distance from cluster x to cluster y
size_x : int
    Size of cluster x
size_y : int
    Size of cluster y
size_i : int
    Size of cluster i

Returns
-------
d_xyi : double
    Distance from the new cluster xy to cluster i
"""
ctypedef double (*linkage_distance_update)(double d_xi, double d_yi,
                                           double d_xy, int size_x,
                                           int size_y, int size_i) noexcept


cdef double _single(double d_xi, double d_yi, double d_xy,
                    int size_x, int size_y, int size_i) noexcept:
    return min(d_xi, d_yi)


cdef double _complete(double d_xi, double d_yi, double d_xy,
                      int size_x, int size_y, int size_i) noexcept:
    return max(d_xi, d_yi)


cdef double _average(double d_xi, double d_yi, double d_xy,
                     int size_x, int size_y, int size_i) noexcept:
    return (size_x * d_xi + size_y * d_yi) / (size_x + size_y)


cdef double _centroid(double d_xi, double d_yi, double d_xy,
                      int size_x, int size_y, int size_i) noexcept:
    return sqrt((((size_x * d_xi * d_xi) + (size_y * d_yi * d_yi)) -
                 (size_x * size_y * d_xy * d_xy) / (size_x + size_y)) /
                (size_x + size_y))


cdef double _median(double d_xi, double d_yi, double d_xy,
                    int size_x, int size_y, int size_i) noexcept:
    return sqrt(0.5 * (d_xi * d_xi + d_yi * d_yi) - 0.25 * d_xy * d_xy)


cdef double _ward(double d_xi, double d_yi, double d_xy,
                  int size_x, int size_y, int size_i) noexcept:
    cdef double t = 1.0 / (size_x + size_y + size_i)
    return sqrt((size_i + size_x) * t * d_xi * d_xi +
                (size_i + size_y) * t * d_yi * d_yi -
                size_i * t * d_xy * d_xy)


cdef double _weighted(double d_xi, double d_yi, double d_xy,
                      int size_x, int size_y, int size_i) noexcept:
    return 0.5 * (d_xi + d_yi)
