# -*-cython-*-
"""
Qhull shared definitions, for use by other Cython modules

"""
#
# Copyright (C)  Pauli Virtanen, 2010.
#
# Distributed under the same BSD license as Scipy.
#

cdef extern from "numpy/ndarrayobject.h":
    cdef enum:
        NPY_MAXDIMS

ctypedef struct DelaunayInfo_t:
    int ndim
    int npoints
    int nsimplex
    double *points
    int *simplices
    int *neighbors
    double *equations
    double *transform
    int *vertex_to_simplex
    double paraboloid_scale
    double paraboloid_shift
    double *max_bound
    double *min_bound
    int *vertex_neighbors_indices
    int *vertex_neighbors_indptr

cdef int _get_delaunay_info(DelaunayInfo_t *, obj,
                            int compute_transform,
                            int compute_vertex_to_simplex,
                            int compute_vertex_neighbors) except -1

#
# N-D geometry
#

cdef int _barycentric_inside(int ndim, double *transform,
                             const double *x, double *c, double eps) noexcept nogil

cdef void _barycentric_coordinate_single(int ndim, double *transform,
                                         const double *x, double *c, int i) noexcept nogil

cdef void _barycentric_coordinates(int ndim, double *transform,
                                   const double *x, double *c) noexcept nogil

#
# N+1-D geometry
#

cdef void _lift_point(DelaunayInfo_t *d, const double *x, double *z) noexcept nogil

cdef double _distplane(DelaunayInfo_t *d, int isimplex, double *point) noexcept nogil

#
# Finding simplices
#

cdef int _is_point_fully_outside(DelaunayInfo_t *d, const double *x, double eps) noexcept nogil

cdef int _find_simplex_bruteforce(DelaunayInfo_t *d, double *c, const double *x,
                                  double eps, double eps_broad) noexcept nogil

cdef int _find_simplex_directed(DelaunayInfo_t *d, double *c, const double *x,
                                int *start, double eps, double eps_broad) noexcept nogil

cdef int _find_simplex(DelaunayInfo_t *d, double *c, const double *x, int *start,
                       double eps, double eps_broad) noexcept nogil
