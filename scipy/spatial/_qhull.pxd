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
                             double *x, double *c, double eps) nogil

cdef void _barycentric_coordinate_single(int ndim, double *transform,
                                         double *x, double *c, int i) nogil

cdef void _barycentric_coordinates(int ndim, double *transform,
                                   double *x, double *c) nogil

#
# N+1-D geometry
#

cdef void _lift_point(DelaunayInfo_t *d, double *x, double *z) nogil

cdef double _distplane(DelaunayInfo_t *d, int isimplex, double *point) nogil

#
# Finding simplices
#

cdef int _is_point_fully_outside(DelaunayInfo_t *d, double *x, double eps) nogil

cdef int _find_simplex_bruteforce(DelaunayInfo_t *d, double *c, double *x,
                                  double eps, double eps_broad) nogil

cdef int _find_simplex_directed(DelaunayInfo_t *d, double *c, double *x,
                                int *start, double eps, double eps_broad) nogil

cdef int _find_simplex(DelaunayInfo_t *d, double *c, double *x, int *start,
                       double eps, double eps_broad) nogil
