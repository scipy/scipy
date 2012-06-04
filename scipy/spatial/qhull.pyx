"""
Wrappers for Qhull triangulation, plus some additional N-D geometry utilities

.. versionadded:: 0.9

"""
#
# Copyright (C)  Pauli Virtanen, 2010.
#
# Distributed under the same BSD license as Scipy.
#

import threading
import numpy as np
cimport numpy as np
cimport cython
cimport qhull

__all__ = ['Delaunay', 'tsearch']

#------------------------------------------------------------------------------
# Qhull interface
#------------------------------------------------------------------------------

cdef extern from "stdio.h":
    extern void *stdin
    extern void *stderr
    extern void *stdout

cdef extern from "math.h":
    double fabs(double x) nogil
    double sqrt(double x) nogil

cdef extern from "qhull/src/qset.h":
    ctypedef union setelemT:
        void *p
        int i

    ctypedef struct setT:
        int maxsize
        setelemT e[1]

cdef extern from "qhull/src/qhull.h":
    ctypedef double realT
    ctypedef double coordT
    ctypedef double pointT
    ctypedef int boolT
    ctypedef unsigned int flagT

    ctypedef struct facetT:
        coordT offset
        coordT *center
        coordT *normal
        facetT *next
        facetT *previous
        unsigned id
        setT *vertices
        setT *neighbors
        flagT simplicial
        flagT flipped
        flagT upperdelaunay

    ctypedef struct vertexT:
        vertexT *next
        vertexT *previous
        unsigned int id, visitid
        pointT *point
        setT *neighbours

    ctypedef struct qhT:
        boolT DELAUNAY
        boolT SCALElast
        boolT KEEPcoplanar
        boolT MERGEexact
        boolT NOerrexit
        boolT PROJECTdelaunay
        boolT ATinfinity
        int normal_size
        char *qhull_command
        facetT *facet_list
        facetT *facet_tail
        int num_facets
        unsigned int facet_id
        pointT *first_point
        pointT *input_points
        realT last_low
        realT last_high
        realT last_newhigh
        realT max_outside
        realT MINoutside
        realT DISTround


    extern qhT qh_qh
    extern int qh_PRINToff
    extern int qh_ALL

    void qh_init_A(void *inp, void *out, void *err, int argc, char **argv) nogil
    void qh_init_B(realT *points, int numpoints, int dim, boolT ismalloc) nogil
    void qh_checkflags(char *, char *) nogil
    void qh_initflags(char *) nogil
    void qh_option(char *, char*, char* ) nogil
    void qh_freeqhull(boolT) nogil
    void qh_memfreeshort(int *curlong, int *totlong) nogil
    void qh_qhull() nogil
    void qh_check_output() nogil
    void qh_produce_output() nogil
    void qh_triangulate() nogil
    void qh_checkpolygon() nogil
    void qh_findgood_all() nogil
    void qh_appendprint(int format) nogil
    realT *qh_readpoints(int* num, int *dim, boolT* ismalloc) nogil
    int qh_new_qhull(int dim, int numpoints, realT *points,
                     boolT ismalloc, char* qhull_cmd, void *outfile,
                     void *errfile) nogil
    int qh_pointid(pointT *point) nogil

# Qhull is not threadsafe: needs locking
_qhull_lock = threading.Lock()


#------------------------------------------------------------------------------
# LAPACK interface
#------------------------------------------------------------------------------

cdef extern from "qhull_blas.h":
    void qh_dgetrf(int *m, int *n, double *a, int *lda, int *ipiv,
                   int *info) nogil
    void qh_dgetrs(char *trans, int *n, int *nrhs, double *a, int *lda,
                   int *ipiv, double *b, int *ldb, int *info) nogil
    void qh_dgecon(char *norm, int *n, double *a, int *lda, double *anorm,
                   double *rcond, double *work, int *iwork, int *info) nogil


#------------------------------------------------------------------------------
# Delaunay triangulation using Qhull
#------------------------------------------------------------------------------

def _construct_delaunay(np.ndarray[np.double_t, ndim=2] points):
    """
    Perform Delaunay triangulation of the given set of points.

    """

    # Run qhull with the options
    #
    # - d: perform delaunay triangulation
    # - Qbb: scale last coordinate for Delaunay
    # - Qz: reduces Delaunay precision errors for cospherical sites
    # - Qt: output only simplical facets (can produce degenerate 0-area ones)
    #
    cdef char *options = "qhull d Qz Qbb Qt"
    cdef int curlong, totlong
    cdef int dim
    cdef int numpoints
    cdef int exitcode

    points = np.ascontiguousarray(points)
    numpoints = points.shape[0]
    dim = points.shape[1]

    if numpoints <= 0:
        raise ValueError("No points to triangulate")

    if dim < 2:
        raise ValueError("Need at least 2-D data to triangulate")

    _qhull_lock.acquire()
    try:
        qh_qh.NOerrexit = 1
        with nogil:
            exitcode = qh_new_qhull(dim, numpoints, <realT*>points.data, 0,
                                    options, NULL, stderr)

        try:
            if exitcode != 0:
                raise RuntimeError("Qhull error")

            with nogil:
                qh_triangulate() # get rid of non-simplical facets

            if qh_qh.SCALElast:
                paraboloid_scale = qh_qh.last_newhigh / (
                    qh_qh.last_high - qh_qh.last_low)
                paraboloid_shift = - qh_qh.last_low * paraboloid_scale
            else:
                paraboloid_scale = 1.0
                paraboloid_shift = 0.0

            vertices, neighbors, equations = \
                      _qhull_get_facet_array(dim, numpoints)

            return (vertices, neighbors, equations,
                    paraboloid_scale, paraboloid_shift)
        finally:
            with nogil:
                qh_freeqhull(0)
                qh_memfreeshort(&curlong, &totlong)
            if curlong != 0 or totlong != 0:
                raise RuntimeError("qhull: did not free %d bytes (%d pieces)" %
                                   (totlong, curlong))
    finally:
        _qhull_lock.release()


@cython.boundscheck(False)
@cython.cdivision(True)
def _qhull_get_facet_array(int ndim, int numpoints):
    """
    Return array of simplical facets currently in Qhull.

    Returns
    -------
    vertices : array of int, shape (nfacets, ndim+1)
        Indices of coordinates of vertices forming the simplical facets
    neighbors : array of int, shape (nfacets, ndim)
        Indices of neighboring facets.  The kth neighbor is opposite
        the kth vertex, and the first neighbor is the horizon facet
        for the first vertex.

        Facets extending to infinity are denoted with index -1.

    """

    cdef facetT* facet
    cdef facetT* neighbor
    cdef vertexT *vertex
    cdef int i, j, point, error_non_simplical
    cdef np.ndarray[np.npy_int, ndim=2] vertices
    cdef np.ndarray[np.npy_int, ndim=2] neighbors
    cdef np.ndarray[np.double_t, ndim=2] equations
    cdef np.ndarray[np.npy_int, ndim=1] id_map

    id_map = np.empty((qh_qh.facet_id,), dtype=np.intc)
    id_map.fill(-1)

    # Compute facet indices
    facet = qh_qh.facet_list
    j = 0
    while facet and facet.next:
        if facet.simplicial and not facet.upperdelaunay:
            id_map[facet.id] = j
            j += 1
        facet = facet.next

    # Allocate output
    vertices = np.zeros((j, ndim+1), dtype=np.intc)
    neighbors = np.zeros((j, ndim+1), dtype=np.intc)
    equations = np.zeros((j, ndim+2), dtype=np.double)

    # Retrieve facet information
    error_non_simplical = 0

    with nogil:
        facet = qh_qh.facet_list
        j = 0
        while facet and facet.next:
            if not facet.simplicial:
                error_non_simplical = 1
                break

            if facet.upperdelaunay:
                facet = facet.next
                continue

            # Save vertex info
            for i in xrange(ndim+1):
                vertex = <vertexT*>facet.vertices.e[i].p
                point = qh_pointid(vertex.point)
                vertices[j, i] = point

            # Save neighbor info
            for i in xrange(ndim+1):
                neighbor = <facetT*>facet.neighbors.e[i].p
                neighbors[j,i] = id_map[neighbor.id]

            # Save simplex equation info
            for i in xrange(ndim+1):
                equations[j,i] = facet.normal[i]
            equations[j,ndim+1] = facet.offset

            j += 1
            facet = facet.next

    if error_non_simplical:
        raise ValueError("non-simplical facet encountered")

    return vertices, neighbors, equations


#------------------------------------------------------------------------------
# Barycentric coordinates
#------------------------------------------------------------------------------

@cython.boundscheck(False)
@cython.cdivision(True)
def _get_barycentric_transforms(np.ndarray[np.double_t, ndim=2] points,
                                np.ndarray[np.npy_int, ndim=2] vertices,
                                double eps):
    """
    Compute barycentric affine coordinate transformations for given
    simplices.

    Returns
    -------
    Tinvs : array, shape (nsimplex, ndim+1, ndim)
        Barycentric transforms for each simplex.

        Tinvs[i,:ndim,:ndim] contains inverse of the matrix ``T``,
        and Tinvs[i,ndim,:] contains the vector ``r_n`` (see below).

    Notes
    -----
    Barycentric transform from ``x`` to ``c`` is defined by::

        T c = x - r_n

    where the ``r_1, ..., r_n`` are the vertices of the simplex.
    The matrix ``T`` is defined by the condition::

        T e_j = r_j - r_n

    where ``e_j`` is the unit axis vector, e.g, ``e_2 = [0,1,0,0,...]``
    This implies that ``T_ij = (r_j - r_n)_i``.

    For the barycentric transforms, we need to compute the inverse
    matrix ``T^-1`` and store the vectors ``r_n`` for each vertex.
    These are stacked into the `Tinvs` returned.

    """
    cdef np.ndarray[np.double_t, ndim=2] T
    cdef np.ndarray[np.double_t, ndim=3] Tinvs
    cdef int isimplex
    cdef int i, j, n, nrhs, lda, ldb, info
    cdef int ipiv[NPY_MAXDIMS+1]
    cdef int ndim, nsimplex
    cdef double centroid[NPY_MAXDIMS], c[NPY_MAXDIMS+1]
    cdef double *transform
    cdef double anorm, rcond
    cdef double nan, rcond_limit

    cdef double work[4*NPY_MAXDIMS]
    cdef int iwork[NPY_MAXDIMS]

    cdef double x1, x2, x3
    cdef double y1, y2, y3
    cdef double det

    nan = np.nan
    ndim = points.shape[1]
    nsimplex = vertices.shape[0]

    T = np.zeros((ndim, ndim), dtype=np.double)
    Tinvs = np.zeros((nsimplex, ndim+1, ndim), dtype=np.double)

    # Maximum inverse condition number to allow: we want at least three
    # of the digits be significant, to be safe
    rcond_limit = 1000*eps

    with nogil:
        for isimplex in xrange(nsimplex):
            for i in xrange(ndim):
                Tinvs[isimplex,ndim,i] = points[vertices[isimplex,ndim],i]
                for j in xrange(ndim):
                    T[i,j] = (points[vertices[isimplex,j],i]
                              - Tinvs[isimplex,ndim,i])
                Tinvs[isimplex,i,i] = 1

            # compute 1-norm for estimating condition number
            anorm = _matrix_norm1(ndim, <double*>T.data)

            # LU decomposition
            n = ndim
            nrhs = ndim
            lda = ndim
            ldb = ndim
            qh_dgetrf(&n, &n, <double*>T.data, &lda, ipiv, &info)

            # Check condition number
            if info == 0:
                qh_dgecon("1", &n, <double*>T.data, &lda, &anorm, &rcond,
                          work, iwork, &info)

                if rcond < rcond_limit:
                    # The transform seems singular
                    info = 1

            # Compute transform
            if info == 0:
                qh_dgetrs("N", &n, &nrhs, <double*>T.data, &lda, ipiv,
                          (<double*>Tinvs.data) + ndim*(ndim+1)*isimplex,
                          &ldb, &info)

            # Deal with degenerate simplices
            if info != 0:
                for i in range(ndim+1):
                    for j in range(ndim):
                        Tinvs[isimplex,i,j] = nan

    return Tinvs

@cython.boundscheck(False)
cdef double _matrix_norm1(int n, double *a) nogil:
    """Compute the 1-norm of a square matrix given in in Fortran order"""
    cdef double maxsum = 0, colsum
    cdef int i, j

    for j in range(n):
        colsum = 0
        for i in range(n):
            colsum += fabs(a[0])
            a += 1
        if maxsum < colsum:
            maxsum = colsum
    return maxsum

cdef int _barycentric_inside(int ndim, double *transform,
                             double *x, double *c, double eps) nogil:
    """
    Check whether point is inside a simplex, using barycentric
    coordinates.  `c` will be filled with barycentric coordinates, if
    the point happens to be inside.

    """
    cdef int i, j
    c[ndim] = 1.0
    for i in xrange(ndim):
        c[i] = 0
        for j in xrange(ndim):
            c[i] += transform[ndim*i + j] * (x[j] - transform[ndim*ndim + j])
        c[ndim] -= c[i]

        if not (-eps <= c[i] <= 1 + eps):
            return 0
    if not (-eps <= c[ndim] <= 1 + eps):
        return 0
    return 1

cdef void _barycentric_coordinate_single(int ndim, double *transform,
                                         double *x, double *c, int i) nogil:
    """
    Compute a single barycentric coordinate.

    Before the ndim+1'th coordinate can be computed, the other must have
    been computed earlier.

    """
    cdef int j

    if i == ndim:
        c[ndim] = 1.0
        for j in xrange(ndim):
            c[ndim] -= c[j]
    else:
        c[i] = 0
        for j in xrange(ndim):
            c[i] += transform[ndim*i + j] * (x[j] - transform[ndim*ndim + j])

cdef void _barycentric_coordinates(int ndim, double *transform,
                                   double *x, double *c) nogil:
    """
    Compute barycentric coordinates.

    """
    cdef int i, j
    c[ndim] = 1.0
    for i in xrange(ndim):
        c[i] = 0
        for j in xrange(ndim):
            c[i] += transform[ndim*i + j] * (x[j] - transform[ndim*ndim + j])
        c[ndim] -= c[i]


#------------------------------------------------------------------------------
# N-D geometry
#------------------------------------------------------------------------------

cdef void _lift_point(DelaunayInfo_t *d, double *x, double *z) nogil:
    cdef int i
    z[d.ndim] = 0
    for i in xrange(d.ndim):
        z[i] = x[i]
        z[d.ndim] += x[i]**2
    z[d.ndim] *= d.paraboloid_scale
    z[d.ndim] += d.paraboloid_shift

cdef double _distplane(DelaunayInfo_t *d, int isimplex, double *point) nogil:
    """
    qh_distplane
    """
    cdef double dist
    cdef int k
    dist = d.equations[isimplex*(d.ndim+2) + d.ndim+1]
    for k in xrange(d.ndim+1):
        dist += d.equations[isimplex*(d.ndim+2) + k] * point[k]
    return dist


#------------------------------------------------------------------------------
# Iterating over ridges connected to a vertex in 2D
#------------------------------------------------------------------------------

cdef void _RidgeIter2D_init(RidgeIter2D_t *it, DelaunayInfo_t *d,
                            int vertex) nogil:
    """
    Start iteration over all triangles connected to the given vertex.

    """

    cdef double c[3]
    cdef int k, ivertex, start

    start = 0
    it.info = d
    it.vertex = vertex
    it.triangle = d.vertex_to_simplex[vertex]
    it.start_triangle = it.triangle
    it.restart = 0

    if it.triangle != -1:
        # find some edge connected to this vertex
        for k in xrange(3):
            ivertex = it.info.vertices[it.triangle*3 + k]
            if ivertex != vertex:
                it.vertex2 = ivertex
                it.index = k
                it.start_index = k
                break
    else:
        it.start_index = -1
        it.index = -1

cdef void _RidgeIter2D_next(RidgeIter2D_t *it) nogil:
    cdef int itri, k, ivertex

    #
    # Remember: k-th edge and k-th neigbour are opposite vertex k;
    #           imagine now we are iterating around vertex `O`
    #
    #         .O------,
    #       ./ |\.    |
    #      ./  | \.   |
    #      \   |  \.  |
    #       \  |k  \. |
    #        \ |    \.|
    #         `+------k
    #

    if it.restart:
        if it.start_index == -1:
            # we already did that -> we have iterated over everything
            it.index = -1
            return

        # restart to opposite direction
        it.triangle = it.start_triangle
        for k in xrange(3):
            ivertex = it.info.vertices[it.triangle*3 + k]
            if ivertex != it.vertex and k != it.start_index:
                it.index = k
                it.vertex2 = ivertex
                break
        it.start_index = -1
        it.restart = 0

        if it.info.neighbors[it.triangle*3 + it.index] == -1:
            it.index = -1
            return
        else:
            _RidgeIter2D_next(it)
            if it.index == -1:
                return

    # jump to the next triangle
    itri = it.info.neighbors[it.triangle*3 + it.index]

    # if it's outside triangulation, take the last edge, and signal
    # restart to the opposite direction
    if itri == -1:
        for k in xrange(3):
            ivertex = it.info.vertices[it.triangle*3 + k]
            if ivertex != it.vertex and k != it.index:
                it.index = k
                it.vertex2 = ivertex
                break

        it.restart = 1
        return

    # Find at which index we are now:
    #
    # it.vertex
    #      O-------k------.
    #      | \-          /
    #      |   \- E  B  /
    #      |     \-    /
    #      | A     \- /
    #      +---------´
    #
    # A = it.triangle
    # B = itri
    # E = it.index
    # O = it.vertex
    #
    for k in xrange(3):
        ivertex = it.info.vertices[itri*3 + k]
        if it.info.neighbors[itri*3 + k] != it.triangle and \
               ivertex != it.vertex:
            it.index = k
            it.vertex2 = ivertex
            break

    it.triangle = itri

    # check termination
    if it.triangle == it.start_triangle:
        it.index = -1
        return

cdef class RidgeIter2D(object):
    cdef RidgeIter2D_t it
    cdef object delaunay
    cdef DelaunayInfo_t info

    def __init__(self, delaunay, ivertex):
        if delaunay.ndim != 2:
            raise ValueError("RidgeIter2D supports only 2-D")
        self.delaunay = delaunay
        _get_delaunay_info(&self.info, delaunay, 0, 1)
        _RidgeIter2D_init(&self.it, &self.info, ivertex)

    def __iter__(self):
        return self

    def __next__(self):
        if self.it.index == -1:
            raise StopIteration()
        ret = (self.it.vertex, self.it.vertex2, self.it.index, self.it.triangle)
        _RidgeIter2D_next(&self.it)
        return ret


#------------------------------------------------------------------------------
# Finding simplices
#------------------------------------------------------------------------------

cdef int _is_point_fully_outside(DelaunayInfo_t *d, double *x,
                                 double eps) nogil:
    """
    Is the point outside the bounding box of the triangulation?

    """

    cdef int i
    for i in xrange(d.ndim):
        if x[i] < d.min_bound[i] - eps or x[i] > d.max_bound[i] + eps:
            return 1
    return 0

cdef int _find_simplex_bruteforce(DelaunayInfo_t *d, double *c,
                                  double *x, double eps,
                                  double eps_broad) nogil:
    """
    Find simplex containing point `x` by going through all simplices.

    """
    cdef int inside, isimplex
    cdef int k, m, ineighbor, iself
    cdef double *transform

    if _is_point_fully_outside(d, x, eps):
        return -1

    for isimplex in xrange(d.nsimplex):
        transform = d.transform + isimplex*d.ndim*(d.ndim+1)

        if transform[0] == transform[0]:
            # transform is valid (non-nan)
            inside = _barycentric_inside(d.ndim, transform, x, c, eps)
            if inside:
                return isimplex
        else:
            # transform is invalid (nan, implying degenerate simplex)

            # we replace this inside-check by a check of the neighbors
            # with a larger epsilon

            for k in xrange(d.ndim+1):
                ineighbor = d.neighbors[(d.ndim+1)*isimplex + k]
                if ineighbor == -1:
                    continue

                transform = d.transform + ineighbor*d.ndim*(d.ndim+1)
                if transform[0] != transform[0]:
                    # another bad simplex
                    continue

                _barycentric_coordinates(d.ndim, transform, x, c)

                # Check that the point lies (almost) inside the
                # neigbor simplex
                inside = 1
                for m in xrange(d.ndim+1):
                    if d.neighbors[(d.ndim+1)*ineighbor + m] == isimplex:
                        # allow extra leeway towards isimplex
                        if not (-eps_broad <= c[m] <= 1 + eps):
                            inside = 0
                            break
                    else:
                        # normal check
                        if not (-eps <= c[m] <= 1 + eps):
                            inside = 0
                            break
                if inside:
                    return ineighbor
    return -1

cdef int _find_simplex_directed(DelaunayInfo_t *d, double *c,
                                double *x, int *start, double eps,
                                double eps_broad) nogil:
    """
    Find simplex containing point `x` via a directed walk in the tesselation.

    If the simplex is found, the array `c` is filled with the corresponding
    barycentric coordinates.

    Notes
    -----

    The idea here is the following:

    1) In a simplex, the k-th neighbour is opposite the k-th vertex.
       Call the ridge between them the k-th ridge.

    2) If the k-th barycentric coordinate of the target point is negative,
       then the k-th vertex and the target point lie on the opposite sides
       of the k-th ridge.

    3) Consequently, the k-th neighbour simplex is *closer* to the target point
       than the present simplex, if projected on the normal of the k-th ridge.

    4) In a regular tesselation, hopping to any such direction is OK.

       Also, if one of the negative-coordinate neighbors happens to be -1,
       then the target point is outside the tesselation (because the
       tesselation is convex!).

    5) If all barycentric coordinates are in [-eps, 1+eps], we have found the
       simplex containing the target point.

    6) If all barycentric coordinates are non-negative but 5) is not true,
       we are in an inconsistent situation -- this should never happen.

    This may however enter an infinite loop due to rounding errors in
    the computation of the barycentric coordinates, so the iteration
    count needs to be limited, and a fallback to brute force provided.

    """
    cdef int k, m, ndim, inside, isimplex, cycle_k
    cdef double *transform
    cdef double v

    ndim = d.ndim
    isimplex = start[0]

    if isimplex < 0 or isimplex >= d.nsimplex:
        isimplex = 0

    # The maximum iteration count: it should be large enough so that
    # the algorithm usually succeeds, but smaller than nsimplex so
    # that for the cases where the algorithm fails, the main cost
    # still comes from the brute force search.

    for cycle_k in range(1 + d.nsimplex//4):
        if isimplex == -1:
            break

        transform = d.transform + isimplex*ndim*(ndim+1)

        inside = 1
        for k in xrange(ndim+1):
            _barycentric_coordinate_single(ndim, transform, x, c, k)

            if c[k] < -eps:
                # The target point is in the direction of neighbor `k`!
                m = d.neighbors[(ndim+1)*isimplex + k]
                if m == -1:
                    # The point is outside the triangulation: bail out
                    start[0] = isimplex
                    return -1

                isimplex = m
                inside = -1
                break
            elif c[k] <= 1 + eps:
                # we're inside this simplex
                pass
            else:
                # we're outside (or the coordinate is nan; a degenerate simplex)
                inside = 0

        if inside == -1:
            # hopped to another simplex
            continue
        elif inside == 1:
            # we've found the right one!
            break
        else:
            # we've failed utterly (degenerate simplices in the way).
            # fall back to brute force
            isimplex = _find_simplex_bruteforce(d, c, x, eps, eps_broad)
            break
    else:
        # the algorithm failed to converge -- fall back to brute force
        isimplex = _find_simplex_bruteforce(d, c, x, eps, eps_broad)

    start[0] = isimplex
    return isimplex

cdef int _find_simplex(DelaunayInfo_t *d, double *c,
                       double *x, int *start, double eps,
                       double eps_broad) nogil:
    """
    Find simplex containing point `x` by walking the triangulation.

    Notes
    -----
    This algorithm is similar as used by ``qh_findbest``.  The idea
    is the following:

    1. Delaunay triangulation is a projection of the lower half of a convex
       hull, of points lifted on a paraboloid.

       Simplices in the triangulation == facets on the convex hull.

    2. If a point belongs to a given simplex in the triangulation,
       its image on the paraboloid is on the positive side of
       the corresponding facet.

    3. However, it is not necessarily the *only* such facet.

    4. Also, it is not necessarily the facet whose hyperplane distance
       to the point on the paraboloid is the largest.

    ..note::

        If I'm not mistaken, `qh_findbestfacet` finds a facet for
        which the plane distance is maximized -- so it doesn't always
        return the simplex containing the point given. For example:

        >>> p = np.array([(1 - 1e-4, 0.1)])
        >>> points = np.array([(0,0), (1, 1), (1, 0), (0.99189033, 0.37674127),
        ...                    (0.99440079, 0.45182168)], dtype=np.double)
        >>> tri = qhull.delaunay(points)
        >>> tri.vertices
        array([[4, 1, 0],
               [4, 2, 1],
               [3, 2, 0],
               [3, 4, 0],
               [3, 4, 2]])
        >>> dist = qhull.plane_distance(tri, p)
        >>> dist
        array([[-0.12231439,  0.00184863,  0.01049659, -0.04714842,
                0.00425905]])
        >>> tri.vertices[dist.argmax()]
        array([3, 2, 0]

        Now, the maximally positive-distant simplex is [3, 2, 0], although
        the simplex containing the point is [4, 2, 1].

    In this algorithm, we walk around the tesselation trying to locate
    a positive-distant facet. After finding one, we fall back to a
    directed search.

    """
    cdef int isimplex, i, j, k, inside, ineigh, neighbor_found
    cdef int ndim
    cdef double z[NPY_MAXDIMS+1]
    cdef double best_dist, dist
    cdef int changed

    if _is_point_fully_outside(d, x, eps):
        return -1
    if d.nsimplex <= 0:
        return -1

    ndim = d.ndim
    isimplex = start[0]

    if isimplex < 0 or isimplex >= d.nsimplex:
        isimplex = 0

    # Lift point to paraboloid
    _lift_point(d, x, z)

    # Walk the tesselation searching for a facet with a positive planar distance
    best_dist = _distplane(d, isimplex, z)
    changed = 1
    while changed:
        if best_dist > 0:
            break
        changed = 0
        for k in xrange(ndim+1):
            ineigh = d.neighbors[(ndim+1)*isimplex + k]
            if ineigh == -1:
                continue
            dist = _distplane(d, ineigh, z)

            # Note addition of eps -- otherwise, this code does not
            # necessarily terminate! The compiler may use extended
            # accuracy of the FPU so that (dist > best_dist), but
            # after storing to double size, dist == best_dist,
            # resulting to non-terminating loop

            if dist > best_dist + eps*(1 + fabs(best_dist)):
                # Note: this is intentional: we jump in the middle of the cycle,
                #       and continue the cycle from the next k.
                #
                #       This apparently sweeps the different directions more
                #       efficiently. We don't need full accuracy, since we do
                #       a directed search afterwards in any case.
                isimplex = ineigh
                best_dist = dist
                changed = 1

    # We should now be somewhere near the simplex containing the point,
    # locate it with a directed search
    start[0] = isimplex
    return _find_simplex_directed(d, c, x, start, eps, eps_broad)


#------------------------------------------------------------------------------
# Delaunay triangulation interface, for Python
#------------------------------------------------------------------------------

class Delaunay(object):
    """
    Delaunay(points)

    Delaunay tesselation in N dimensions.

    Parameters
    ----------
    points : ndarray of floats, shape (npoints, ndim)
        Coordinates of points to triangulate

    Attributes
    ----------
    points : ndarray of double, shape (npoints, ndim)
        Points in the triangulation.
    vertices : ndarray of ints, shape (nsimplex, ndim+1)
        Indices of vertices forming simplices in the triangulation.
    neighbors : ndarray of ints, shape (nsimplex, ndim+1)
        Indices of neighbor simplices for each simplex.
        The kth neighbor is opposite to the kth vertex.
        For simplices at the boundary, -1 denotes no neighbor.
    equations : ndarray of double, shape (nsimplex, ndim+2)
        [normal, offset] forming the hyperplane equation of the facet
        on the paraboloid (see [Qhull]_ documentation for more).
    paraboloid_scale, paraboloid_shift : float
        Scale and shift for the extra paraboloid dimension
        (see [Qhull]_ documentation for more).
    transform : ndarray of double, shape (nsimplex, ndim+1, ndim)
        Affine transform from ``x`` to the barycentric coordinates ``c``.
        This is defined by::

            T c = x - r

        At vertex ``j``, ``c_j = 1`` and the other coordinates zero.

        For simplex ``i``, ``transform[i,:ndim,:ndim]`` contains
        inverse of the matrix ``T``, and ``transform[i,ndim,:]``
        contains the vector ``r``.
    vertex_to_simplex : ndarray of int, shape (npoints,)
        Lookup array, from a vertex, to some simplex which it is a part of.
    convex_hull : ndarray of int, shape (nfaces, ndim)
        Vertices of facets forming the convex hull of the point set.
        The array contains the indices of the points belonging to
        the (N-1)-dimensional facets that form the convex hull
        of the triangulation.

    Notes
    -----
    The tesselation is computed using the Qhull libary [Qhull]_.

    .. versionadded:: 0.9

    References
    ----------
    .. [Qhull] http://www.qhull.org/

    """
    def __init__(self, points):
        points = np.ascontiguousarray(points).astype(np.double)
        vertices, neighbors, equations, paraboloid_scale, paraboloid_shift = \
                  _construct_delaunay(points)

        self.ndim = points.shape[1]
        self.npoints = points.shape[0]
        self.nsimplex = vertices.shape[0]
        self.points = points
        self.vertices = vertices
        self.neighbors = neighbors
        self.equations = equations
        self.paraboloid_scale = paraboloid_scale
        self.paraboloid_shift = paraboloid_shift
        self.min_bound = self.points.min(axis=0)
        self.max_bound = self.points.max(axis=0)
        self._transform = None
        self._vertex_to_simplex = None

    @property
    def transform(self):
        """
        Affine transform from ``x`` to the barycentric coordinates ``c``.

        :type: ndarray of double, shape (nsimplex, ndim+1, ndim)

        This is defined by::

            T c = x - r

        At vertex ``j``, ``c_j = 1`` and the other coordinates zero.

        For simplex ``i``, ``transform[i,:ndim,:ndim]`` contains
        inverse of the matrix ``T``, and ``transform[i,ndim,:]``
        contains the vector ``r``.

        """
        if self._transform is None:
            self._transform = _get_barycentric_transforms(self.points,
                                                          self.vertices,
                                                          np.finfo(float).eps)
        return self._transform

    @property
    @cython.boundscheck(False)
    def vertex_to_simplex(self):
        """
        Lookup array, from a vertex, to some simplex which it is a part of.

        :type: ndarray of int, shape (npoints,)
        """
        cdef int isimplex, k, ivertex, nsimplex, ndim
        cdef np.ndarray[np.npy_int, ndim=2] vertices
        cdef np.ndarray[np.npy_int, ndim=1] arr

        if self._vertex_to_simplex is None:
            self._vertex_to_simplex = np.empty((self.npoints,), dtype=np.intc)
            self._vertex_to_simplex.fill(-1)

            arr = self._vertex_to_simplex
            vertices = self.vertices

            nsimplex = self.nsimplex
            ndim = self.ndim

            with nogil:
                for isimplex in xrange(nsimplex):
                    for k in xrange(ndim+1):
                        ivertex = vertices[isimplex, k]
                        if arr[ivertex] == -1:
                            arr[ivertex] = isimplex

        return self._vertex_to_simplex

    @property
    @cython.boundscheck(False)
    def convex_hull(self):
        """
        Vertices of facets forming the convex hull of the point set.

        :type: ndarray of int, shape (nfaces, ndim)

        The array contains the indices of the points
        belonging to the (N-1)-dimensional facets that form the convex
        hull of the triangulation.

        """
        cdef int isimplex, k, j, ndim, nsimplex, m, msize
        cdef object out
        cdef np.ndarray[np.npy_int, ndim=2] arr
        cdef np.ndarray[np.npy_int, ndim=2] neighbors
        cdef np.ndarray[np.npy_int, ndim=2] vertices

        neighbors = self.neighbors
        vertices = self.vertices
        ndim = self.ndim
        nsimplex = self.nsimplex

        msize = 10
        out = np.empty((msize, ndim), dtype=np.intc)
        arr = out

        m = 0
        for isimplex in xrange(nsimplex):
            for k in xrange(ndim+1):
                if neighbors[isimplex,k] == -1:
                    for j in xrange(ndim+1):
                        if j < k:
                            arr[m,j] = vertices[isimplex,j]
                        elif j > k:
                            arr[m,j-1] = vertices[isimplex,j]
                    m += 1

                    if m >= msize:
                        arr = None
                        msize = 2*msize + 1
                        out.resize(msize, ndim)
                        arr = out

        arr = None
        try:
            out.resize(m, ndim)
        except ValueError:
            # XXX: work around a Cython bug on Python 2.4
            #      still leaks memory, though
            return np.resize(out, (m, ndim))
        return out

    @cython.boundscheck(False)
    def find_simplex(self, xi, bruteforce=False, tol=None):
        """
        find_simplex(self, xi, bruteforce=False, tol=None)

        Find the simplices containing the given points.

        Parameters
        ----------
        tri : DelaunayInfo
            Delaunay triangulation
        xi : ndarray of double, shape (..., ndim)
            Points to locate
        bruteforce : bool, optional
            Whether to only perform a brute-force search
        tol : float, optional
            Tolerance allowed in the inside-triangle check.
            Default is ``100*eps``.

        Returns
        -------
        i : ndarray of int, same shape as `xi`
            Indices of simplices containing each point.
            Points outside the triangulation get the value -1.

        Notes
        -----
        This uses an algorithm adapted from Qhull's ``qh_findbestfacet``,
        which makes use of the connection between a convex hull and a
        Delaunay triangulation. After finding the simplex closest to
        the point in N+1 dimensions, the algorithm falls back to
        directed search in N dimensions.

        """
        cdef DelaunayInfo_t info
        cdef int isimplex
        cdef double c[NPY_MAXDIMS]
        cdef double eps, eps_broad
        cdef int start
        cdef int k
        cdef np.ndarray[np.double_t, ndim=2] x
        cdef np.ndarray[np.npy_int, ndim=1] out_

        xi = np.asanyarray(xi)

        if xi.shape[-1] != self.ndim:
            raise ValueError("wrong dimensionality in xi")

        xi_shape = xi.shape
        xi = xi.reshape(np.prod(xi.shape[:-1]), xi.shape[-1])
        x = np.ascontiguousarray(xi.astype(np.double))

        start = 0

        if tol is None:
            eps = 100 * np.finfo(np.double).eps
        else:
            eps = tol
        eps_broad = np.sqrt(eps)
        out = np.zeros((xi.shape[0],), dtype=np.intc)
        out_ = out
        _get_delaunay_info(&info, self, 1, 0)

        if bruteforce:
            with nogil:
                for k in xrange(x.shape[0]):
                    isimplex = _find_simplex_bruteforce(
                        &info, c,
                        <double*>x.data + info.ndim*k,
                        eps, eps_broad)
                    out_[k] = isimplex
        else:
            with nogil:
                for k in xrange(x.shape[0]):
                    isimplex = _find_simplex(&info, c,
                                             <double*>x.data + info.ndim*k,
                                             &start, eps, eps_broad)
                    out_[k] = isimplex

        return out.reshape(xi_shape[:-1])

    @cython.boundscheck(False)
    def plane_distance(self, xi):
        """
        plane_distance(self, xi)

        Compute hyperplane distances to the point `xi` from all simplices.

        """
        cdef np.ndarray[np.double_t, ndim=2] x
        cdef np.ndarray[np.double_t, ndim=2] out_
        cdef DelaunayInfo_t info
        cdef double z[NPY_MAXDIMS+1]
        cdef int i, j, k

        if xi.shape[-1] != self.ndim:
            raise ValueError("xi has different dimensionality than "
                             "triangulation")

        xi_shape = xi.shape
        xi = xi.reshape(np.prod(xi.shape[:-1]), xi.shape[-1])
        x = np.ascontiguousarray(xi.astype(np.double))

        _get_delaunay_info(&info, self, 0, 0)

        out = np.zeros((x.shape[0], info.nsimplex), dtype=np.double)
        out_ = out

        with nogil:
            for i in xrange(x.shape[0]):
                for j in xrange(info.nsimplex):
                    _lift_point(&info, (<double*>x.data) + info.ndim*i, z)
                    out_[i,j] = _distplane(&info, j, z)

        return out.reshape(xi_shape[:-1] + (self.nsimplex,))

    def lift_points(self, x):
        """
        lift_points(self, x)

        Lift points to the Qhull paraboloid.

        """
        z = np.zeros(x.shape[:-1] + (x.shape[-1]+1,), dtype=np.double)
        z[...,:-1] = x
        z[...,-1] = (x**2).sum(axis=-1)
        z[...,-1] *= self.paraboloid_scale
        z[...,-1] += self.paraboloid_shift
        return z

# Alias familiar from other environments
def tsearch(tri, xi):
    """
    tsearch(tri, xi)

    Find simplices containing the given points. This function does the
    same thing as `Delaunay.find_simplex`.

    .. versionadded:: 0.9

    See Also
    --------
    Delaunay.find_simplex

    """
    return tri.find_simplex(xi)


#------------------------------------------------------------------------------
# Delaunay triangulation interface, for low-level C
#------------------------------------------------------------------------------

cdef int _get_delaunay_info(DelaunayInfo_t *info,
                            obj,
                            int compute_transform,
                            int compute_vertex_to_simplex) except -1:
    cdef np.ndarray[np.double_t, ndim=3] transform
    cdef np.ndarray[np.npy_int, ndim=1] vertex_to_simplex
    cdef np.ndarray[np.double_t, ndim=2] points = obj.points
    cdef np.ndarray[np.npy_int, ndim=2] vertices = obj.vertices
    cdef np.ndarray[np.npy_int, ndim=2] neighbors = obj.neighbors
    cdef np.ndarray[np.double_t, ndim=2] equations = obj.equations
    cdef np.ndarray[np.double_t, ndim=1] min_bound = obj.min_bound
    cdef np.ndarray[np.double_t, ndim=1] max_bound = obj.max_bound

    info.ndim = points.shape[1]
    info.npoints = points.shape[0]
    info.nsimplex = vertices.shape[0]
    info.points = <double*>points.data
    info.vertices = <int*>vertices.data
    info.neighbors = <int*>neighbors.data
    info.equations = <double*>equations.data
    info.paraboloid_scale = obj.paraboloid_scale
    info.paraboloid_shift = obj.paraboloid_shift
    if compute_transform:
        transform = obj.transform
        info.transform = <double*>transform.data
    else:
        info.transform = NULL
    if compute_vertex_to_simplex:
        vertex_to_simplex = obj.vertex_to_simplex
        info.vertex_to_simplex = <int*>vertex_to_simplex.data
    else:
        info.vertex_to_simplex = NULL
    info.min_bound = <double*>min_bound.data
    info.max_bound = <double*>max_bound.data

    return 0
