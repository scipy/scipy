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
cimport setlist

from numpy.compat import asbytes

cdef extern from "numpy/npy_math.h":
    double nan "NPY_NAN"

__all__ = ['Delaunay', 'ConvexHull', 'Voronoi', 'tsearch']

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

cdef extern from "setjmp.h" nogil:
    ctypedef struct jmp_buf:
        pass
    int setjmp(jmp_buf STATE) nogil
    void longjmp(jmp_buf STATE, int VALUE) nogil

# Define the clockwise constant
cdef extern from "qhull/src/user.h":
    cdef enum:
        qh_ORIENTclock

cdef extern from "qhull/src/qset.h":
    ctypedef union setelemT:
        void *p
        int i

    ctypedef struct setT:
        int maxsize
        setelemT e[1]

    int qh_setsize(setT *set) nogil
    void qh_setappend(setT **setp, void *elem) nogil

cdef extern from "qhull/src/libqhull.h":
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
        setT *ridges
        setT *coplanarset
        flagT simplicial
        flagT flipped
        flagT upperdelaunay
        flagT toporient
        unsigned visitid

    ctypedef struct vertexT:
        vertexT *next
        vertexT *previous
        unsigned int id, visitid
        pointT *point
        setT *neighbors

    ctypedef struct ridgeT:
        setT *vertices
        facetT *top
        facetT *bottom

    ctypedef struct qhT:
        boolT DELAUNAY
        boolT SCALElast
        boolT KEEPcoplanar
        boolT MERGEexact
        boolT NOerrexit
        boolT PROJECTdelaunay
        boolT ATinfinity
        boolT UPPERdelaunay
        boolT hasTriangulation
        int normal_size
        char *qhull_command
        facetT *facet_list
        facetT *facet_tail
        vertexT *vertex_list
        vertexT *vertex_tail
        int num_facets
        int num_points
        int center_size
        unsigned int facet_id
        pointT *first_point
        pointT *input_points
        realT last_low
        realT last_high
        realT last_newhigh
        realT max_outside
        realT MINoutside
        realT DISTround
        realT totvol
        realT totarea
        jmp_buf errexit
        setT *other_points
        unsigned int visit_id
        unsigned int vertex_visit

    extern qhT *qh_qh
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
    setT *qh_pointvertex() nogil
    realT *qh_readpoints(int* num, int *dim, boolT* ismalloc) nogil
    int qh_new_qhull(int dim, int numpoints, realT *points,
                     boolT ismalloc, char* qhull_cmd, void *outfile,
                     void *errfile) nogil
    int qh_pointid(pointT *point) nogil
    vertexT *qh_nearvertex(facetT *facet, pointT *point, double *dist) nogil
    boolT qh_addpoint(pointT *furthest, facetT *facet, boolT checkdist) nogil
    facetT *qh_findbestfacet(pointT *point, boolT bestoutside,
                             realT *bestdist, boolT *isoutside) nogil
    void qh_setdelaunay(int dim, int count, pointT *points) nogil
    void qh_restore_qhull(qhT **oldqh) nogil
    qhT *qh_save_qhull() nogil

cdef extern from "qhull/src/io.h":
    ctypedef enum qh_RIDGE:
        qh_RIDGEall
        qh_RIDGEinner
        qh_RIDGEouter

    ctypedef void printvridgeT(void *fp, vertexT *vertex, vertexT *vertexA,
                               setT *centers, boolT unbounded)
    int qh_eachvoronoi_all(void *fp, void* printvridge,
                           boolT isUpper, qh_RIDGE innerouter,
                           boolT inorder) nogil

    void qh_order_vertexneighbors(vertexT *vertex) nogil
    int qh_compare_facetvisit(void *p1, void *p2) nogil

cdef extern from "qhull/src/geom.h":
    pointT *qh_facetcenter(setT *vertices) nogil

cdef extern from "qhull/src/geom.h":
    double qh_getarea(facetT *facetlist) nogil

cdef extern from "qhull/src/poly.h":
    void qh_check_maxout() nogil

cdef extern from "qhull/src/mem.h":
    void qh_memfree(void *object, int insize)

from libc.string cimport memcpy
from libc.stdlib cimport qsort

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
# Qhull wrapper
#------------------------------------------------------------------------------

# Qhull is not threadsafe: needs locking
_qhull_lock = threading.Lock()

# Qhull has (swappable) global state: keep track which Qhull instance is active
# and how many instances are alive
cdef _Qhull _active_qhull = None
cdef int _qhull_count = 0

# Qhull objects pending cleanup
#
# Python's garbage collector can trigger a call to a destructor while
# the qhull lock is held.  Destructors, for instance that of
# Voronoi/etc, can call _Qhull methods that require the lock, which
# causes a deadlock also in a single-threaded code.
#
# We ensure that _Qhull.close is safe to call from a destructor, by
# postponing the cleanup if the lock happens to be held. The other
# methods are not safe to call.
#
cdef list _qhull_pending_cleanup = []

class QhullError(RuntimeError):
    pass

@cython.final
cdef class _Qhull:
    cdef qhT *_saved_qh
    cdef list _point_arrays
    cdef public bytes options
    cdef public bytes mode_option
    cdef public object furthest_site

    cdef readonly int ndim
    cdef int numpoints, _is_delaunay
    cdef np.ndarray _ridge_points

    cdef list _ridge_vertices
    cdef object _ridge_error
    cdef int _nridges

    cdef np.ndarray _ridge_equations

    @cython.final
    def __init__(self,
                 bytes mode_option,
                 np.ndarray[np.double_t, ndim=2] points,
                 bytes options=None,
                 bytes required_options=None,
                 furthest_site=False,
                 incremental=False):
        global _active_qhull, _qhull_count
        cdef int exitcode

        points = np.ascontiguousarray(points, dtype=np.double)

        self.numpoints = points.shape[0]
        self.ndim = points.shape[1]

        if self.numpoints <= 0:
            raise ValueError("No points given")
        if self.ndim < 2:
            raise ValueError("Need at least 2-D data")
        if np.isnan(points).any():
            raise ValueError("Points cannot contain NaN")

        # Process options
        option_set = set()
        if options is not None:
            option_set.update(options.split())
        if furthest_site:
            if b"Qz" in option_set:
                option_set.remove(b"Qz")
            option_set.add(b"Qu")
        if required_options is not None:
            required_option_set = set(required_options.split())
            if b"QJ" in option_set and b"Qt" in required_option_set:
                # safe to remove, QJ always produces simplical output
                required_option_set.remove(b"Qt")
            option_set.update(required_option_set)

        if incremental:
            incremental_bad_ops = set([b'Qbb', b'Qbk', b'QBk', b'QbB', b'Qz'])
            bad_opts = []
            for bad_opt in incremental_bad_ops:
                if bad_opt in options:
                    bad_opts.append(bad_opt)
            if bad_opts:
                raise ValueError("Qhull options %r are incompatible "
                                 "with incremental mode" % (bad_opts,))

            if b"Qt" in option_set:
                # Qhull wants this
                option_set.add(b"Q11")

            # We need to own the copy of the points in incremental mode
            points = points.copy()

        if mode_option in (b"d", b"v"):
            self._is_delaunay = 1
        else:
            self._is_delaunay = 0

        self._point_arrays = [points]
        self.options = b" ".join(option_set)
        self.mode_option = mode_option
        self.furthest_site = furthest_site

        options = b"qhull "  + mode_option +  b" " + self.options

        _qhull_lock.acquire()
        try:
            if _active_qhull is not None:
                _active_qhull._deactivate()

            _active_qhull = self
            _qhull_count += 1

            options_c = <char*>options
            with nogil:
                exitcode = qh_new_qhull(self.ndim, self.numpoints,
                                        <realT*>points.data, 0,
                                        options_c, NULL, stderr)

            if exitcode != 0:
                self._uninit()
                raise QhullError("Qhull error")
        finally:
            _qhull_lock.release()

    @cython.final
    def volume_area(self):
        cdef double volume
        cdef double area

        _qhull_lock.acquire()
        try:
            self._activate()
            qh_getarea(qh_qh.facet_list)
            volume = qh_qh.totvol
            area = qh_qh.totarea
        finally:
            _qhull_lock.release()

        return volume, area

    @cython.final
    def close(self):
        if _qhull_lock.acquire(False):
            try:
                self._cleanup_pending()
                self._uninit()
            finally:
                _qhull_lock.release()
        else:
            # Failed to acquire the lock. 
            _qhull_pending_cleanup.append(self)

    @cython.final
    cdef int _cleanup_pending(self) except -1:
        """
        Process any pending cleanups (_qhull_lock MUST be held when calling this)
        """
        cdef _Qhull qh
        cdef int k

        for k in range(len(_qhull_pending_cleanup)):
            try:
                qh = _qhull_pending_cleanup.pop()
            except IndexError:
                break
            qh._uninit()

        return 0

    @cython.final
    cdef int _activate(self) except -1:
        """
        Activate this instance (_qhull_lock MUST be held when calling this)
        """
        global _active_qhull

        if _active_qhull is self:
            return 0
        elif _active_qhull is not None:
            _active_qhull._deactivate()

        assert _active_qhull is None

        if self._saved_qh == NULL:
            raise RuntimeError("Qhull instance is closed")

        qh_restore_qhull(&self._saved_qh)
        self._saved_qh = NULL
        _active_qhull = self

        return 0

    @cython.final
    cdef int _deactivate(self) except -1:
        """
        Deactivate this instance (_qhull_lock MUST be held when calling this)
        """
        global _active_qhull

        assert _active_qhull is self
        assert self._saved_qh == NULL

        self._saved_qh = qh_save_qhull()
        _active_qhull = None

    @cython.final
    cdef int _uninit(self) except -1:
        """
        Uninitialize this instance (_qhull_lock MUST be held when calling this)
        """
        global _active_qhull, _qhull_count
        cdef int curlong, totlong

        if not (_active_qhull is self or self._saved_qh != NULL):
            # already freed
            return 0

        self._activate()

        qh_freeqhull(qh_ALL)

        _qhull_count -= 1
        _active_qhull = None
        self._saved_qh = NULL

        if _qhull_count == 0:
            # last one out cleans the house
            qh_memfreeshort(&curlong, &totlong)
            if curlong != 0 or totlong != 0:
                raise QhullError(
                    "qhull: did not free %d bytes (%d pieces)" %
                    (totlong, curlong))

        return 0

    @cython.final
    def get_points(self):
        if len(self._point_arrays) == 1:
            return self._point_arrays[0]
        else:
            return np.concatenate(
                [x[:,:self.ndim] for x in self._point_arrays],
                axis=0)

    @cython.final
    def add_points(self, points):
        cdef int j
        cdef realT *p
        cdef facetT *facet
        cdef double bestdist
        cdef boolT isoutside
        cdef np.ndarray arr

        points = np.asarray(points)
        if points.ndim!=2 or points.shape[1] != self._point_arrays[0].shape[1]:
            raise ValueError("invalid size for new points array")
        if points.size == 0:
            return

        if self._is_delaunay:
            arr = np.empty((points.shape[0], self.ndim+1), dtype=np.double)
            arr[:,:-1] = points
        else:
            arr = np.array(points, dtype=np.double, order="C", copy=True)

        _qhull_lock.acquire()
        try:
            self._activate()

            # nonlocal error handling
            exitcode = setjmp(qh_qh.errexit)
            if exitcode != 0:
                raise QhullError("Qhull error")
            qh_qh.NOerrexit = 0

            # add points to triangulation
            if self._is_delaunay:
                # lift to paraboloid
                qh_setdelaunay(arr.shape[1], arr.shape[0], <realT*>arr.data)

            p = <realT*>arr.data

            for j in xrange(arr.shape[0]):
                facet = qh_findbestfacet(p, 0, &bestdist, &isoutside)
                if isoutside:
                    if not qh_addpoint(p, facet, 0):
                        break
                else:
                    # append the point to the "other points" list, to
                    # maintain the point IDs
                    qh_setappend(&qh_qh.other_points, p)

                p += arr.shape[1]

            qh_check_maxout()
            qh_qh.hasTriangulation = 0

            self._point_arrays.append(arr)
            self.numpoints += arr.shape[0]
        finally:
            qh_qh.NOerrexit = 1
            _qhull_lock.release()

    @cython.final
    def get_paraboloid_shift_scale(self):
        cdef double paraboloid_scale
        cdef double paraboloid_shift

        _qhull_lock.acquire()
        try:
            self._activate()

            if qh_qh.SCALElast:
                paraboloid_scale = qh_qh.last_newhigh / (
                    qh_qh.last_high - qh_qh.last_low)
                paraboloid_shift = - qh_qh.last_low * paraboloid_scale
            else:
                paraboloid_scale = 1.0
                paraboloid_shift = 0.0
        finally:
            _qhull_lock.release()

        return paraboloid_scale, paraboloid_shift

    @cython.final
    def triangulate(self):
        _qhull_lock.acquire()
        try:
            self._activate()

            with nogil:
                qh_triangulate() # get rid of non-simplical facets
        finally:
            _qhull_lock.release()

    @cython.final
    def get_simplex_facet_array(self):
        _qhull_lock.acquire()
        try:
            self._activate()
            return self._get_simplex_facet_array()
        finally:
            _qhull_lock.release()

    @cython.final
    @cython.boundscheck(False)
    @cython.cdivision(True)
    cdef _get_simplex_facet_array(self):
        """
        Return array of simplical facets currently in Qhull.

        Returns
        -------
        facets : array of int, shape (nfacets, ndim+1)
            Indices of coordinates of vertices forming the simplical facets
        neighbors : array of int, shape (nfacets, ndim)
            Indices of neighboring facets.  The kth neighbor is opposite
            the kth vertex, and the first neighbor is the horizon facet
            for the first vertex.

            Facets extending to infinity are denoted with index -1.
        equations : array of double, shape (nfacets, ndim+2)

        """
        cdef facetT* facet
        cdef facetT* neighbor
        cdef vertexT *vertex
        cdef pointT *point
        cdef int i, j, ipoint, ipoint2, ncoplanar
        cdef object tmp
        cdef np.ndarray[np.npy_int, ndim=2] facets
        cdef np.ndarray[np.npy_int, ndim=2] neighbors
        cdef np.ndarray[np.npy_int, ndim=2] coplanar
        cdef np.ndarray[np.double_t, ndim=2] equations
        cdef np.ndarray[np.npy_int, ndim=1] id_map
        cdef double dist
        cdef int facet_ndim
        cdef int numpoints
        cdef unsigned int lower_bound
        cdef unsigned int swapped_index

        facet_ndim = self.ndim
        numpoints = self.numpoints

        if self._is_delaunay:
            facet_ndim += 1

        id_map = np.empty(qh_qh.facet_id, dtype=np.intc)

        # Compute facet indices
        with nogil:
            for i in range(qh_qh.facet_id):
                id_map[i] = -1

            facet = qh_qh.facet_list
            j = 0
            while facet and facet.next:
                if not self._is_delaunay or facet.upperdelaunay == qh_qh.UPPERdelaunay:
                    if not facet.simplicial and ( \
                           qh_setsize(facet.vertices) != facet_ndim or \
                           qh_setsize(facet.neighbors) != facet_ndim):
                        with gil:
                            raise QhullError(
                                "non-simplical facet encountered: %r vertices"
                                % (qh_setsize(facet.vertices),))

                    id_map[facet.id] = j
                    j += 1

                facet = facet.next

        # Allocate output
        facets = np.zeros((j, facet_ndim), dtype=np.intc)
        neighbors = np.zeros((j, facet_ndim), dtype=np.intc)
        equations = np.zeros((j, facet_ndim+1), dtype=np.double)

        ncoplanar = 0
        coplanar = np.zeros((10, 3), dtype=np.intc)

        # Retrieve facet information
        with nogil:
            facet = qh_qh.facet_list
            j = 0
            while facet and facet.next:
                if self._is_delaunay and facet.upperdelaunay != qh_qh.UPPERdelaunay:
                    facet = facet.next
                    continue

                # Use a lower bound so that the tight loop in high dimensions
                # is not affected by the conditional below
                lower_bound = 0
                if (self._is_delaunay and
                    facet.toporient == qh_ORIENTclock and facet_ndim == 3):
                    # Swap the first and second indices to maintain a
                    # counter-clockwise orientation.
                    for i in xrange(2):
                        # Save the vertex info
                        swapped_index = 1 ^ i
                        vertex = <vertexT*>facet.vertices.e[i].p
                        ipoint = qh_pointid(vertex.point)
                        facets[j, swapped_index] = ipoint

                        # Save the neighbor info
                        neighbor = <facetT*>facet.neighbors.e[i].p
                        neighbors[j, swapped_index] = id_map[neighbor.id]

                    lower_bound = 2

                for i in xrange(lower_bound, facet_ndim):
                    # Save the vertex info
                    vertex = <vertexT*>facet.vertices.e[i].p
                    ipoint = qh_pointid(vertex.point)
                    facets[j, i] = ipoint

                    # Save the neighbor info
                    neighbor = <facetT*>facet.neighbors.e[i].p
                    neighbors[j, i] = id_map[neighbor.id]

                # Save simplex equation info
                for i in xrange(facet_ndim):
                    equations[j, i] = facet.normal[i]
                equations[j, facet_ndim] = facet.offset

                # Save coplanar info
                if facet.coplanarset:
                    for i in range(qh_setsize(facet.coplanarset)):
                        point = <pointT*>facet.coplanarset.e[i].p
                        vertex = qh_nearvertex(facet, point, &dist)

                        if ncoplanar >= coplanar.shape[0]:
                            with gil:
                                tmp = coplanar
                                coplanar = None
                                try:
                                    tmp.resize(2 * ncoplanar + 1, 3)
                                except ValueError:
                                    # Work around Cython issue on Python 2.4
                                    tmp = np.resize(tmp, (2*ncoplanar+1, 3))
                                coplanar = tmp

                        coplanar[ncoplanar, 0] = qh_pointid(point)
                        coplanar[ncoplanar, 1] = id_map[facet.id]
                        coplanar[ncoplanar, 2] = qh_pointid(vertex.point)
                        ncoplanar += 1

                j += 1
                facet = facet.next

        return facets, neighbors, equations, coplanar[:ncoplanar]

    @cython.final
    def get_voronoi_diagram(_Qhull self):
        _qhull_lock.acquire()
        try:
            self._activate()
            return self._get_voronoi_diagram()
        finally:
            _qhull_lock.release()

    @cython.final
    @cython.boundscheck(False)
    @cython.cdivision(True)
    cdef _get_voronoi_diagram(_Qhull self):
        """
        Return the voronoi diagram currently in Qhull.

        Returns
        -------
        voronoi_vertices : array of double, shape (nvoronoi_vertices, ndim)
            Coordinates of the Voronoi vertices

        ridge_points : array of double, shape (nridges, 2)
            Voronoi ridges, as indices to the points array.

        ridge_vertices : list of lists, shape (nridges, *)
            Voronoi vertices for each Voronoi ridge, as indices to
            the Voronoi vertices array.
            Infinity is indicated by index ``-1``.

        regions : list of lists, shape (nregion, *)
            Voronoi vertices of all regions.

        point_region : array of int, shape (npoint,)
            Index of the Voronoi region for each input point.

        """
        cdef int i, j, k
        cdef vertexT *vertex
        cdef facetT *neighbor
        cdef facetT *facet

        cdef object tmp
        cdef np.ndarray[np.double_t, ndim=2] voronoi_vertices
        cdef np.ndarray[np.intp_t, ndim=1] point_region
        cdef int nvoronoi_vertices
        cdef pointT infty_point[NPY_MAXDIMS+1]
        cdef pointT *point
        cdef pointT *center
        cdef double dist
        cdef int inf_seen

        cdef list regions
        cdef list cur_region

        # -- Grab Voronoi ridges
        self._nridges = 0
        self._ridge_error = None
        self._ridge_points = np.empty((10, 2), np.intc)
        self._ridge_vertices = []

        qh_eachvoronoi_all(<void*>self, &_visit_voronoi, qh_qh.UPPERdelaunay,
                           qh_RIDGEall, 1)

        self._ridge_points = self._ridge_points[:self._nridges]

        if self._ridge_error is not None:
            raise self._ridge_error

        # Now, qh_eachvoronoi_all has initialized the visitids of facets
        # to correspond do the Voronoi vertex indices.

        # -- Grab Voronoi regions
        regions = []

        point_region = np.empty(self.numpoints, np.intp)
        for i in range(self.numpoints):
            point_region[i] = -1

        vertex = qh_qh.vertex_list
        while vertex and vertex.next:
            qh_order_vertexneighbors_nd(self.ndim+1, vertex)

            i = qh_pointid(vertex.point)
            if i < self.numpoints:
                # Qz results to one extra point
                point_region[i] = len(regions)

            inf_seen = 0
            cur_region = []
            for k in xrange(qh_setsize(vertex.neighbors)):
                neighbor = <facetT*>vertex.neighbors.e[k].p
                i = neighbor.visitid - 1
                if i == -1:
                    if not inf_seen:
                        inf_seen = 1
                    else:
                        continue
                cur_region.append(int(i))
            if len(cur_region) == 1 and cur_region[0] == -1:
                # report similarly as qvoronoi o
                cur_region = []
            regions.append(cur_region)

            vertex = vertex.next

        # -- Grab Voronoi vertices and point-to-region map
        nvoronoi_vertices = 0
        voronoi_vertices = np.empty((10, self.ndim), np.double)

        facet = qh_qh.facet_list
        while facet and facet.next:
            if facet.visitid > 0:
                # finite Voronoi vertex

                center = qh_facetcenter(facet.vertices)

                nvoronoi_vertices = max(facet.visitid, nvoronoi_vertices)
                if nvoronoi_vertices >= voronoi_vertices.shape[0]:
                    tmp = voronoi_vertices
                    voronoi_vertices = None
                    try:
                        tmp.resize(2*nvoronoi_vertices + 1, self.ndim)
                    except ValueError:
                        tmp = np.resize(tmp, (2*nvoronoi_vertices+1, self.ndim))
                    voronoi_vertices = tmp

                for k in range(self.ndim):
                    voronoi_vertices[facet.visitid-1, k] = center[k]

                qh_memfree(center, qh_qh.center_size)

                if facet.coplanarset:
                    for k in range(qh_setsize(facet.coplanarset)):
                        point = <pointT*>facet.coplanarset.e[k].p
                        vertex = qh_nearvertex(facet, point, &dist)

                        i = qh_pointid(point)
                        j = qh_pointid(vertex.point)

                        if i < self.numpoints:
                            # Qz can result to one extra point
                            point_region[i] = point_region[j]

            facet = facet.next

        voronoi_vertices = voronoi_vertices[:nvoronoi_vertices]

        return voronoi_vertices, self._ridge_points, self._ridge_vertices, \
               regions, point_region

    @cython.final
    def get_extremes_2d(_Qhull self):
        if self._is_delaunay:
            raise ValueError("Cannot compute for Delaunay")

        _qhull_lock.acquire()
        try:
            self._activate()
            return self._get_extremes_2d()
        finally:
            _qhull_lock.release()

    @cython.final
    @cython.boundscheck(False)
    @cython.cdivision(True)
    cdef _get_extremes_2d(_Qhull self):
        """
        Compute the extremal points in a 2-D convex hull, i.e. the
        vertices of the convex hull, ordered counterclockwise.

        See qhull/io.c:qh_printextremes_2d

        """
        cdef facetT *facet
        cdef facetT *startfacet
        cdef facetT *nextfacet
        cdef vertexT *vertexA
        cdef vertexT *vertexB
        cdef int[:] extremes
        cdef int nextremes

        nextremes = 0
        extremes_arr = np.zeros(100, dtype=np.intc)
        extremes = extremes_arr

        qh_qh.visit_id += 1
        qh_qh.vertex_visit += 1

        facet = qh_qh.facet_list
        startfacet = facet
        while facet:
            if facet.visitid == qh_qh.visit_id:
                raise QhullError("Qhull internal error: loop in facet list")

            if facet.toporient:
                vertexA = <vertexT*>facet.vertices.e[0].p
                vertexB = <vertexT*>facet.vertices.e[1].p
                nextfacet = <facetT*>facet.neighbors.e[0].p
            else:
                vertexB = <vertexT*>facet.vertices.e[0].p
                vertexA = <vertexT*>facet.vertices.e[1].p
                nextfacet = <facetT*>facet.neighbors.e[1].p

            if nextremes + 2 >= extremes.shape[0]:
                extremes = None
                extremes_arr.resize(2*extremes_arr.shape[0]+1)
                extremes = extremes_arr

            if vertexA.visitid != qh_qh.vertex_visit:
                vertexA.visitid = qh_qh.vertex_visit
                extremes[nextremes] = qh_pointid(vertexA.point)
                nextremes += 1

            if vertexB.visitid != qh_qh.vertex_visit:
                vertexB.visitid = qh_qh.vertex_visit
                extremes[nextremes] = qh_pointid(vertexB.point)
                nextremes += 1

            facet.visitid = qh_qh.visit_id
            facet = nextfacet

            if facet == startfacet:
                break

        extremes = None
        extremes_arr.resize(nextremes)
        return extremes_arr


cdef void _visit_voronoi(void *ptr, vertexT *vertex, vertexT *vertexA,
                         setT *centers, boolT unbounded):
    cdef _Qhull qh = <_Qhull>ptr
    cdef int point_1, point_2, ix
    cdef list cur_vertices

    if qh._ridge_error is not None:
        return

    if qh._nridges >= qh._ridge_points.shape[0]:
        try:
            qh._ridge_points.resize(2*qh._nridges + 1, 2)
        except Exception, e:
            qh._ridge_error = e
            return

    # Record which points the ridge is between
    point_1 = qh_pointid(vertex.point)
    point_2 = qh_pointid(vertexA.point)

    p = <int*>qh._ridge_points.data
    p[2*qh._nridges + 0] = point_1
    p[2*qh._nridges + 1] = point_2

    # Record which voronoi vertices constitute the ridge
    cur_vertices = []
    for i in xrange(qh_setsize(centers)):
        ix = (<facetT*>centers.e[i].p).visitid - 1
        cur_vertices.append(ix)
    qh._ridge_vertices.append(cur_vertices)

    qh._nridges += 1

    return

cdef void qh_order_vertexneighbors_nd(int nd, vertexT *vertex):
    if nd == 3:
        qh_order_vertexneighbors(vertex)
    elif nd >= 4:
        qsort(<facetT**>&vertex.neighbors.e[0].p, qh_setsize(vertex.neighbors),
              sizeof(facetT*), qh_compare_facetvisit)


#------------------------------------------------------------------------------
# Barycentric coordinates
#------------------------------------------------------------------------------

@cython.boundscheck(False)
@cython.cdivision(True)
def _get_barycentric_transforms(np.ndarray[np.double_t, ndim=2] points,
                                np.ndarray[np.npy_int, ndim=2] simplices,
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
    cdef double centroid[NPY_MAXDIMS]
    cdef double c[NPY_MAXDIMS+1]
    cdef double *transform
    cdef double anorm, rcond
    cdef double rcond_limit

    cdef double work[4*NPY_MAXDIMS]
    cdef int iwork[NPY_MAXDIMS]

    cdef double x1, x2, x3
    cdef double y1, y2, y3
    cdef double det

    ndim = points.shape[1]
    nsimplex = simplices.shape[0]

    T = np.zeros((ndim, ndim), dtype=np.double)
    Tinvs = np.zeros((nsimplex, ndim+1, ndim), dtype=np.double)

    # Maximum inverse condition number to allow: we want at least three
    # of the digits be significant, to be safe
    rcond_limit = 1000*eps

    with nogil:
        for isimplex in xrange(nsimplex):
            for i in xrange(ndim):
                Tinvs[isimplex,ndim,i] = points[simplices[isimplex,ndim],i]
                for j in xrange(ndim):
                    T[i,j] = (points[simplices[isimplex,j],i]
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
        >>> tri.simplices
        array([[4, 1, 0],
               [4, 2, 1],
               [3, 2, 0],
               [3, 4, 0],
               [3, 4, 2]])
        >>> dist = qhull.plane_distance(tri, p)
        >>> dist
        array([[-0.12231439,  0.00184863,  0.01049659, -0.04714842,
                0.00425905]])
        >>> tri.simplices[dist.argmax()]
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

class _QhullUser(object):
    """
    Takes care of basic dealings with the Qhull objects
    """

    _qhull = None

    def __init__(self, qhull, incremental=False):
        self._qhull = None
        try:
            self._update(qhull)
            if incremental:
                # last, to deal with exceptions
                self._qhull = qhull
        finally:
            if qhull is not self._qhull:
                qhull.close()

    def close(self):
        """
        close()

        Finish incremental processing.

        Call this to free resources taken up by Qhull, when using the
        incremental mode. After calling this, adding more points is no
        longer possible.
        """
        if self._qhull is not None:
            self._qhull.close()
            self._qhull = None

    def __del__(self):
        self.close()

    def _update(self, qhull):
        self.points = qhull.get_points()
        self.ndim = self.points.shape[1]
        self.npoints = self.points.shape[0]
        self.points = self.points
        self.min_bound = self.points.min(axis=0)
        self.max_bound = self.points.max(axis=0)

    def add_points(self, points, restart=False):
        """
        add_points(points, restart=False)

        Process a set of additional new points.

        Parameters
        ----------
        points : ndarray
            New points to add. The dimensionality should match that of the
            initial points.
        restart : bool, optional
            Whether to restart processing from scratch, rather than
            adding points incrementally.

        Raises
        ------
        QhullError
            Raised when Qhull encounters an error condition, such as
            geometrical degeneracy when options to resolve are not enabled.

        See Also
        --------
        close

        Notes
        -----
        You need to specify ``incremental=True`` when constructing the
        object to be able to add points incrementally. Incremental addition
        of points is also not possible after `close` has been called.

        """
        if self._qhull is None:
            raise RuntimeError("incremental mode not enabled or already closed")

        if restart:
            points = np.concatenate([self.points, points], axis=0)
            qhull = _Qhull(self._qhull.mode_option, points,
                           options=self._qhull.options,
                           furthest_site=self._qhull.furthest_site,
                           incremental=True)
            try:
                self._update(qhull)
                self._qhull = qhull
            finally:
                if qhull is not self._qhull:
                    qhull.close()
            return

        self._qhull.add_points(points)
        self._update(self._qhull)


class Delaunay(_QhullUser):
    """
    Delaunay(points, furthest_site=False, incremental=False, qhull_options=None)

    Delaunay tesselation in N dimensions.

    .. versionadded:: 0.9

    Parameters
    ----------
    points : ndarray of floats, shape (npoints, ndim)
        Coordinates of points to triangulate
    furthest_site : bool, optional
        Whether to compute a furthest-site Delaunay triangulation.
        Default: False

        .. versionadded:: 0.12.0
    incremental : bool, optional
        Allow adding new points incrementally. This takes up some additional
        resources.
    qhull_options : str, optional
        Additional options to pass to Qhull. See Qhull manual for
        details. Option "Qt" is always enabled.
        Default:"Qbb Qc Qz Qx" for ndim > 4 and "Qbb Qc Qz" otherwise.
        Incremental mode omits "Qz".

        .. versionadded:: 0.12.0

    Attributes
    ----------
    points : ndarray of double, shape (npoints, ndim)
        Coordinates of input points.
    simplices : ndarray of ints, shape (nsimplex, ndim+1)
        Indices of the points forming the simplices in the triangulation.
        For 2-D, the points are oriented counterclockwise.
    neighbors : ndarray of ints, shape (nsimplex, ndim+1)
        Indices of neighbor simplices for each simplex.
        The kth neighbor is opposite to the kth vertex.
        For simplices at the boundary, -1 denotes no neighbor.
    equations : ndarray of double, shape (nsimplex, ndim+2)
        [normal, offset] forming the hyperplane equation of the facet
        on the paraboloid
        (see `Qhull documentation <http://www.qhull.org/>`__ for more).
    paraboloid_scale, paraboloid_shift : float
        Scale and shift for the extra paraboloid dimension
        (see `Qhull documentation <http://www.qhull.org/>`__ for more).
    transform : ndarray of double, shape (nsimplex, ndim+1, ndim)
        Affine transform from ``x`` to the barycentric coordinates ``c``.
        This is defined by::

            T c = x - r

        At vertex ``j``, ``c_j = 1`` and the other coordinates zero.

        For simplex ``i``, ``transform[i,:ndim,:ndim]`` contains
        inverse of the matrix ``T``, and ``transform[i,ndim,:]``
        contains the vector ``r``.

        If the simplex is degenerate or nearly degenerate, its
        barycentric transform contains NaNs.
    vertex_to_simplex : ndarray of int, shape (npoints,)
        Lookup array, from a vertex, to some simplex which it is a part of.
        If qhull option "Qc" was not specified, the list will contain -1
        for points that are not vertices of the tesselation.
    convex_hull : ndarray of int, shape (nfaces, ndim)
        Vertices of facets forming the convex hull of the point set.
        The array contains the indices of the points belonging to
        the (N-1)-dimensional facets that form the convex hull
        of the triangulation.

        .. note::

           Computing convex hulls via the Delaunay triangulation is
           inefficient and subject to increased numerical instability.
           Use `ConvexHull` instead.
    coplanar : ndarray of int, shape (ncoplanar, 3)
        Indices of coplanar points and the corresponding indices of
        the nearest facet and the nearest vertex.  Coplanar
        points are input points which were *not* included in the
        triangulation due to numerical precision issues.

        If option "Qc" is not specified, this list is not computed.

        .. versionadded:: 0.12.0
    vertices
        Same as `simplices`, but deprecated.
    vertex_neighbor_vertices : tuple of two ndarrays of int; (indices, indptr)
        Neighboring vertices of vertices. The indices of neighboring
        vertices of vertex `k` are ``indptr[indices[k]:indices[k+1]]``.

    Raises
    ------
    QhullError
        Raised when Qhull encounters an error condition, such as
        geometrical degeneracy when options to resolve are not enabled.
    ValueError
        Raised if an incompatible array is given as input.

    Notes
    -----
    The tesselation is computed using the Qhull library 
    `Qhull library <http://www.qhull.org/>`__.

    .. note::

       Unless you pass in the Qhull option "QJ", Qhull does not
       guarantee that each input point appears as a vertex in the
       Delaunay triangulation. Omitted points are listed in the
       `coplanar` attribute.

    Do not call the ``add_points`` method from a ``__del__``
    destructor.

    Examples
    --------
    Triangulation of a set of points:

    >>> points = np.array([[0, 0], [0, 1.1], [1, 0], [1, 1]])
    >>> from scipy.spatial import Delaunay
    >>> tri = Delaunay(points)

    We can plot it:

    >>> import matplotlib.pyplot as plt
    >>> plt.triplot(points[:,0], points[:,1], tri.simplices.copy())
    >>> plt.plot(points[:,0], points[:,1], 'o')
    >>> plt.show()

    Point indices and coordinates for the two triangles forming the
    triangulation:

    >>> tri.simplices
    array([[2, 3, 0],                 # may vary
           [3, 1, 0]], dtype=int32)

    Note that depending on how rounding errors go, the simplices may
    be in a different order than above.

    >>> points[tri.simplices]
    array([[[ 1. ,  0. ],            # may vary
            [ 1. ,  1. ],
            [ 0. ,  0. ]],
           [[ 1. ,  1. ],
            [ 0. ,  1.1],
            [ 0. ,  0. ]]])

    Triangle 0 is the only neighbor of triangle 1, and it's opposite to
    vertex 1 of triangle 1:

    >>> tri.neighbors[1]
    array([-1,  0, -1], dtype=int32)
    >>> points[tri.simplices[1,1]]
    array([ 0. ,  1.1])

    We can find out which triangle points are in:

    >>> p = np.array([(0.1, 0.2), (1.5, 0.5)])
    >>> tri.find_simplex(p)
    array([ 1, -1], dtype=int32)

    We can also compute barycentric coordinates in triangle 1 for
    these points:

    >>> b = tri.transform[1,:2].dot(p - tri.transform[1,2])
    >>> np.c_[b, 1 - b.sum(axis=1)]
    array([[ 0.1       ,  0.2       ,  0.7       ],
           [ 1.27272727,  0.27272727, -0.54545455]])

    The coordinates for the first point are all positive, meaning it
    is indeed inside the triangle.

    """

    def __init__(self, points, furthest_site=False, incremental=False,
                 qhull_options=None):
        if np.ma.isMaskedArray(points):
            raise ValueError('Input points cannot be a masked array')
        points = np.ascontiguousarray(points, dtype=np.double)

        if qhull_options is None:
            if not incremental:
                qhull_options = b"Qbb Qc Qz"
            else:
                qhull_options = b"Qc"
            if points.shape[1] >= 5:
                qhull_options += b" Qx"
        else:
            qhull_options = asbytes(qhull_options)

        # Run qhull
        qhull = _Qhull(b"d", points, qhull_options, required_options=b"Qt",
                       furthest_site=furthest_site, incremental=incremental)
        _QhullUser.__init__(self, qhull, incremental=incremental)

    def _update(self, qhull):
        qhull.triangulate()

        self.paraboloid_scale, self.paraboloid_shift = \
                               qhull.get_paraboloid_shift_scale()

        self.simplices, self.neighbors, self.equations, self.coplanar = \
                       qhull.get_simplex_facet_array()

        self.nsimplex = self.simplices.shape[0]
        self._transform = None
        self._vertex_to_simplex = None
        self._vertex_neighbor_vertices = None

        # Backwards compatibility (Scipy < 0.12.0)
        self.vertices = self.simplices

        _QhullUser._update(self, qhull)

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
                                                          self.simplices,
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
        cdef np.ndarray[np.npy_int, ndim=2] simplices
        cdef np.ndarray[np.npy_int, ndim=1] arr

        if self._vertex_to_simplex is None:
            self._vertex_to_simplex = np.empty((self.npoints,), dtype=np.intc)
            self._vertex_to_simplex.fill(-1)

            # include coplanar points
            self._vertex_to_simplex[self.coplanar[:,0]] = self.coplanar[:,2]

            # include other points
            arr = self._vertex_to_simplex
            simplices = self.simplices

            coplanar = self.coplanar
            ncoplanar = coplanar.shape[0]

            nsimplex = self.nsimplex
            ndim = self.ndim

            with nogil:
                for isimplex in xrange(nsimplex):
                    for k in xrange(ndim+1):
                        ivertex = simplices[isimplex, k]
                        if arr[ivertex] == -1:
                            arr[ivertex] = isimplex

        return self._vertex_to_simplex

    @property
    @cython.boundscheck(False)
    def vertex_neighbor_vertices(self):
        """
        Neighboring vertices of vertices.

        Tuple of two ndarrays of int: (indices, indptr). The indices of
        neighboring vertices of vertex `k` are
        ``indptr[indices[k]:indices[k+1]]``.

        """
        cdef int i, j, k, m, is_neighbor, is_missing, ndata, idata
        cdef int nsimplex, npoints, ndim
        cdef np.ndarray[np.npy_int, ndim=2] simplices
        cdef setlist.setlist_t sets

        if self._vertex_neighbor_vertices is None:
            ndim = self.ndim
            npoints = self.npoints
            nsimplex = self.nsimplex
            simplices = self.simplices

            setlist.init(&sets, npoints, ndim+1)

            try:
                with nogil:
                    for i in xrange(nsimplex):
                        for j in xrange(ndim+1):
                            for k in xrange(ndim+1):
                                if simplices[i,j] != simplices[i,k]:
                                    if setlist.add(&sets, simplices[i,j], simplices[i,k]):
                                        with gil:
                                            raise MemoryError

                self._vertex_neighbor_vertices = setlist.tocsr(&sets)
            finally:
                setlist.free(&sets)

        return self._vertex_neighbor_vertices

    @property
    @cython.boundscheck(False)
    def convex_hull(self):
        """
        Vertices of facets forming the convex hull of the point set.

        :type: ndarray of int, shape (nfaces, ndim)

        The array contains the indices of the points
        belonging to the (N-1)-dimensional facets that form the convex
        hull of the triangulation.

        .. note::

           Computing convex hulls via the Delaunay triangulation is
           inefficient and subject to increased numerical instability.
           Use `ConvexHull` instead.

        """
        cdef int isimplex, k, j, ndim, nsimplex, m, msize
        cdef object out
        cdef np.ndarray[np.npy_int, ndim=2] arr
        cdef np.ndarray[np.npy_int, ndim=2] neighbors
        cdef np.ndarray[np.npy_int, ndim=2] simplices

        neighbors = self.neighbors
        simplices = self.simplices
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
                            arr[m,j] = simplices[isimplex,j]
                        elif j > k:
                            arr[m,j-1] = simplices[isimplex,j]
                    m += 1

                    if m >= msize:
                        arr = None
                        msize = 2*msize + 1
                        try:
                            out.resize(msize, ndim)
                        except ValueError:
                            # Work around Cython bug on Python 2.4
                            out = np.resize(out, (msize, ndim))
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
        xi = xi.reshape(-1, xi.shape[-1])
        x = np.ascontiguousarray(xi.astype(np.double))

        start = 0

        if tol is None:
            eps = 100 * np.finfo(np.double).eps
        else:
            eps = tol
        eps_broad = sqrt(eps)
        out = np.zeros((xi.shape[0],), dtype=np.intc)
        out_ = out
        _get_delaunay_info(&info, self, 1, 0, 0)

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
        xi = xi.reshape(-1, xi.shape[-1])
        x = np.ascontiguousarray(xi.astype(np.double))

        _get_delaunay_info(&info, self, 0, 0, 0)

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
                            int compute_vertex_to_simplex,
                            int compute_vertex_neighbor_vertices) except -1:
    cdef np.ndarray[np.double_t, ndim=3] transform
    cdef np.ndarray[np.npy_int, ndim=1] vertex_to_simplex
    cdef np.ndarray[np.npy_int, ndim=1] vn_indices, vn_indptr
    cdef np.ndarray[np.double_t, ndim=2] points = obj.points
    cdef np.ndarray[np.npy_int, ndim=2] simplices = obj.simplices
    cdef np.ndarray[np.npy_int, ndim=2] neighbors = obj.neighbors
    cdef np.ndarray[np.double_t, ndim=2] equations = obj.equations
    cdef np.ndarray[np.double_t, ndim=1] min_bound = obj.min_bound
    cdef np.ndarray[np.double_t, ndim=1] max_bound = obj.max_bound

    info.ndim = points.shape[1]
    info.npoints = points.shape[0]
    info.nsimplex = simplices.shape[0]
    info.points = <double*>points.data
    info.simplices = <int*>simplices.data
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
    if compute_vertex_neighbor_vertices:
        vn_indices, vn_indptr = obj.vertex_neighbor_vertices
        info.vertex_neighbors_indices = <int*>vn_indices.data
        info.vertex_neighbors_indptr = <int*>vn_indptr.data
    else:
        info.vertex_neighbors_indices = NULL
        info.vertex_neighbors_indptr = NULL
    info.min_bound = <double*>min_bound.data
    info.max_bound = <double*>max_bound.data

    return 0


#------------------------------------------------------------------------------
# Convex hulls
#------------------------------------------------------------------------------

class ConvexHull(_QhullUser):
    """
    ConvexHull(points, incremental=False, qhull_options=None)

    Convex hulls in N dimensions.

    .. versionadded:: 0.12.0

    Parameters
    ----------
    points : ndarray of floats, shape (npoints, ndim)
        Coordinates of points to construct a convex hull from
    incremental : bool, optional
        Allow adding new points incrementally. This takes up some additional
        resources.
    qhull_options : str, optional
        Additional options to pass to Qhull. See Qhull manual
        for details. (Default: "Qx" for ndim > 4 and "" otherwise)
        Option "Qt" is always enabled.

    Attributes
    ----------
    points : ndarray of double, shape (npoints, ndim)
        Coordinates of input points.
    vertices : ndarray of ints, shape (nvertices,)
        Indices of points forming the vertices of the convex hull.
        For 2-D convex hulls, the vertices are in counterclockwise order.
        For other dimensions, they are in input order.
    simplices : ndarray of ints, shape (nfacet, ndim)
        Indices of points forming the simplical facets of the convex hull.
    neighbors : ndarray of ints, shape (nfacet, ndim)
        Indices of neighbor facets for each facet.
        The kth neighbor is opposite to the kth vertex.
        -1 denotes no neighbor.
    equations : ndarray of double, shape (nfacet, ndim+1)
        [normal, offset] forming the hyperplane equation of the facet
        (see `Qhull documentation <http://www.qhull.org/>`__  for more).
    coplanar : ndarray of int, shape (ncoplanar, 3)
        Indices of coplanar points and the corresponding indices of
        the nearest facets and nearest vertex indices.  Coplanar
        points are input points which were *not* included in the
        triangulation due to numerical precision issues.

        If option "Qc" is not specified, this list is not computed.
    area : float
        Area of the convex hull
    volume : float
        Volume of the convex hull

    Raises
    ------
    QhullError
        Raised when Qhull encounters an error condition, such as
        geometrical degeneracy when options to resolve are not enabled.
    ValueError
        Raised if an incompatible array is given as input.

    Notes
    -----
    The convex hull is computed using the 
    `Qhull library <http://www.qhull.org/>`__.

    Do not call the ``add_points`` method from a ``__del__``
    destructor.

    Examples
    --------

    Convex hull of a random set of points:

    >>> from scipy.spatial import ConvexHull
    >>> points = np.random.rand(30, 2)   # 30 random points in 2-D
    >>> hull = ConvexHull(points)

    Plot it:

    >>> import matplotlib.pyplot as plt
    >>> plt.plot(points[:,0], points[:,1], 'o')
    >>> for simplex in hull.simplices:
    ...     plt.plot(points[simplex, 0], points[simplex, 1], 'k-')

    We could also have directly used the vertices of the hull, which
    for 2-D are guaranteed to be in counterclockwise order:

    >>> plt.plot(points[hull.vertices,0], points[hull.vertices,1], 'r--', lw=2)
    >>> plt.plot(points[hull.vertices[0],0], points[hull.vertices[0],1], 'ro')
    >>> plt.show()

    References
    ----------
    .. [Qhull] http://www.qhull.org/

    """

    def __init__(self, points, incremental=False, qhull_options=None):
        if np.ma.isMaskedArray(points):
            raise ValueError('Input points cannot be a masked array')
        points = np.ascontiguousarray(points, dtype=np.double)

        if qhull_options is None:
            qhull_options = b""
            if points.shape[1] >= 5:
                qhull_options += b"Qx"
        else:
            qhull_options = asbytes(qhull_options)

        # Run qhull
        qhull = _Qhull(b"i", points, qhull_options, required_options=b"Qt",
                       incremental=incremental)
        _QhullUser.__init__(self, qhull, incremental=incremental)

    def _update(self, qhull):
        qhull.triangulate()

        self.simplices, self.neighbors, self.equations, self.coplanar = \
                       qhull.get_simplex_facet_array()

        self.volume, self.area = qhull.volume_area()

        if qhull.ndim == 2:
            self._vertices = qhull.get_extremes_2d()
        else:
            self._vertices = None

        self.nsimplex = self.simplices.shape[0]

        _QhullUser._update(self, qhull)

    @property
    def vertices(self):
        if self._vertices is None:
            self._vertices = np.unique(self.simplices)
        return self._vertices


#------------------------------------------------------------------------------
# Voronoi diagrams
#------------------------------------------------------------------------------

class Voronoi(_QhullUser):
    """
    Voronoi(points, furthest_site=False, incremental=False, qhull_options=None)

    Voronoi diagrams in N dimensions.

    .. versionadded:: 0.12.0

    Parameters
    ----------
    points : ndarray of floats, shape (npoints, ndim)
        Coordinates of points to construct a convex hull from
    furthest_site : bool, optional
        Whether to compute a furthest-site Voronoi diagram. Default: False
    incremental : bool, optional
        Allow adding new points incrementally. This takes up some additional
        resources.
    qhull_options : str, optional
        Additional options to pass to Qhull. See Qhull manual
        for details. (Default: "Qbb Qc Qz Qx" for ndim > 4 and
        "Qbb Qc Qz" otherwise. Incremental mode omits "Qz".)

    Attributes
    ----------
    points : ndarray of double, shape (npoints, ndim)
        Coordinates of input points.
    vertices : ndarray of double, shape (nvertices, ndim)
        Coordinates of the Voronoi vertices.
    ridge_points : ndarray of ints, shape ``(nridges, 2)``
        Indices of the points between which each Voronoi ridge lies.
    ridge_vertices : list of list of ints, shape ``(nridges, *)``
        Indices of the Voronoi vertices forming each Voronoi ridge.
    regions : list of list of ints, shape ``(nregions, *)``
        Indices of the Voronoi vertices forming each Voronoi region.
        -1 indicates vertex outside the Voronoi diagram.
    point_region : list of ints, shape (npoints)
        Index of the Voronoi region for each input point.
        If qhull option "Qc" was not specified, the list will contain -1
        for points that are not associated with a Voronoi region.

    Raises
    ------
    QhullError
        Raised when Qhull encounters an error condition, such as
        geometrical degeneracy when options to resolve are not enabled.
    ValueError
        Raised if an incompatible array is given as input.

    Notes
    -----
    The Voronoi diagram is computed using the 
    `Qhull library <http://www.qhull.org/>`__.

    Do not call the ``add_points`` method from a ``__del__``
    destructor.

    Examples
    --------
    Voronoi diagram for a set of point:

    >>> points = np.array([[0, 0], [0, 1], [0, 2], [1, 0], [1, 1], [1, 2],
    ...                    [2, 0], [2, 1], [2, 2]])
    >>> from scipy.spatial import Voronoi, voronoi_plot_2d
    >>> vor = Voronoi(points)

    Plot it:

    >>> import matplotlib.pyplot as plt
    >>> voronoi_plot_2d(vor)
    >>> plt.show()

    The Voronoi vertices:

    >>> vor.vertices
    array([[ 0.5,  0.5],
           [ 1.5,  0.5],
           [ 0.5,  1.5],
           [ 1.5,  1.5]])

    There is a single finite Voronoi region, and four finite Voronoi
    ridges:

    >>> vor.regions
    [[], [-1, 0], [-1, 1], [1, -1, 0], [3, -1, 2], [-1, 3], [-1, 2], [3, 2, 0, 1], [2, -1, 0], [3, -1, 1]]
    >>> vor.ridge_vertices
    [[-1, 0], [-1, 0], [-1, 1], [-1, 1], [0, 1], [-1, 3], [-1, 2], [2, 3], [-1, 3], [-1, 2], [0, 2], [1, 3]]

    The ridges are perpendicular between lines drawn between the following
    input points:

    >>> vor.ridge_points
    array([[0, 1],
           [0, 3],
           [6, 3],
           [6, 7],
           [3, 4],
           [5, 8],
           [5, 2],
           [5, 4],
           [8, 7],
           [2, 1],
           [4, 1],
           [4, 7]], dtype=int32)

    """
    def __init__(self, points, furthest_site=False, incremental=False,
                 qhull_options=None):
        if np.ma.isMaskedArray(points):
            raise ValueError('Input points cannot be a masked array')
        points = np.ascontiguousarray(points, dtype=np.double)

        if qhull_options is None:
            if not incremental:
                qhull_options = b"Qbb Qc Qz"
            else:
                qhull_options = b"Qc"
            if points.shape[1] >= 5:
                qhull_options += b" Qx"
        else:
            qhull_options = asbytes(qhull_options)

        # Run qhull
        qhull = _Qhull(b"v", points, qhull_options, furthest_site=furthest_site,
                       incremental=incremental)
        _QhullUser.__init__(self, qhull, incremental=incremental)

    def _update(self, qhull):
        self.vertices, self.ridge_points, self.ridge_vertices, \
                       self.regions, self.point_region = \
                       qhull.get_voronoi_diagram()

        self._ridge_dict = None

        _QhullUser._update(self, qhull)

    @property
    def ridge_dict(self):
        if self._ridge_dict is None:
            self._ridge_dict = dict(zip(map(tuple, self.ridge_points.tolist()),
                                        self.ridge_vertices))
        return self._ridge_dict
