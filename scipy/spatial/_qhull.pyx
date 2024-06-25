# cython: cpow=True

"""
Wrappers for Qhull triangulation, plus some additional N-D geometry utilities

.. versionadded:: 0.9

"""
#
# Copyright (C)  Pauli Virtanen, 2010.
#
# Distributed under the same BSD license as Scipy.
#


import numpy as np
cimport numpy as np
cimport cython
from cpython.pythread cimport (
    PyThread_type_lock, PyThread_allocate_lock, PyThread_free_lock,
    PyThread_acquire_lock, PyThread_release_lock)

from . cimport _qhull
from . cimport setlist
from libc cimport stdlib
from libc.math cimport NAN
from scipy._lib.messagestream cimport MessageStream
from libc.stdio cimport FILE

from scipy.linalg.cython_lapack cimport dgetrf, dgetrs, dgecon

np.import_array()


__all__ = ['Delaunay', 'ConvexHull', 'QhullError', 'Voronoi', 'HalfspaceIntersection', 'tsearch']

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
cdef extern from "qhull_src/src/user_r.h":
    cdef enum:
        qh_ORIENTclock

cdef extern from "qhull_src/src/qset_r.h":
    ctypedef union setelemT:
        void *p
        int i

    ctypedef struct setT:
        int maxsize
        setelemT e[1]

    int qh_setsize(qhT *, setT *set) nogil
    void qh_setappend(qhT *, setT **setp, void *elem) nogil

cdef extern from "qhull_src/src/libqhull_r.h":
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
        flagT good
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
        boolT hasAreaVolume
        int normal_size
        char *qhull_command
        facetT *facet_list
        facetT *facet_tail
        vertexT *vertex_list
        vertexT *vertex_tail
        int num_facets
        int num_visible
        int num_vertices
        int center_size
        unsigned int facet_id
        int hull_dim
        int num_points
        pointT *first_point
        pointT *input_points
        coordT* feasible_point
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

    extern int qh_PRINToff
    extern int qh_ALL

    void qh_init_A(qhT *, void *inp, void *out, void *err, int argc, char **argv) nogil
    void qh_init_B(qhT *, realT *points, int numpoints, int dim, boolT ismalloc) nogil
    void qh_checkflags(qhT *, char *, char *) nogil
    void qh_initflags(qhT *, char *) nogil
    void qh_option(qhT *, char *, char*, char* ) nogil
    void qh_freeqhull(qhT *, boolT) nogil
    void qh_memfreeshort(qhT *, int *curlong, int *totlong) nogil
    void qh_qhull(qhT *) nogil
    void qh_check_output(qhT *) nogil
    void qh_produce_output(qhT *) nogil
    void qh_triangulate(qhT *) nogil
    void qh_checkpolygon(qhT *) nogil
    void qh_findgood_all(qhT *, facetT *facetlist) nogil
    void qh_appendprint(qhT *, int format) nogil
    setT *qh_pointvertex(qhT *) nogil
    realT *qh_readpoints(qhT *, int* num, int *dim, boolT* ismalloc) nogil
    void qh_zero(qhT *, void *errfile) nogil
    int qh_pointid(qhT *, pointT *point) nogil
    vertexT *qh_nearvertex(qhT *, facetT *facet, pointT *point, double *dist) nogil
    boolT qh_addpoint(qhT *, pointT *furthest, facetT *facet, boolT checkdist) nogil
    facetT *qh_findbestfacet(qhT *, pointT *point, boolT bestoutside,
                             realT *bestdist, boolT *isoutside) nogil
    void qh_setdelaunay(qhT *, int dim, int count, pointT *points) nogil
    coordT* qh_sethalfspace_all(qhT *, int dim, int count, coordT* halfspaces, pointT *feasible)

cdef extern from "qhull_misc.h":
    ctypedef int CBLAS_INT   # actual type defined in the header file
    void qhull_misc_lib_check()
    int qh_new_qhull_scipy(qhT *, int dim, int numpoints, realT *points,
                           boolT ismalloc, char* qhull_cmd, void *outfile,
                           void *errfile, coordT* feaspoint) nogil

cdef extern from "qhull_src/src/io_r.h":
    ctypedef enum qh_RIDGE:
        qh_RIDGEall
        qh_RIDGEinner
        qh_RIDGEouter

    ctypedef void printvridgeT(qhT *, void *fp, vertexT *vertex, vertexT *vertexA,
                               setT *centers, boolT unbounded)
    int qh_eachvoronoi_all(qhT *, FILE *fp, void* printvridge,
                           boolT isUpper, qh_RIDGE innerouter,
                           boolT inorder) nogil

    void qh_order_vertexneighbors(qhT *, vertexT *vertex) nogil
    int qh_compare_facetvisit(const void *p1, const void *p2) nogil

cdef extern from "qhull_src/src/geom_r.h":
    pointT *qh_facetcenter(qhT *, setT *vertices) nogil
    double qh_getarea(qhT *, facetT *facetlist) nogil

cdef extern from "qhull_src/src/poly_r.h":
    void qh_check_maxout(qhT *) nogil

cdef extern from "qhull_src/src/mem_r.h":
    void qh_memfree(qhT *, void *object, int insize)

from libc.stdlib cimport qsort


#------------------------------------------------------------------------------
# Qhull wrapper
#------------------------------------------------------------------------------

# Check Qhull library compatibility at import time
qhull_misc_lib_check()


class QhullError(RuntimeError):
    pass


@cython.final
cdef class _Qhull:
    # Note that the qhT struct is allocated separately --- otherwise
    # it may end up allocated in a way not compatible with the CRT
    # (on Windows)
    cdef qhT *_qh

    cdef list _point_arrays
    cdef list _dual_point_arrays
    cdef MessageStream _messages

    cdef public bytes options
    cdef public bytes mode_option
    cdef public object furthest_site

    cdef readonly int ndim
    cdef int numpoints, _is_delaunay, _is_halfspaces
    cdef np.ndarray _ridge_points

    cdef list _ridge_vertices
    cdef object _ridge_error
    cdef int _nridges

    cdef np.ndarray _ridge_equations
    cdef PyThread_type_lock _lock

    @cython.final
    def __init__(self,
                 bytes mode_option,
                 np.ndarray[np.double_t, ndim=2] points,
                 bytes options=None,
                 bytes required_options=None,
                 furthest_site=False,
                 incremental=False,
                 np.ndarray[np.double_t, ndim=1] interior_point=None):
        cdef int exitcode

        self._lock = PyThread_allocate_lock()
        if self._lock == NULL:
            raise MemoryError("thread lock allocation failed")

        self._qh = NULL
        self._messages = MessageStream()

        points = np.ascontiguousarray(points, dtype=np.double)

        self.ndim = points.shape[1]

        self.numpoints = points.shape[0]

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
            bad_opts = [bad_opt for bad_opt in incremental_bad_ops
                        if bad_opt in options]
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

        if mode_option.startswith(b"H"):
            self._is_halfspaces = 1
        else:
            self._is_halfspaces = 0

        self._point_arrays = [points]
        self._dual_point_arrays = []
        self.options = b" ".join(option_set)
        self.mode_option = mode_option
        self.furthest_site = furthest_site

        options = b"qhull "  + mode_option +  b" " + self.options

        options_c = <char*>options

        self._messages.clear()

        cdef coordT* coord
        with nogil:
            self._qh = <qhT*>stdlib.malloc(sizeof(qhT))
            if self._qh == NULL:
                with gil:
                    raise MemoryError("memory allocation failed")
            qh_zero(self._qh, self._messages.handle)
            if interior_point is not None:
                coord = <coordT*>interior_point.data
            else:
                coord = NULL
            exitcode = qh_new_qhull_scipy(self._qh, self.ndim, self.numpoints,
                                          <realT*>points.data, 0,
                                          options_c, NULL, self._messages.handle, coord)

        if exitcode != 0:
            msg = self._messages.get()
            self.close()
            raise QhullError(msg)

    cdef void acquire_lock(self):
        if not PyThread_acquire_lock(self._lock, 0):
            PyThread_acquire_lock(self._lock, 1)

    cdef void release_lock(self):
        PyThread_release_lock(self._lock)

    def check_active(self):
        if self._qh == NULL:
            raise RuntimeError("Qhull instance is closed")

    @cython.final
    def __dealloc__(self):
        cdef int curlong, totlong

        self.acquire_lock()
        try:
            if self._qh != NULL:
                qh_freeqhull(self._qh, qh_ALL)
                qh_memfreeshort(self._qh, &curlong, &totlong)
                stdlib.free(self._qh)
                self._qh = NULL

                if curlong != 0 or totlong != 0:
                    raise QhullError(
                        "qhull: did not free %d bytes (%d pieces)" %
                        (totlong, curlong))

            if self._messages is not None:
                self._messages.close()
        finally:
            self.release_lock()

        PyThread_free_lock(self._lock)

    @cython.final
    def close(self):
        """
        Uninitialize this instance
        """
        # Note: this is direct copypaste from __dealloc__(), keep it
        # in sync with that.  The code must be written directly in
        # __dealloc__, because otherwise the generated C code tries to
        # call PyObject_GetAttrStr(self, "close") which on Pypy
        # crashes.

        cdef int curlong, totlong

        self.acquire_lock()
        try:
            if self._qh != NULL:
                qh_freeqhull(self._qh, qh_ALL)
                qh_memfreeshort(self._qh, &curlong, &totlong)
                stdlib.free(self._qh)
                self._qh = NULL

                if curlong != 0 or totlong != 0:
                    raise QhullError(
                        "qhull: did not free %d bytes (%d pieces)" %
                        (totlong, curlong))

            if self._messages is not None:
                self._messages.close()
        finally:
            self.release_lock()

    @cython.final
    def get_points(self):
        if len(self._point_arrays) == 1:
            return self._point_arrays[0]
        else:
            return np.concatenate(
                [x[:,:self.ndim] for x in self._point_arrays],
                axis=0)

    @cython.final
    def add_points(self, points, interior_point=None):
        cdef int j
        cdef realT *p
        cdef facetT *facet
        cdef double bestdist
        cdef boolT isoutside
        cdef np.ndarray arr

        self.acquire_lock()

        try:
            self.check_active()

            points = np.asarray(points)
            if points.ndim!=2 or points.shape[1] != self._point_arrays[0].shape[1]:
                raise ValueError("invalid size for new points array")
            if points.size == 0:
                return

            if self._is_delaunay:
                arr = np.empty((points.shape[0], self.ndim+1), dtype=np.double)
                arr[:,:-1] = points
            elif self._is_halfspaces:
                #Store the halfspaces in _points and the dual points in _dual_points later
                self._point_arrays.append(np.array(points, copy=True))
                dists = points[:, :-1].dot(interior_point)+points[:, -1]
                arr = np.array(-points[:, :-1]/dists, dtype=np.double, order="C", copy=True)
            else:
                arr = np.array(points, dtype=np.double, order="C", copy=True)

            self._messages.clear()

            try:
                # nonlocal error handling
                exitcode = setjmp(self._qh[0].errexit)
                if exitcode != 0:
                    raise QhullError(self._messages.get())
                self._qh[0].NOerrexit = 0

                # add points to triangulation
                if self._is_delaunay:
                    # lift to paraboloid
                    qh_setdelaunay(self._qh, arr.shape[1], arr.shape[0], <realT*>arr.data)

                p = <realT*>arr.data

                for j in range(arr.shape[0]):
                    facet = qh_findbestfacet(self._qh, p, 0, &bestdist, &isoutside)
                    if isoutside:
                        if not qh_addpoint(self._qh, p, facet, 0):
                            break
                    else:
                        # append the point to the "other points" list, to
                        # maintain the point IDs
                        qh_setappend(self._qh, &self._qh[0].other_points, p)

                    p += arr.shape[1]

                qh_check_maxout(self._qh)
                self._qh[0].hasTriangulation = 0

                if self._is_halfspaces:
                    self._dual_point_arrays.append(arr)
                else:
                    self._point_arrays.append(arr)
                self.numpoints += arr.shape[0]

                # update facet visibility
                with nogil:
                    qh_findgood_all(self._qh, self._qh[0].facet_list)
            finally:
                self._qh[0].NOerrexit = 1
        finally:
            self.release_lock()

    @cython.final
    def get_paraboloid_shift_scale(self):
        cdef double paraboloid_scale
        cdef double paraboloid_shift

        self.acquire_lock()
        try:
            self.check_active()

            if self._qh[0].SCALElast:
                paraboloid_scale = self._qh[0].last_newhigh / (
                    self._qh[0].last_high - self._qh[0].last_low)
                paraboloid_shift = - self._qh[0].last_low * paraboloid_scale
            else:
                paraboloid_scale = 1.0
                paraboloid_shift = 0.0

            return paraboloid_scale, paraboloid_shift
        finally:
            self.release_lock()

    @cython.final
    def volume_area(self):
        cdef double volume
        cdef double area

        self.acquire_lock()
        try:
            self.check_active()

            self._qh.hasAreaVolume = 0
            with nogil:
                qh_getarea(self._qh, self._qh[0].facet_list)

            volume = self._qh[0].totvol
            area = self._qh[0].totarea

            return volume, area
        finally:
            self.release_lock()

    @cython.final
    def triangulate(self):
        self.acquire_lock()
        try:
            self.check_active()

            with nogil:
                qh_triangulate(self._qh) # get rid of non-simplical facets
        finally:
            self.release_lock()

    @cython.final
    @cython.boundscheck(False)
    @cython.cdivision(True)
    def get_simplex_facet_array(self):
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
        cdef int i, j, ipoint, ncoplanar
        cdef object tmp
        cdef np.ndarray[np.npy_int, ndim=2] facets
        cdef np.ndarray[np.npy_int, ndim=1] good
        cdef np.ndarray[np.npy_int, ndim=2] neighbors
        cdef np.ndarray[np.npy_int, ndim=2] coplanar
        cdef np.ndarray[np.double_t, ndim=2] equations
        cdef np.ndarray[np.npy_int, ndim=1] id_map
        cdef double dist
        cdef int facet_ndim
        cdef unsigned int lower_bound
        cdef unsigned int swapped_index

        self.acquire_lock()
        try:
            self.check_active()

            facet_ndim = self.ndim

            if self._is_halfspaces:
                facet_ndim = self.ndim - 1

            if self._is_delaunay:
                facet_ndim += 1

            id_map = np.empty(self._qh[0].facet_id, dtype=np.intc)

            # Compute facet indices
            with nogil:
                for i in range(self._qh[0].facet_id):
                    id_map[i] = -1

                facet = self._qh[0].facet_list
                j = 0
                while facet and facet.next:
                    if not self._is_delaunay or facet.upperdelaunay == self._qh[0].UPPERdelaunay:
                        if not facet.simplicial and ( \
                            qh_setsize(self._qh, facet.vertices) != facet_ndim or \
                            qh_setsize(self._qh, facet.neighbors) != facet_ndim):
                            with gil:
                                raise QhullError(
                                    "non-simplical facet encountered: %r vertices"
                                    % (qh_setsize(self._qh, facet.vertices),))

                        id_map[facet.id] = j
                        j += 1

                    facet = facet.next

            # Allocate output
            facets = np.zeros((j, facet_ndim), dtype=np.intc)
            good = np.zeros(j, dtype=np.intc)
            neighbors = np.zeros((j, facet_ndim), dtype=np.intc)
            equations = np.zeros((j, facet_ndim+1), dtype=np.double)

            ncoplanar = 0
            coplanar = np.zeros((10, 3), dtype=np.intc)
            coplanar_shape = coplanar.shape

            # Retrieve facet information
            with nogil:
                facet = self._qh[0].facet_list
                j = 0
                while facet and facet.next:
                    if self._is_delaunay and facet.upperdelaunay != self._qh[0].UPPERdelaunay:
                        facet = facet.next
                        continue

                    # Use a lower bound so that the tight loop in high dimensions
                    # is not affected by the conditional below
                    lower_bound = 0
                    if (self._is_delaunay and
                        facet.toporient == qh_ORIENTclock and facet_ndim == 3):
                        # Swap the first and second indices to maintain a
                        # counter-clockwise orientation.
                        for i in range(2):
                            # Save the vertex info
                            swapped_index = 1 ^ i
                            vertex = <vertexT*>facet.vertices.e[i].p
                            ipoint = qh_pointid(self._qh, vertex.point)
                            facets[j, swapped_index] = ipoint

                            # Save the neighbor info
                            neighbor = <facetT*>facet.neighbors.e[i].p
                            neighbors[j, swapped_index] = id_map[neighbor.id]

                        lower_bound = 2

                    for i in range(lower_bound, facet_ndim):
                        # Save the vertex info
                        vertex = <vertexT*>facet.vertices.e[i].p
                        ipoint = qh_pointid(self._qh, vertex.point)
                        facets[j, i] = ipoint

                        # Save the neighbor info
                        neighbor = <facetT*>facet.neighbors.e[i].p
                        neighbors[j, i] = id_map[neighbor.id]

                    # Save simplex equation info
                    for i in range(facet_ndim):
                        equations[j, i] = facet.normal[i]
                    equations[j, facet_ndim] = facet.offset

                    # Save coplanar info
                    if facet.coplanarset:
                        for i in range(qh_setsize(self._qh, facet.coplanarset)):
                            point = <pointT*>facet.coplanarset.e[i].p
                            vertex = qh_nearvertex(self._qh, facet, point, &dist)

                            if ncoplanar >= coplanar_shape[0]:
                                with gil:
                                    tmp = coplanar
                                    coplanar = None
                                    # The array is always safe to resize
                                    tmp.resize(2 * ncoplanar + 1, 3, refcheck=False)
                                    coplanar = tmp

                            coplanar[ncoplanar, 0] = qh_pointid(self._qh, point)
                            coplanar[ncoplanar, 1] = id_map[facet.id]
                            coplanar[ncoplanar, 2] = qh_pointid(self._qh, vertex.point)
                            ncoplanar += 1

                    # Save good info
                    good[j] = facet.good

                    j += 1
                    facet = facet.next

            return facets, neighbors, equations, coplanar[:ncoplanar], good
        finally:
            self.release_lock()

    @cython.final
    @cython.boundscheck(False)
    @cython.cdivision(True)
    def get_hull_points(self):
        """Returns all points currently contained in Qhull.
        It is equivalent to retrieving the input in most cases, except in
        halfspace mode, where the points are in fact the points of the dual
        hull.

        Returns
        -------
        points: array of double, shape (nrpoints, ndim)
            The array of points contained in Qhull.

        """
        cdef pointT *point
        cdef int i, j, numpoints, point_ndim
        cdef np.ndarray[np.npy_double, ndim=2] points

        self.acquire_lock()

        try:
            self.check_active()

            point_ndim = self.ndim

            if self._is_halfspaces:
                point_ndim -= 1

            if self._is_delaunay:
                point_ndim += 1

            numpoints = self._qh.num_points
            points = np.empty((numpoints, point_ndim))

            with nogil:
                point = self._qh.first_point
                for i in range(numpoints):
                    for j in range(point_ndim):
                        points[i,j] = point[j]
                    point += self._qh.hull_dim

            return points
        finally:
            self.release_lock()

    @cython.final
    @cython.boundscheck(False)
    @cython.cdivision(True)
    def get_hull_facets(self):
        """Returns the facets contained in the current Qhull.
        This function does not assume that the hull is simplicial,
        meaning that facets will have different number of vertices.
        It is thus less efficient but more general than get_simplex_facet_array.

        Returns
        -------
        facets: list of lists of ints
            The indices of the vertices forming each facet.
        """
        cdef facetT *facet
        cdef vertexT* vertex
        cdef int i, j, numfacets, facet_ndim
        cdef np.ndarray[np.double_t, ndim=2] equations
        cdef list facets, facetsi

        self.acquire_lock()
        try:
            self.check_active()

            facet_ndim = self.ndim

            if self._is_halfspaces:
                facet_ndim -= 1

            if self._is_delaunay:
                facet_ndim += 1

            numfacets = self._qh.num_facets - self._qh.num_visible

            facet = self._qh.facet_list
            equations = np.empty((numfacets, facet_ndim+1))

            facets = []

            i = 0
            while facet and facet.next:
                facetsi = []
                j = 0
                for j in range(facet_ndim):
                    equations[i, j] = facet.normal[j]
                equations[i, facet_ndim] = facet.offset

                j = 0
                vertex = <vertexT*>facet.vertices.e[0].p
                while vertex:
                    # Save the vertex info
                    ipoint = qh_pointid(self._qh, vertex.point)
                    facetsi.append(ipoint)
                    j += 1
                    vertex = <vertexT*>facet.vertices.e[j].p

                i += 1
                facets.append(facetsi)
                facet = facet.next

            return facets, equations
        finally:
            self.release_lock()

    @cython.final
    @cython.boundscheck(False)
    @cython.cdivision(True)
    def get_voronoi_diagram(_Qhull self):
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
        cdef pointT *point
        cdef pointT *center
        cdef double dist
        cdef int inf_seen

        cdef list regions
        cdef list cur_region

        self.acquire_lock()
        try:
            self.check_active()

            # -- Grab Voronoi ridges
            self._nridges = 0
            self._ridge_error = None
            self._ridge_points = np.empty((10, 2), np.intc)
            self._ridge_vertices = []

            qh_eachvoronoi_all(self._qh, <FILE*>self, &_visit_voronoi, self._qh[0].UPPERdelaunay,
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

            vertex = self._qh[0].vertex_list
            while vertex and vertex.next:
                qh_order_vertexneighbors_nd(self._qh, self.ndim+1, vertex)

                i = qh_pointid(self._qh, vertex.point)
                if i < self.numpoints:
                    # Qz results to one extra point
                    point_region[i] = len(regions)

                inf_seen = 0
                cur_region = []
                for k in range(qh_setsize(self._qh, vertex.neighbors)):
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

            facet = self._qh[0].facet_list
            while facet and facet.next:
                if facet.visitid > 0:
                    # finite Voronoi vertex

                    center = qh_facetcenter(self._qh, facet.vertices)

                    nvoronoi_vertices = max(facet.visitid, nvoronoi_vertices)
                    if nvoronoi_vertices >= voronoi_vertices.shape[0]:
                        tmp = voronoi_vertices
                        voronoi_vertices = None
                        # Array is safe to resize
                        tmp.resize(2*nvoronoi_vertices + 1, self.ndim, refcheck=False)
                        voronoi_vertices = tmp

                    for k in range(self.ndim):
                        voronoi_vertices[facet.visitid-1, k] = center[k]

                    qh_memfree(self._qh, center, self._qh[0].center_size)

                    if facet.coplanarset:
                        for k in range(qh_setsize(self._qh, facet.coplanarset)):
                            point = <pointT*>facet.coplanarset.e[k].p
                            vertex = qh_nearvertex(self._qh, facet, point, &dist)

                            i = qh_pointid(self._qh, point)
                            j = qh_pointid(self._qh, vertex.point)

                            if i < self.numpoints:
                                # Qz can result to one extra point
                                point_region[i] = point_region[j]

                facet = facet.next

            voronoi_vertices = voronoi_vertices[:nvoronoi_vertices]

            return voronoi_vertices, self._ridge_points, self._ridge_vertices, \
                regions, point_region
        finally:
            self.release_lock()

    @cython.final
    @cython.boundscheck(False)
    @cython.cdivision(True)
    def get_extremes_2d(_Qhull self):
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

        self.acquire_lock()
        try:
            self.check_active()

            if self._is_delaunay:
                raise ValueError("Cannot compute for Delaunay")

            nextremes = 0
            extremes_arr = np.zeros(100, dtype=np.intc)
            extremes = extremes_arr

            self._qh[0].visit_id += 1
            self._qh[0].vertex_visit += 1

            facet = self._qh[0].facet_list
            startfacet = facet
            while facet:
                if facet.visitid == self._qh[0].visit_id:
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
                    # Array is safe to resize
                    extremes_arr.resize(2*extremes_arr.shape[0]+1, refcheck=False)
                    extremes = extremes_arr

                if vertexA.visitid != self._qh[0].vertex_visit:
                    vertexA.visitid = self._qh[0].vertex_visit
                    extremes[nextremes] = qh_pointid(self._qh, vertexA.point)
                    nextremes += 1

                if vertexB.visitid != self._qh[0].vertex_visit:
                    vertexB.visitid = self._qh[0].vertex_visit
                    extremes[nextremes] = qh_pointid(self._qh, vertexB.point)
                    nextremes += 1

                facet.visitid = self._qh[0].visit_id
                facet = nextfacet

                if facet == startfacet:
                    break

            extremes = None
            # This array is always safe to resize
            extremes_arr.resize(nextremes, refcheck=False)
            return extremes_arr
        finally:
            self.release_lock()


cdef void _visit_voronoi(qhT *_qh, FILE *ptr, vertexT *vertex, vertexT *vertexA,
                         setT *centers, boolT unbounded) noexcept:
    cdef _Qhull qh = <_Qhull>ptr
    cdef int point_1, point_2, ix
    cdef list cur_vertices

    if qh._ridge_error is not None:
        return

    if qh._nridges >= qh._ridge_points.shape[0]:
        try:
            # The array is guaranteed to be safe to resize
            qh._ridge_points.resize(2*qh._nridges + 1, 2, refcheck=False)
        except Exception, e:
            qh._ridge_error = e
            return

    # Record which points the ridge is between
    point_1 = qh_pointid(_qh, vertex.point)
    point_2 = qh_pointid(_qh, vertexA.point)

    p = <int*>qh._ridge_points.data
    p[2*qh._nridges + 0] = point_1
    p[2*qh._nridges + 1] = point_2

    # Record which voronoi vertices constitute the ridge
    cur_vertices = []
    for i in range(qh_setsize(_qh, centers)):
        ix = (<facetT*>centers.e[i].p).visitid - 1
        cur_vertices.append(ix)
    qh._ridge_vertices.append(cur_vertices)

    qh._nridges += 1

    return


cdef void qh_order_vertexneighbors_nd(qhT *qh, int nd, vertexT *vertex) noexcept:
    if nd == 3:
        qh_order_vertexneighbors(qh, vertex)
    elif nd >= 4:
        qsort(<facetT**>&vertex.neighbors.e[0].p, qh_setsize(qh, vertex.neighbors),
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
    cdef int i, j
    cdef CBLAS_INT n, nrhs, lda, ldb
    cdef CBLAS_INT info = 0
    cdef CBLAS_INT ipiv[NPY_MAXDIMS+1]
    cdef int ndim, nsimplex
    cdef double anorm
    cdef double rcond = 0.0
    cdef double rcond_limit

    cdef double work[4*NPY_MAXDIMS]
    cdef CBLAS_INT iwork[NPY_MAXDIMS]

    ndim = points.shape[1]
    nsimplex = simplices.shape[0]

    T = np.empty((ndim, ndim), dtype=np.double)
    Tinvs = np.zeros((nsimplex, ndim+1, ndim), dtype=np.double)

    # Maximum inverse condition number to allow: we want at least three
    # of the digits be significant, to be safe
    rcond_limit = 1000*eps

    with nogil:
        for isimplex in range(nsimplex):
            for i in range(ndim):
                Tinvs[isimplex,ndim,i] = points[simplices[isimplex,ndim],i]
                for j in range(ndim):
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
            dgetrf(&n, &n, <double*>T.data, &lda, ipiv, &info)

            # Check condition number
            if info == 0:
                dgecon("1", &n, <double*>T.data, &lda, &anorm, &rcond,
                       work, iwork, &info)

                if rcond < rcond_limit:
                    # The transform seems singular
                    info = 1

            # Compute transform
            if info == 0:
                dgetrs("N", &n, &nrhs, <double*>T.data, &lda, ipiv,
                          (<double*>Tinvs.data) + ndim*(ndim+1)*isimplex,
                          &ldb, &info)

            # Deal with degenerate simplices
            if info != 0:
                for i in range(ndim+1):
                    for j in range(ndim):
                        Tinvs[isimplex,i,j] = NAN

    return Tinvs

@cython.boundscheck(False)
cdef double _matrix_norm1(int n, double *a) noexcept nogil:
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
                             const double *x, double *c, double eps) noexcept nogil:
    """
    Check whether point is inside a simplex, using barycentric
    coordinates.  `c` will be filled with barycentric coordinates, if
    the point happens to be inside.

    """
    cdef int i, j
    c[ndim] = 1.0
    for i in range(ndim):
        c[i] = 0
        for j in range(ndim):
            c[i] += transform[ndim*i + j] * (x[j] - transform[ndim*ndim + j])
        c[ndim] -= c[i]

        if not (-eps <= c[i] <= 1 + eps):
            return 0
    if not (-eps <= c[ndim] <= 1 + eps):
        return 0
    return 1

cdef void _barycentric_coordinate_single(int ndim, double *transform,
                                         const double *x, double *c, int i) noexcept nogil:
    """
    Compute a single barycentric coordinate.

    Before the ndim+1'th coordinate can be computed, the other must have
    been computed earlier.

    """
    cdef int j

    if i == ndim:
        c[ndim] = 1.0
        for j in range(ndim):
            c[ndim] -= c[j]
    else:
        c[i] = 0
        for j in range(ndim):
            c[i] += transform[ndim*i + j] * (x[j] - transform[ndim*ndim + j])

cdef void _barycentric_coordinates(int ndim, double *transform,
                                   const double *x, double *c) noexcept nogil:
    """
    Compute barycentric coordinates.

    """
    cdef int i, j
    c[ndim] = 1.0
    for i in range(ndim):
        c[i] = 0
        for j in range(ndim):
            c[i] += transform[ndim*i + j] * (x[j] - transform[ndim*ndim + j])
        c[ndim] -= c[i]


#------------------------------------------------------------------------------
# N-D geometry
#------------------------------------------------------------------------------

cdef void _lift_point(DelaunayInfo_t *d, const double *x, double *z) noexcept nogil:
    cdef int i
    z[d.ndim] = 0
    for i in range(d.ndim):
        z[i] = x[i]
        z[d.ndim] += x[i]**2
    z[d.ndim] *= d.paraboloid_scale
    z[d.ndim] += d.paraboloid_shift

cdef double _distplane(DelaunayInfo_t *d, int isimplex, double *point) noexcept nogil:
    """
    qh_distplane
    """
    cdef double dist
    cdef int k
    dist = d.equations[isimplex*(d.ndim+2) + d.ndim+1]
    for k in range(d.ndim+1):
        dist += d.equations[isimplex*(d.ndim+2) + k] * point[k]
    return dist


#------------------------------------------------------------------------------
# Finding simplices
#------------------------------------------------------------------------------

cdef int _is_point_fully_outside(DelaunayInfo_t *d, const double *x,
                                 double eps) noexcept nogil:
    """
    Is the point outside the bounding box of the triangulation?

    """

    cdef int i
    for i in range(d.ndim):
        if x[i] < d.min_bound[i] - eps or x[i] > d.max_bound[i] + eps:
            return 1
    return 0

cdef int _find_simplex_bruteforce(DelaunayInfo_t *d, double *c,
                                  const double *x, double eps,
                                  double eps_broad) noexcept nogil:
    """
    Find simplex containing point `x` by going through all simplices.

    """
    cdef int inside, isimplex
    cdef int k, m, ineighbor
    cdef double *transform

    if _is_point_fully_outside(d, x, eps):
        return -1

    for isimplex in range(d.nsimplex):
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

            for k in range(d.ndim+1):
                ineighbor = d.neighbors[(d.ndim+1)*isimplex + k]
                if ineighbor == -1:
                    continue

                transform = d.transform + ineighbor*d.ndim*(d.ndim+1)
                if transform[0] != transform[0]:
                    # another bad simplex
                    continue

                _barycentric_coordinates(d.ndim, transform, x, c)

                # Check that the point lies (almost) inside the
                # neighbor simplex
                inside = 1
                for m in range(d.ndim+1):
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
                                const double *x, int *start, double eps,
                                double eps_broad) noexcept nogil:
    """
    Find simplex containing point `x` via a directed walk in the tessellation.

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

    4) In a regular tessellation, hopping to any such direction is OK.

       Also, if one of the negative-coordinate neighbors happens to be -1,
       then the target point is outside the tessellation (because the
       tessellation is convex!).

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
        for k in range(ndim+1):
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
                       const double *x, int *start, double eps,
                       double eps_broad) noexcept nogil:
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

    In this algorithm, we walk around the tessellation trying to locate
    a positive-distant facet. After finding one, we fall back to a
    directed search.

    """
    cdef int isimplex, k, ineigh
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

    # Walk the tessellation searching for a facet with a positive planar distance
    best_dist = _distplane(d, isimplex, z)
    changed = 1
    while changed:
        if best_dist > 0:
            break
        changed = 0
        for k in range(ndim+1):
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

class _QhullUser:
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
        self._points = qhull.get_points()
        self.ndim = self._points.shape[1]
        self.npoints = self._points.shape[0]
        self.min_bound = self._points.min(axis=0)
        self.max_bound = self._points.max(axis=0)

    def _add_points(self, points, restart=False, interior_point=None):
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
            points = np.concatenate([self._points, points], axis=0)
            qhull = _Qhull(self._qhull.mode_option, points,
                           options=self._qhull.options,
                           furthest_site=self._qhull.furthest_site,
                           incremental=True, interior_point=interior_point)
            try:
                self._update(qhull)
                self._qhull = qhull
            finally:
                if qhull is not self._qhull:
                    qhull.close()
            return

        self._qhull.add_points(points, interior_point)
        self._update(self._qhull)

class Delaunay(_QhullUser):
    """
    Delaunay(points, furthest_site=False, incremental=False, qhull_options=None)

    Delaunay tessellation in N dimensions.

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
        Default:"Qbb Qc Qz Qx Q12" for ndim > 4 and "Qbb Qc Qz Q12" otherwise.
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
        for points that are not vertices of the tessellation.
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
    vertex_neighbor_vertices : tuple of two ndarrays of int; (indptr, indices)
        Neighboring vertices of vertices. The indices of neighboring
        vertices of vertex `k` are ``indices[indptr[k]:indptr[k+1]]``.
    furthest_site
        True if this was a furthest site triangulation and False if not.

        .. versionadded:: 1.4.0

    Raises
    ------
    QhullError
        Raised when Qhull encounters an error condition, such as
        geometrical degeneracy when options to resolve are not enabled.
    ValueError
        Raised if an incompatible array is given as input.

    Notes
    -----
    The tessellation is computed using the Qhull library
    `Qhull library <http://www.qhull.org/>`__.

    .. note::

       Unless you pass in the Qhull option "QJ", Qhull does not
       guarantee that each input point appears as a vertex in the
       Delaunay triangulation. Omitted points are listed in the
       `coplanar` attribute.

    Examples
    --------
    Triangulation of a set of points:

    >>> import numpy as np
    >>> points = np.array([[0, 0], [0, 1.1], [1, 0], [1, 1]])
    >>> from scipy.spatial import Delaunay
    >>> tri = Delaunay(points)

    We can plot it:

    >>> import matplotlib.pyplot as plt
    >>> plt.triplot(points[:,0], points[:,1], tri.simplices)
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

    >>> p = np.array([(0.1, 0.2), (1.5, 0.5), (0.5, 1.05)])
    >>> tri.find_simplex(p)
    array([ 1, -1, 1], dtype=int32)

    The returned integers in the array are the indices of the simplex the
    corresponding point is in. If -1 is returned, the point is in no simplex.
    Be aware that the shortcut in the following example only works correctly
    for valid points as invalid points result in -1 which is itself a valid
    index for the last simplex in the list.

    >>> p_valids = np.array([(0.1, 0.2), (0.5, 1.05)])
    >>> tri.simplices[tri.find_simplex(p_valids)]
    array([[3, 1, 0],                 # may vary
           [3, 1, 0]], dtype=int32)

    We can also compute barycentric coordinates in triangle 1 for
    these points:

    >>> b = tri.transform[1,:2].dot(np.transpose(p - tri.transform[1,2]))
    >>> np.c_[np.transpose(b), 1 - b.sum(axis=0)]
    array([[ 0.1       ,  0.09090909,  0.80909091],
           [ 1.5       , -0.90909091,  0.40909091],
           [ 0.5       ,  0.5       ,  0.        ]])

    The coordinates for the first point are all positive, meaning it
    is indeed inside the triangle. The third point is on an edge,
    hence its null third coordinate.

    """

    def __init__(self, points, furthest_site=False, incremental=False,
                 qhull_options=None):
        if np.ma.isMaskedArray(points):
            raise ValueError('Input points cannot be a masked array')
        points = np.ascontiguousarray(points, dtype=np.double)
        if points.ndim != 2:
            raise ValueError("Input points array must have 2 dimensions.")

        if qhull_options is None:
            if not incremental:
                qhull_options = b"Qbb Qc Qz Q12"
            else:
                qhull_options = b"Qc"
            if points.shape[1] >= 5:
                qhull_options += b" Qx"
        else:
            qhull_options = qhull_options.encode('latin1')

        # Run qhull
        qhull = _Qhull(b"d", points, qhull_options, required_options=b"Qt",
                       furthest_site=furthest_site, incremental=incremental)
        _QhullUser.__init__(self, qhull, incremental=incremental)

        self.furthest_site = furthest_site

    def _update(self, qhull):
        qhull.triangulate()

        self.paraboloid_scale, self.paraboloid_shift = \
                               qhull.get_paraboloid_shift_scale()

        self.simplices, self.neighbors, self.equations, self.coplanar, self.good = \
                       qhull.get_simplex_facet_array()

        self.nsimplex = self.simplices.shape[0]
        self._transform = None
        self._vertex_to_simplex = None
        self._vertex_neighbor_vertices = None

        _QhullUser._update(self, qhull)

    def add_points(self, points, restart=False):
        self._add_points(points, restart)

    @property
    def points(self):
        return self._points

    @property
    def transform(self):
        """
        Affine transform from ``x`` to the barycentric coordinates ``c``.

        :type: *ndarray of double, shape (nsimplex, ndim+1, ndim)*

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

        :type: *ndarray of int, shape (npoints,)*
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

            nsimplex = self.nsimplex
            ndim = self.ndim

            with nogil:
                for isimplex in range(nsimplex):
                    for k in range(ndim+1):
                        ivertex = simplices[isimplex, k]
                        if arr[ivertex] == -1:
                            arr[ivertex] = isimplex

        return self._vertex_to_simplex

    @property
    @cython.boundscheck(False)
    def vertex_neighbor_vertices(self):
        """
        Neighboring vertices of vertices.

        Tuple of two ndarrays of int: (indptr, indices). The indices of
        neighboring vertices of vertex `k` are
        ``indices[indptr[k]:indptr[k+1]]``.

        """
        cdef int i, j, k
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
                    for i in range(nsimplex):
                        for j in range(ndim+1):
                            for k in range(ndim+1):
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

        :type: *ndarray of int, shape (nfaces, ndim)*

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
        for isimplex in range(nsimplex):
            for k in range(ndim+1):
                if neighbors[isimplex,k] == -1:
                    for j in range(ndim+1):
                        if j < k:
                            arr[m,j] = simplices[isimplex,j]
                        elif j > k:
                            arr[m,j-1] = simplices[isimplex,j]
                    m += 1

                    if m >= msize:
                        arr = None
                        msize = 2*msize + 1
                        # Array is safe to resize
                        out.resize(msize, ndim, refcheck=False)
                        arr = out

        arr = None
        # Array is safe to resize
        out.resize(m, ndim, refcheck=False)
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
        x_shape = x.shape

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
                for k in range(x_shape[0]):
                    isimplex = _find_simplex_bruteforce(
                        &info, c,
                        <double*>x.data + info.ndim*k,
                        eps, eps_broad)
                    out_[k] = isimplex
        else:
            with nogil:
                for k in range(x_shape[0]):
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
        cdef int i, j

        if xi.shape[-1] != self.ndim:
            raise ValueError("xi has different dimensionality than "
                             "triangulation")

        xi_shape = xi.shape
        xi = xi.reshape(-1, xi.shape[-1])
        x = np.ascontiguousarray(xi.astype(np.double))
        x_shape = x.shape

        _get_delaunay_info(&info, self, 0, 0, 0)

        out = np.zeros((x.shape[0], info.nsimplex), dtype=np.double)
        out_ = out

        with nogil:
            for i in range(x_shape[0]):
                for j in range(info.nsimplex):
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

    Parameters
    ----------
    tri : DelaunayInfo
        Delaunay triangulation
    xi : ndarray of double, shape (..., ndim)
        Points to locate

    Returns
    -------
    i : ndarray of int, same shape as `xi`
        Indices of simplices containing each point.
        Points outside the triangulation get the value -1.

    See Also
    --------
    Delaunay.find_simplex

    Notes
    -----
    .. versionadded:: 0.9

    Examples
    --------
    >>> import numpy as np
    >>> import matplotlib.pyplot as plt
    >>> from scipy.spatial import Delaunay, delaunay_plot_2d, tsearch
    >>> rng = np.random.default_rng()

    The Delaunay triangulation of a set of random points:

    >>> pts = rng.random((20, 2))
    >>> tri = Delaunay(pts)
    >>> _ = delaunay_plot_2d(tri)

    Find the simplices containing a given set of points:

    >>> loc = rng.uniform(0.2, 0.8, (5, 2))
    >>> s = tsearch(tri, loc)
    >>> plt.triplot(pts[:, 0], pts[:, 1], tri.simplices[s], 'b-', mask=s==-1)
    >>> plt.scatter(loc[:, 0], loc[:, 1], c='r', marker='x')
    >>> plt.show()

    """
    return tri.find_simplex(xi)

# Set docstring for foo to docstring of bar, working around change in Cython 0.28
# See https://github.com/scipy/scipy/pull/8581
def _copy_docstr(dst, src):
    try:
        dst.__doc__ = src.__doc__
    except AttributeError:
        dst.__func__.__doc__ = src.__func__.__doc__


_copy_docstr(Delaunay.add_points, _QhullUser._add_points)

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
        vn_indptr, vn_indices = obj.vertex_neighbor_vertices
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
    good : ndarray of bool or None
        A one-dimensional Boolean array indicating which facets are
        good. Used with options that compute good facets, e.g. QGn
        and QG-n. Good facets are defined as those that are
        visible (n) or invisible (-n) from point n, where
        n is the nth point in 'points'. The 'good' attribute may be
        used as an index into 'simplices' to return the good (visible)
        facets: simplices[good]. A facet is visible from the outside
        of the hull only, and neither coplanarity nor degeneracy count
        as cases of visibility.

        If a "QGn" or "QG-n" option is not specified, None is returned.

        .. versionadded:: 1.3.0
    area : float
        Surface area of the convex hull when input dimension > 2.
        When input `points` are 2-dimensional, this is the perimeter of the convex hull.

        .. versionadded:: 0.17.0
    volume : float
        Volume of the convex hull when input dimension > 2.
        When input `points` are 2-dimensional, this is the area of the convex hull.

        .. versionadded:: 0.17.0

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

    Examples
    --------

    Convex hull of a random set of points:

    >>> from scipy.spatial import ConvexHull, convex_hull_plot_2d
    >>> import numpy as np
    >>> rng = np.random.default_rng()
    >>> points = rng.random((30, 2))   # 30 random points in 2-D
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

    Facets visible from a point:

    Create a square and add a point above the square.

    >>> generators = np.array([[0.2, 0.2],
    ...                        [0.2, 0.4],
    ...                        [0.4, 0.4],
    ...                        [0.4, 0.2],
    ...                        [0.3, 0.6]])

    Call ConvexHull with the QG option. QG4 means
    compute the portions of the hull not including
    point 4, indicating the facets that are visible
    from point 4.

    >>> hull = ConvexHull(points=generators,
    ...                   qhull_options='QG4')

    The "good" array indicates which facets are
    visible from point 4.

    >>> print(hull.simplices)
        [[1 0]
         [1 2]
         [3 0]
         [3 2]]
    >>> print(hull.good)
        [False  True False False]

    Now plot it, highlighting the visible facets.

    >>> fig = plt.figure()
    >>> ax = fig.add_subplot(1,1,1)
    >>> for visible_facet in hull.simplices[hull.good]:
    ...     ax.plot(hull.points[visible_facet, 0],
    ...             hull.points[visible_facet, 1],
    ...             color='violet',
    ...             lw=6)
    >>> convex_hull_plot_2d(hull, ax=ax)
        <Figure size 640x480 with 1 Axes> # may vary
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
            qhull_options = qhull_options.encode('latin1')

        # Run qhull
        qhull = _Qhull(b"i", points, qhull_options, required_options=b"Qt",
                       incremental=incremental)
        _QhullUser.__init__(self, qhull, incremental=incremental)

    def _update(self, qhull):
        qhull.triangulate()

        self.simplices, self.neighbors, self.equations, self.coplanar, self.good = \
                       qhull.get_simplex_facet_array()

        # only populate self.good with a QG option set
        option_set = set()
        if qhull.options is not None:
            option_set.update(qhull.options.split())

        QG_option_present = 0
        for option in option_set:
            if b"QG" in option:
                QG_option_present += 1
                break

        if not QG_option_present:
            self.good = None
        else:
            self.good = self.good.astype(bool)

        self.volume, self.area = qhull.volume_area()

        if qhull.ndim == 2:
            self._vertices = qhull.get_extremes_2d()
        else:
            self._vertices = None

        self.nsimplex = self.simplices.shape[0]

        _QhullUser._update(self, qhull)

    def add_points(self, points, restart=False):
        self._add_points(points, restart)

    @property
    def points(self):
        return self._points

    @property
    def vertices(self):
        if self._vertices is None:
            self._vertices = np.unique(self.simplices)
        return self._vertices

_copy_docstr(ConvexHull.add_points, _QhullUser._add_points)

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
        Coordinates of points to construct a Voronoi diagram from
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
        When qhull option "Qz" was specified, an empty sublist
        represents the Voronoi region for a point at infinity that
        was added internally.
    point_region : array of ints, shape (npoints)
        Index of the Voronoi region for each input point.
        If qhull option "Qc" was not specified, the list will contain -1
        for points that are not associated with a Voronoi region.
        If qhull option "Qz" was specified, there will be one less
        element than the number of regions because an extra point
        at infinity is added internally to facilitate computation.
    furthest_site
        True if this was a furthest site triangulation and False if not.

        .. versionadded:: 1.4.0

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

    Examples
    --------
    Voronoi diagram for a set of point:

    >>> import numpy as np
    >>> points = np.array([[0, 0], [0, 1], [0, 2], [1, 0], [1, 1], [1, 2],
    ...                    [2, 0], [2, 1], [2, 2]])
    >>> from scipy.spatial import Voronoi, voronoi_plot_2d
    >>> vor = Voronoi(points)

    Plot it:

    >>> import matplotlib.pyplot as plt
    >>> fig = voronoi_plot_2d(vor)
    >>> plt.show()

    The Voronoi vertices:

    >>> vor.vertices
    array([[0.5, 0.5],
           [0.5, 1.5],
           [1.5, 0.5],
           [1.5, 1.5]])

    There is a single finite Voronoi region, and four finite Voronoi
    ridges:

    >>> vor.regions
    [[], [-1, 0], [-1, 1], [1, -1, 0], [3, -1, 2], [-1, 3], [-1, 2], [0, 1, 3, 2], [2, -1, 0], [3, -1, 1]]
    >>> vor.ridge_vertices
    [[-1, 0], [-1, 0], [-1, 1], [-1, 1], [0, 1], [-1, 3], [-1, 2], [2, 3], [-1, 3], [-1, 2], [1, 3], [0, 2]]

    The ridges are perpendicular between lines drawn between the following
    input points:

    >>> vor.ridge_points
    array([[0, 3],
           [0, 1],
           [2, 5],
           [2, 1],
           [1, 4],
           [7, 8],
           [7, 6],
           [7, 4],
           [8, 5],
           [6, 3],
           [4, 5],
           [4, 3]], dtype=int32)

    """
    def __init__(self, points, furthest_site=False, incremental=False,
                 qhull_options=None):
        if np.ma.isMaskedArray(points):
            raise ValueError('Input points cannot be a masked array')
        points = np.ascontiguousarray(points, dtype=np.double)
        if points.ndim != 2:
            raise ValueError("Input points array must have 2 dimensions.")

        if qhull_options is None:
            if not incremental:
                qhull_options = b"Qbb Qc Qz"
            else:
                qhull_options = b"Qc"
            if points.shape[1] >= 5:
                qhull_options += b" Qx"
        else:
            qhull_options = qhull_options.encode('latin1')

        # Run qhull
        qhull = _Qhull(b"v", points, qhull_options, furthest_site=furthest_site,
                       incremental=incremental)
        _QhullUser.__init__(self, qhull, incremental=incremental)

        self.furthest_site = furthest_site

    def _update(self, qhull):
        self.vertices, self.ridge_points, self.ridge_vertices, \
                       self.regions, self.point_region = \
                       qhull.get_voronoi_diagram()

        self._ridge_dict = None

        _QhullUser._update(self, qhull)

    def add_points(self, points, restart=False):
        self._add_points(points, restart)

    @property
    def points(self):
        return self._points

    @property
    def ridge_dict(self):
        if self._ridge_dict is None:
            self._ridge_dict = dict(zip(map(tuple, self.ridge_points.tolist()),
                                        self.ridge_vertices))
        return self._ridge_dict


_copy_docstr(Voronoi.add_points, _QhullUser._add_points)

#------------------------------------------------------------------------------
# Halfspace Intersection
#------------------------------------------------------------------------------

class HalfspaceIntersection(_QhullUser):
    """
    HalfspaceIntersection(halfspaces, interior_point, incremental=False, qhull_options=None)

    Halfspace intersections in N dimensions.

    .. versionadded:: 0.19.0

    Parameters
    ----------
    halfspaces : ndarray of floats, shape (nineq, ndim+1)
        Stacked Inequalities of the form Ax + b <= 0 in format [A; b]
    interior_point : ndarray of floats, shape (ndim,)
        Point clearly inside the region defined by halfspaces. Also called a feasible
        point, it can be obtained by linear programming.
    incremental : bool, optional
        Allow adding new halfspaces incrementally. This takes up some additional
        resources.
    qhull_options : str, optional
        Additional options to pass to Qhull. See Qhull manual
        for details. (Default: "Qx" for ndim > 4 and "" otherwise)
        Option "H" is always enabled.

    Attributes
    ----------
    halfspaces : ndarray of double, shape (nineq, ndim+1)
        Input halfspaces.
    interior_point :ndarray of floats, shape (ndim,)
        Input interior point.
    intersections : ndarray of double, shape (ninter, ndim)
        Intersections of all halfspaces.
    dual_points : ndarray of double, shape (nineq, ndim)
        Dual points of the input halfspaces.
    dual_facets : list of lists of ints
        Indices of points forming the (non necessarily simplicial) facets of
        the dual convex hull.
    dual_vertices : ndarray of ints, shape (nvertices,)
        Indices of halfspaces forming the vertices of the dual convex hull.
        For 2-D convex hulls, the vertices are in counterclockwise order.
        For other dimensions, they are in input order.
    dual_equations : ndarray of double, shape (nfacet, ndim+1)
        [normal, offset] forming the hyperplane equation of the dual facet
        (see `Qhull documentation <http://www.qhull.org/>`__  for more).
    dual_area : float
        Area of the dual convex hull
    dual_volume : float
        Volume of the dual convex hull

    Raises
    ------
    QhullError
        Raised when Qhull encounters an error condition, such as
        geometrical degeneracy when options to resolve are not enabled.
    ValueError
        Raised if an incompatible array is given as input.

    Notes
    -----
    The intersections are computed using the
    `Qhull library <http://www.qhull.org/>`__.
    This reproduces the "qhalf" functionality of Qhull.

    Examples
    --------

    Halfspace intersection of planes forming some polygon

    >>> from scipy.spatial import HalfspaceIntersection
    >>> import numpy as np
    >>> halfspaces = np.array([[-1, 0., 0.],
    ...                        [0., -1., 0.],
    ...                        [2., 1., -4.],
    ...                        [-0.5, 1., -2.]])
    >>> feasible_point = np.array([0.5, 0.5])
    >>> hs = HalfspaceIntersection(halfspaces, feasible_point)

    Plot halfspaces as filled regions and intersection points:

    >>> import matplotlib.pyplot as plt
    >>> fig = plt.figure()
    >>> ax = fig.add_subplot(1, 1, 1, aspect='equal')
    >>> xlim, ylim = (-1, 3), (-1, 3)
    >>> ax.set_xlim(xlim)
    >>> ax.set_ylim(ylim)
    >>> x = np.linspace(-1, 3, 100)
    >>> symbols = ['-', '+', 'x', '*']
    >>> signs = [0, 0, -1, -1]
    >>> fmt = {"color": None, "edgecolor": "b", "alpha": 0.5}
    >>> for h, sym, sign in zip(halfspaces, symbols, signs):
    ...     hlist = h.tolist()
    ...     fmt["hatch"] = sym
    ...     if h[1]== 0:
    ...         ax.axvline(-h[2]/h[0], label='{}x+{}y+{}=0'.format(*hlist))
    ...         xi = np.linspace(xlim[sign], -h[2]/h[0], 100)
    ...         ax.fill_between(xi, ylim[0], ylim[1], **fmt)
    ...     else:
    ...         ax.plot(x, (-h[2]-h[0]*x)/h[1], label='{}x+{}y+{}=0'.format(*hlist))
    ...         ax.fill_between(x, (-h[2]-h[0]*x)/h[1], ylim[sign], **fmt)
    >>> x, y = zip(*hs.intersections)
    >>> ax.plot(x, y, 'o', markersize=8)

    By default, qhull does not provide with a way to compute an interior point.
    This can easily be computed using linear programming. Considering halfspaces
    of the form :math:`Ax + b \leq 0`, solving the linear program:

    .. math::

        max \: y

        s.t. Ax + y ||A_i|| \leq -b

    With :math:`A_i` being the rows of A, i.e. the normals to each plane.

    Will yield a point x that is furthest inside the convex polyhedron. To
    be precise, it is the center of the largest hypersphere of radius y
    inscribed in the polyhedron. This point is called the Chebyshev center
    of the polyhedron (see [1]_ 4.3.1, pp148-149). The
    equations outputted by Qhull are always normalized.

    >>> from scipy.optimize import linprog
    >>> from matplotlib.patches import Circle
    >>> norm_vector = np.reshape(np.linalg.norm(halfspaces[:, :-1], axis=1),
    ...     (halfspaces.shape[0], 1))
    >>> c = np.zeros((halfspaces.shape[1],))
    >>> c[-1] = -1
    >>> A = np.hstack((halfspaces[:, :-1], norm_vector))
    >>> b = - halfspaces[:, -1:]
    >>> res = linprog(c, A_ub=A, b_ub=b, bounds=(None, None))
    >>> x = res.x[:-1]
    >>> y = res.x[-1]
    >>> circle = Circle(x, radius=y, alpha=0.3)
    >>> ax.add_patch(circle)
    >>> plt.legend(bbox_to_anchor=(1.6, 1.0))
    >>> plt.show()

    References
    ----------
    .. [Qhull] http://www.qhull.org/
    .. [1] S. Boyd, L. Vandenberghe, Convex Optimization, available
           at http://stanford.edu/~boyd/cvxbook/

    """

    def __init__(self, halfspaces, interior_point,
                    incremental=False, qhull_options=None):
        if np.ma.isMaskedArray(halfspaces):
            raise ValueError('Input halfspaces cannot be a masked array')
        if np.ma.isMaskedArray(interior_point):
            raise ValueError('Input interior point cannot be a masked array')
        if interior_point.shape != (halfspaces.shape[1]-1,):
            raise ValueError('Feasible point must be a (ndim-1,) array')
        halfspaces = np.ascontiguousarray(halfspaces, dtype=np.double)
        self.interior_point = np.ascontiguousarray(interior_point, dtype=np.double)

        if qhull_options is None:
            qhull_options = b""
            if halfspaces.shape[1] >= 6:
                qhull_options += b"Qx"
        else:
            qhull_options = qhull_options.encode('latin1')

        # Run qhull
        mode_option = "H"
        qhull = _Qhull(mode_option.encode(), halfspaces, qhull_options,
                       required_options=None, incremental=incremental,
                       interior_point=self.interior_point)

        _QhullUser.__init__(self, qhull, incremental=incremental)

    def _update(self, qhull):
        self.dual_facets, self.dual_equations = qhull.get_hull_facets()

        self.dual_points = qhull.get_hull_points()

        self.dual_volume, self.dual_area = qhull.volume_area()

        self.intersections = self.dual_equations[:, :-1]/-self.dual_equations[:, -1:] + self.interior_point

        if qhull.ndim == 2:
            self._vertices = qhull.get_extremes_2d()
        else:
            self._vertices = None

        _QhullUser._update(self, qhull)

        self.ndim = self.halfspaces.shape[1] - 1
        self.nineq = self.halfspaces.shape[0]

    def add_halfspaces(self, halfspaces, restart=False):
        """
        add_halfspaces(halfspaces, restart=False)

        Process a set of additional new halfspaces.

        Parameters
        ----------
        halfspaces : ndarray
            New halfspaces to add. The dimensionality should match that of the
            initial halfspaces.
        restart : bool, optional
            Whether to restart processing from scratch, rather than
            adding halfspaces incrementally.

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
        object to be able to add halfspaces incrementally. Incremental addition
        of halfspaces is also not possible after `close` has been called.

        """
        self._add_points(halfspaces, restart, self.interior_point)

    @property
    def halfspaces(self):
        return self._points

    @property
    def dual_vertices(self):
        if self._vertices is None:
            self._vertices = np.unique(np.array(self.dual_facets))
        return self._vertices
