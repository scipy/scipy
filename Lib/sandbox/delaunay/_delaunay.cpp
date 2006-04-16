#include "Python.h"
#include <stdlib.h>
#include <map>

#include "VoronoiDiagramGenerator.h"
#include "delaunay_utils.h"
#include "natneighbors.h"
#include "numpy/arrayobject.h"

using namespace std;

extern "C" {

static void reorder_edges(int npoints, int ntriangles, 
    double *x, double *y, 
    int *edge_db, int *tri_edges, int *tri_nbrs)
{
    int neighbors[3], nodes[3];
    int i, tmp;
    int case1, case2;

    for (i=0; i<ntriangles; i++) {
        nodes[0] = INDEX2(edge_db, INDEX3(tri_edges,i,0), 0);
        nodes[1] = INDEX2(edge_db, INDEX3(tri_edges,i,0), 1);
        tmp = INDEX2(edge_db, INDEX3(tri_edges,i,1), 0);
        if (tmp == nodes[0]) {
            case1 = 1;
            nodes[2] = INDEX2(edge_db, INDEX3(tri_edges,i,1), 1);
        } else if (tmp == nodes[1]) {
            case1 = 0;
            nodes[2] = INDEX2(edge_db, INDEX3(tri_edges,i,1), 1);
        } else if (INDEX2(edge_db, INDEX3(tri_edges,i,1), 1) == nodes[0]) {
            case1 = 1;
            nodes[2] = tmp;
        } else {
            case1 = 0;
            nodes[2] = tmp;
        }

        if (ONRIGHT(x[nodes[0]], y[nodes[0]],
                    x[nodes[1]], y[nodes[1]],
                    x[nodes[2]], y[nodes[2]])) {
            // flip to make counter-clockwise
            tmp = nodes[2];
            nodes[2] = nodes[1];
            nodes[1] = tmp;
            case2 = 1;
        } else case2 = 0;

        // I worked it out on paper. You're just gonna have to trust me on this.
        if (!case1 && !case2) {
            neighbors[0] = INDEX3(tri_nbrs, i, 1);
            neighbors[1] = INDEX3(tri_nbrs, i, 2);
            neighbors[2] = INDEX3(tri_nbrs, i, 0);
        } else if (case1 && !case2) {
            neighbors[0] = INDEX3(tri_nbrs, i, 2);
            neighbors[1] = INDEX3(tri_nbrs, i, 1);
            neighbors[2] = INDEX3(tri_nbrs, i, 0);
        } else if (!case1 && case2) {
            neighbors[0] = INDEX3(tri_nbrs, i, 1);
            neighbors[1] = INDEX3(tri_nbrs, i, 0);
            neighbors[2] = INDEX3(tri_nbrs, i, 2);
        } else {
            neighbors[0] = INDEX3(tri_nbrs, i, 2);
            neighbors[1] = INDEX3(tri_nbrs, i, 0);
            neighbors[2] = INDEX3(tri_nbrs, i, 1);
        }

        // Not trusting me? Okay, let's go through it:
        // We have three edges to deal with and three nodes. Without loss
        // of generality, let's label the nodes A, B, and C with (A, B) 
        // forming the first edge in the order they arrive on input.
        // Then there are eight possibilities as to how the other edge-tuples
        // may be labeled, but only two variations that are going to affect the
        // output:
        //
        // AB         AB
        // BC (CB)    AC (CA)
        // CA (AC)    BC (CB)
        //
        // The distinction is whether A is in the second edge or B is.
        // This is the test "case1" above.
        //
        // The second test we need to perform is for counter-clockwiseness.
        // Again, there are only two variations that will affect the outcome:
        // either ABC is counter-clockwise, or it isn't. In the former case,
        // we're done setting the node order, we just need to associate the 
        // appropriate neighbor triangles with their opposite nodes, something
        // which can be done by inspection. In the latter case, to order the
        // nodes counter-clockwise, we only have to switch B and C to get
        // nodes ACB. Then we simply set the neighbor list by inspection again.
        //
        //          CCW     CW
        // AB
        // BC       120     102 -+
        // CA                    |
        //                       +- neighbor order
        // AB                    |
        // AC       210     201 -+
        // BC
        //          ABC     ACB -+- node order


        INDEX3(tri_edges,i,0) = nodes[0];
        INDEX3(tri_edges,i,1) = nodes[1];
        INDEX3(tri_edges,i,2) = nodes[2];
        INDEX3(tri_nbrs,i,0) = neighbors[0];
        INDEX3(tri_nbrs,i,1) = neighbors[1];
        INDEX3(tri_nbrs,i,2) = neighbors[2];
    }
}

static PyObject* getMesh(int npoints, double *x, double *y)
{
    PyObject *vertices, *edge_db, *tri_edges, *tri_nbrs;
    PyObject *temp;
    int tri0, tri1, reg0, reg1;
    double tri0x, tri0y, tri1x, tri1y;
    int length, numtri, i, j;
    intp dim[MAX_DIMS];
    int *edge_db_ptr, *tri_edges_ptr, *tri_nbrs_ptr;
    double *vertices_ptr;
    VoronoiDiagramGenerator vdg;

    vdg.generateVoronoi(x, y, npoints, -100, 100, -100, 100, 0);
    vdg.getNumbers(length, numtri);

    dim[0] = length;
    dim[1] = 2;
    edge_db = PyArray_SimpleNew(2, dim, PyArray_INT);
    if (!edge_db) goto fail;
    edge_db_ptr = (int*)PyArray_DATA(edge_db);
    
    dim[0] = numtri;
    vertices = PyArray_SimpleNew(2, dim, PyArray_DOUBLE);
    if (!vertices) goto fail;
    vertices_ptr = (double*)PyArray_DATA(vertices);

    dim[1] = 3;
    tri_edges = PyArray_SimpleNew(2, dim, PyArray_INT);
    if (!tri_edges) goto fail;
    tri_edges_ptr = (int*)PyArray_DATA(tri_edges);

    tri_nbrs = PyArray_SimpleNew(2, dim, PyArray_INT);
    if (!tri_nbrs) goto fail;
    tri_nbrs_ptr = (int*)PyArray_DATA(tri_nbrs);

    for (i=0; i<(3*numtri); i++) {
        tri_edges_ptr[i] = tri_nbrs_ptr[i] = -1;
    }

    vdg.resetEdgeListIter();
    i = -1;
    while (vdg.getNextDelaunay(tri0, tri0x, tri0y, tri1, tri1x, tri1y, reg0, reg1)) {
        i++;
        INDEX2(edge_db_ptr,i,0) = reg0;
        INDEX2(edge_db_ptr,i,1) = reg1;
        if (tri0 > -1) {
            INDEX2(vertices_ptr,tri0,0) = tri0x;
            INDEX2(vertices_ptr,tri0,1) = tri0y;
            for (j=0; j<3; j++) {
                if (INDEX3(tri_edges_ptr,tri0,j) == i) break;
                if (INDEX3(tri_edges_ptr,tri0,j) == -1) {
                    INDEX3(tri_edges_ptr,tri0,j) = i;
                    INDEX3(tri_nbrs_ptr,tri0,j) = tri1;
                    break;
                }
            }
        }
        if (tri1 > -1) {
            INDEX2(vertices_ptr,tri1,0) = tri1x;
            INDEX2(vertices_ptr,tri1,1) = tri1y;
            for (j=0; j<3; j++) {
                if (INDEX3(tri_edges_ptr,tri1,j) == i) break;
                if (INDEX3(tri_edges_ptr,tri1,j) == -1) {
                    INDEX3(tri_edges_ptr,tri1,j) = i;
                    INDEX3(tri_nbrs_ptr,tri1,j) = tri0;
                    break;
                }
            }
        }
    }

    // tri_edges contains lists of edges; convert to lists of nodes in
    // counterclockwise order and reorder tri_nbrs to match. Each node
    // corresponds to the edge opposite it in the triangle.
    reorder_edges(npoints, numtri, x, y, edge_db_ptr, tri_edges_ptr, 
        tri_nbrs_ptr);

    temp = Py_BuildValue("(OOOO)", vertices, edge_db, tri_edges, tri_nbrs);
    if (!temp) goto fail;

    Py_DECREF(vertices);
    Py_DECREF(edge_db);
    Py_DECREF(tri_edges);
    Py_DECREF(tri_nbrs);

    return temp;

fail:
    Py_XDECREF(vertices);
    Py_XDECREF(edge_db);
    Py_XDECREF(tri_edges);
    Py_XDECREF(tri_nbrs);
    return NULL;
}

static PyObject *linear_planes(int ntriangles, double *x, double *y, double *z,
    int *nodes)
{
    intp dims[2];
    PyObject *planes;
    int i;
    double *planes_ptr;
    double x02, y02, z02, x12, y12, z12, xy0212;
    
    dims[0] = ntriangles;
    dims[1] = 3;
    planes = PyArray_SimpleNew(2, dims, PyArray_DOUBLE);
    if (!planes) return NULL;
    planes_ptr = (double *)PyArray_DATA(planes);

    for (i=0; i<ntriangles; i++) {
        x02 = x[INDEX3(nodes,i,0)] - x[INDEX3(nodes,i,2)];
        y02 = y[INDEX3(nodes,i,0)] - y[INDEX3(nodes,i,2)];
        z02 = z[INDEX3(nodes,i,0)] - z[INDEX3(nodes,i,2)];
        x12 = x[INDEX3(nodes,i,1)] - x[INDEX3(nodes,i,2)];
        y12 = y[INDEX3(nodes,i,1)] - y[INDEX3(nodes,i,2)];
        z12 = z[INDEX3(nodes,i,1)] - z[INDEX3(nodes,i,2)];

        if (y12 != 0.0) {
            xy0212 = y02/y12;
            INDEX3(planes_ptr,i,0) = (z02 - z12 * xy0212) / (x02 - x12 * xy0212);
            INDEX3(planes_ptr,i,1) = (z12 - INDEX3(planes_ptr,i,0)*x12) / y12;
            INDEX3(planes_ptr,i,2) = (z[INDEX3(nodes,i,2)] - 
                                      INDEX3(planes_ptr,i,0)*x[INDEX3(nodes,i,2)] - 
                                      INDEX3(planes_ptr,i,1)*y[INDEX3(nodes,i,2)]);
        } else {
            xy0212 = x02/x12;
            INDEX3(planes_ptr,i,1) = (z02 - z12 * xy0212) / (y02 - y12 * xy0212);
            INDEX3(planes_ptr,i,0) = (z12 - INDEX3(planes_ptr,i,1)*y12) / x12;
            INDEX3(planes_ptr,i,2) = (z[INDEX3(nodes,i,2)] - 
                                      INDEX3(planes_ptr,i,0)*x[INDEX3(nodes,i,2)] - 
                                      INDEX3(planes_ptr,i,1)*y[INDEX3(nodes,i,2)]);
        }
    }

    return (PyObject*)planes;
}

static double linear_interpolate_single(double targetx, double targety, 
    double *x, double *y, int *nodes, int *neighbors,
    PyObject *planes, double defvalue, int start_triangle, int *end_triangle)
{
    double *planes_ptr;
    planes_ptr = (double*)PyArray_DATA(planes);
    if (start_triangle == -1) start_triangle = 0;
    *end_triangle = walking_triangles(start_triangle, targetx, targety, 
        x, y, nodes, neighbors);
    if (*end_triangle == -1) return defvalue;
    return (targetx*INDEX3(planes_ptr,*end_triangle,0) + 
            targety*INDEX3(planes_ptr,*end_triangle,1) +
                    INDEX3(planes_ptr,*end_triangle,2));
}

static PyObject *linear_interpolate_grid(double x0, double x1, int xsteps, 
    double y0, double y1, int ysteps,
    PyObject *planes, double defvalue, 
    int npoints, double *x, double *y, int *nodes, int *neighbors)
{
    int ix, iy;
    double dx, dy, targetx, targety;
    int rowtri, coltri, tri;
    PyObject *z;
    double *z_ptr;
    intp dims[2];

    dims[0] = ysteps;
    dims[1] = xsteps;
    z = PyArray_SimpleNew(2, dims, PyArray_DOUBLE);
    if (!z) return NULL;
    z_ptr = (double*)PyArray_DATA(z);

    dx = (x1 - x0) / (xsteps-1);
    dy = (y1 - y0) / (ysteps-1);

    rowtri = 0;
    for (iy=0; iy<ysteps; iy++) {
        targety = y0 + dy*iy;
        rowtri = walking_triangles(rowtri, x0, targety, x, y, nodes, neighbors);
        tri = rowtri;
        for (ix=0; ix<xsteps; ix++) {
            targetx = x0 + dx*ix;
            INDEXN(z_ptr, xsteps, iy, ix) = linear_interpolate_single(
                targetx, targety,
                x, y, nodes, neighbors, planes, defvalue, tri, &coltri);
            if (coltri != -1) tri = coltri;
        }
    }

    return z;
}

static PyObject *compute_planes_method(PyObject *self, PyObject *args)
{
    PyObject *pyx, *pyy, *pyz, *pynodes;
    PyArrayObject *x, *y, *z, *nodes;
    int npoints, ntriangles;

    PyObject *planes;

    if (!PyArg_ParseTuple(args, "OOOO", &pyx, &pyy, &pyz, &pynodes)) {
        return NULL;
    }
    x = (PyArrayObject*)PyArray_ContiguousFromObject(pyx, PyArray_DOUBLE, 1, 1);
    y = (PyArrayObject*)PyArray_ContiguousFromObject(pyy, PyArray_DOUBLE, 1, 1);
    z = (PyArrayObject*)PyArray_ContiguousFromObject(pyz, PyArray_DOUBLE, 1, 1);
    if ((!x || !y) || !z) {
        PyErr_SetString(PyExc_ValueError, "x,y,z must be 1-D arrays of floats");
        goto fail;
    }
    npoints = PyArray_DIM(x, 0);
    if ((PyArray_DIM(y, 0) != npoints) || (PyArray_DIM(z, 0) != npoints)) {
        PyErr_SetString(PyExc_ValueError, "x,y,z arrays must be of equal length");
        goto fail;
    }
    nodes = (PyArrayObject*)PyArray_ContiguousFromObject(pynodes, PyArray_INT, 2, 2);
    if (!nodes) {
        PyErr_SetString(PyExc_ValueError, "nodes must be a 2-D array of ints");
        goto fail;
    }
    ntriangles = PyArray_DIM(nodes, 0);
    if (PyArray_DIM(nodes, 1) != 3) {
        PyErr_SetString(PyExc_ValueError, "nodes must have shape (ntriangles, 3)");
        goto fail;
    }

    planes = linear_planes(ntriangles, (double*)PyArray_DATA(x), 
        (double*)PyArray_DATA(y), (double*)PyArray_DATA(z), (int*)PyArray_DATA(nodes));

    Py_DECREF(pyx);
    Py_DECREF(pyy);
    Py_DECREF(pyz);
    Py_DECREF(pynodes);

    return planes;

fail:
    Py_XDECREF(pyx);
    Py_XDECREF(pyy);
    Py_XDECREF(pyz);
    Py_XDECREF(pynodes);
    return NULL;
}

static PyObject *linear_interpolate_method(PyObject *self, PyObject *args)
{
    double x0, x1, y0, y1, defvalue;
    int xsteps, ysteps;
    PyObject *pyplanes, *pyx, *pyy, *pynodes, *pyneighbors, *grid;
    int npoints;

    PyArrayObject *planes, *x, *y, *nodes, *neighbors;

    if (!PyArg_ParseTuple(args, "ddlddldOOOOO", &x0, &x1, &xsteps, &y0, &y1, &ysteps,
           &defvalue, &pyplanes, &pyx, &pyy, &pynodes, &pyneighbors)) {
        return NULL;
    }
    x = (PyArrayObject*)PyArray_ContiguousFromObject(pyx, PyArray_DOUBLE, 1, 1);
    y = (PyArrayObject*)PyArray_ContiguousFromObject(pyy, PyArray_DOUBLE, 1, 1);
    if (!x || !y) {
        PyErr_SetString(PyExc_ValueError, "x,y must be 1-D arrays of floats");
        goto fail;
    }
    npoints = PyArray_DIM(x, 0);
    if (PyArray_DIM(y, 0) != npoints) {
        PyErr_SetString(PyExc_ValueError, "x,y arrays must be of equal length");
        goto fail;
    }
    planes = (PyArrayObject*)PyArray_ContiguousFromObject(pyplanes, PyArray_DOUBLE, 2, 2);
    if (!planes) {
        PyErr_SetString(PyExc_ValueError, "planes must be a 2-D array of floats");
        goto fail;
    }
    nodes = (PyArrayObject*)PyArray_ContiguousFromObject(pynodes, PyArray_INT, 2, 2);
    if (!nodes) {
        PyErr_SetString(PyExc_ValueError, "nodes must be a 2-D array of ints");
        goto fail;
    }
    neighbors = (PyArrayObject*)PyArray_ContiguousFromObject(pyneighbors, PyArray_INT, 2, 2);
    if (!neighbors) {
        PyErr_SetString(PyExc_ValueError, "neighbors must be a 2-D array of ints");
        goto fail;
    }

    grid = linear_interpolate_grid(x0, x1, xsteps, y0, y1, ysteps,
        (PyObject*)planes, defvalue, npoints,
        (double*)PyArray_DATA(x), (double*)PyArray_DATA(y),
        (int*)PyArray_DATA(nodes), (int*)PyArray_DATA(neighbors));
    
    Py_DECREF(x);
    Py_DECREF(y);
    Py_DECREF(planes);
    Py_DECREF(nodes);
    Py_DECREF(neighbors);

    return grid;

fail:
    Py_XDECREF(x);
    Py_XDECREF(y);
    Py_XDECREF(planes);
    Py_XDECREF(nodes);
    Py_XDECREF(neighbors);
    return NULL;
}

#define CLEANUP \
    Py_XDECREF(x);\
    Py_XDECREF(y);\
    Py_XDECREF(z);\
    Py_XDECREF(intx);\
    Py_XDECREF(inty);\
    Py_XDECREF(centers);\
    Py_XDECREF(nodes);\
    Py_XDECREF(neighbors);\
    Py_XDECREF(intz);\
    return NULL;

static PyObject *nn_interpolate_unstructured_method(PyObject *self, PyObject *args)
{
    PyObject *pyx, *pyy, *pyz, *pycenters, *pynodes, *pyneighbors, *pyintx, *pyinty;
    PyArrayObject *x, *y, *z, *centers, *nodes, *neighbors, *intx, *inty, *intz;
    double defvalue;
    int size, npoints, ntriangles;

    if (!PyArg_ParseTuple(args, "OOdOOOOOO", &pyintx, &pyinty, &defvalue,
        &pyx, &pyy, &pyz, &pycenters, &pynodes, &pyneighbors)) {
        return NULL;
    }
    x = (PyArrayObject*)PyArray_ContiguousFromObject(pyx, PyArray_DOUBLE, 1, 1);
    y = (PyArrayObject*)PyArray_ContiguousFromObject(pyy, PyArray_DOUBLE, 1, 1);
    z = (PyArrayObject*)PyArray_ContiguousFromObject(pyz, PyArray_DOUBLE, 1, 1);
    if ((!x || !y) || !z) {
        PyErr_SetString(PyExc_ValueError, "x,y must be 1-D arrays of floats");
        //goto fail;
        CLEANUP
    }
    npoints = PyArray_DIM(x, 0);
    if ((PyArray_DIM(y, 0) != npoints) || (PyArray_DIM(z, 0) != npoints)) {
        PyErr_SetString(PyExc_ValueError, "x,y,z arrays must be of equal length");
        //goto fail;
        CLEANUP
    }
    centers = (PyArrayObject*)PyArray_ContiguousFromObject(pycenters, PyArray_DOUBLE, 2, 2);
    if (!centers) {
        PyErr_SetString(PyExc_ValueError, "centers must be a 2-D array of ints");
        //goto fail;
        CLEANUP
    }
    nodes = (PyArrayObject*)PyArray_ContiguousFromObject(pynodes, PyArray_INT, 2, 2);
    if (!nodes) {
        PyErr_SetString(PyExc_ValueError, "nodes must be a 2-D array of ints");
        //goto fail;
        CLEANUP
    }
    neighbors = (PyArrayObject*)PyArray_ContiguousFromObject(pyneighbors, PyArray_INT, 2, 2);
    if (!neighbors) {
        PyErr_SetString(PyExc_ValueError, "neighbors must be a 2-D array of ints");
        //goto fail;
        CLEANUP
    }
    ntriangles = PyArray_DIM(neighbors, 0);
    if ((PyArray_DIM(nodes, 0) != ntriangles)  || 
        (PyArray_DIM(centers, 0) != ntriangles)) {
        PyErr_SetString(PyExc_ValueError, "centers,nodes,neighbors must be of equal length");
        //goto fail;
        CLEANUP
    }
    intx = (PyArrayObject*)PyArray_ContiguousFromObject(pyintx, PyArray_DOUBLE, 0, 0);
    inty = (PyArrayObject*)PyArray_ContiguousFromObject(pyinty, PyArray_DOUBLE, 0, 0);
    if (!intx || !inty) {
        PyErr_SetString(PyExc_ValueError, "intx,inty must be arrays of floats");
        CLEANUP
    }
    if (intx->nd != inty->nd) {
        PyErr_SetString(PyExc_ValueError, "intx,inty must have same shapes");
        CLEANUP
    }
    for (int i=0; i<intx->nd; i++) {
        if (intx->dimensions[i] != inty->dimensions[i]) {
            PyErr_SetString(PyExc_ValueError, "intx,inty must have same shapes");
            CLEANUP
        }
    }
    intz = (PyArrayObject*)PyArray_SimpleNew(intx->nd, intx->dimensions, PyArray_DOUBLE);
    if (!intz) {CLEANUP}

    NaturalNeighbors nn(npoints, ntriangles, 
        (double*)PyArray_DATA(x), (double*)PyArray_DATA(y),
        (double*)PyArray_DATA(centers), (int*)PyArray_DATA(nodes), 
        (int*)PyArray_DATA(neighbors));
    size = PyArray_Size((PyObject*)intx);
    nn.interpolate_unstructured((double*)PyArray_DATA(z), size, 
        (double*)PyArray_DATA(intx), (double*)PyArray_DATA(inty), 
        (double*)PyArray_DATA(intz), defvalue);

    return (PyObject*)intz;
}

#undef CLEANUP

#define CLEANUP \
    Py_XDECREF(x);\
    Py_XDECREF(y);\
    Py_XDECREF(z);\
    Py_XDECREF(centers);\
    Py_XDECREF(nodes);\
    Py_XDECREF(neighbors);\
    return NULL;

static PyObject *nn_interpolate_method(PyObject *self, PyObject *args)
{
    PyObject *pyx, *pyy, *pyz, *pycenters, *pynodes, *pyneighbors, *grid;
    PyArrayObject *x, *y, *z, *centers, *nodes, *neighbors;
    double x0, x1, y0, y1, defvalue;
    double *grid_ptr;
    int xsteps, ysteps;
    int npoints, ntriangles;
    intp dims[2];

    if (!PyArg_ParseTuple(args, "ddlddldOOOOOO", &x0, &x1, &xsteps, 
        &y0, &y1, &ysteps, &defvalue, &pyx, &pyy, &pyz, &pycenters, &pynodes,
        &pyneighbors)) {
        return NULL;
    }
    x = (PyArrayObject*)PyArray_ContiguousFromObject(pyx, PyArray_DOUBLE, 1, 1);
    y = (PyArrayObject*)PyArray_ContiguousFromObject(pyy, PyArray_DOUBLE, 1, 1);
    z = (PyArrayObject*)PyArray_ContiguousFromObject(pyz, PyArray_DOUBLE, 1, 1);
    if ((!x || !y) || !z) {
        PyErr_SetString(PyExc_ValueError, "x,y must be 1-D arrays of floats");
        //goto fail;
        CLEANUP
    }
    npoints = PyArray_DIM(x, 0);
    if (PyArray_DIM(y, 0) != npoints) {
        PyErr_SetString(PyExc_ValueError, "x,y arrays must be of equal length");
        //goto fail;
        CLEANUP
    }
    centers = (PyArrayObject*)PyArray_ContiguousFromObject(pycenters, PyArray_DOUBLE, 2, 2);
    if (!centers) {
        PyErr_SetString(PyExc_ValueError, "centers must be a 2-D array of ints");
        //goto fail;
        CLEANUP
    }
    nodes = (PyArrayObject*)PyArray_ContiguousFromObject(pynodes, PyArray_INT, 2, 2);
    if (!nodes) {
        PyErr_SetString(PyExc_ValueError, "nodes must be a 2-D array of ints");
        //goto fail;
        CLEANUP
    }
    neighbors = (PyArrayObject*)PyArray_ContiguousFromObject(pyneighbors, PyArray_INT, 2, 2);
    if (!neighbors) {
        PyErr_SetString(PyExc_ValueError, "neighbors must be a 2-D array of ints");
        //goto fail;
        CLEANUP
    }
    ntriangles = PyArray_DIM(neighbors, 0);
    if ((PyArray_DIM(nodes, 0) != ntriangles)  || 
        (PyArray_DIM(centers, 0) != ntriangles)) {
        PyErr_SetString(PyExc_ValueError, "centers,nodes,neighbors must be of equal length");
        //goto fail;
        CLEANUP
    }
    //goto succeed; // XXX: Can't cross NaturalNeighbors instantiation with goto

    dims[0] = ysteps;
    dims[1] = xsteps;
    grid = PyArray_SimpleNew(2, dims, PyArray_DOUBLE);
    if (!grid) {CLEANUP} // goto fail;
    grid_ptr = (double*)PyArray_DATA(grid);

    NaturalNeighbors nn(npoints, ntriangles, 
        (double*)PyArray_DATA(x), (double*)PyArray_DATA(y),
        (double*)PyArray_DATA(centers), (int*)PyArray_DATA(nodes), 
        (int*)PyArray_DATA(neighbors));
    nn.interpolate_grid((double*)PyArray_DATA(z), 
        x0, x1, xsteps,
        y0, y1, ysteps,
        (double*)PyArray_DATA(grid),
        defvalue, 0);

    Py_DECREF(x);
    Py_DECREF(y);
    Py_DECREF(z);
    Py_DECREF(centers);
    Py_DECREF(nodes);
    Py_DECREF(neighbors);

    return grid;

}
#undef CLEANUP

static PyObject *delaunay_method(PyObject *self, PyObject *args)
{
    PyObject *pyx, *pyy, *mesh;
    PyArrayObject *x, *y;
    int npoints;

    if (!PyArg_ParseTuple(args, "OO", &pyx, &pyy)) {
        return NULL;
    }
    x = (PyArrayObject*)PyArray_ContiguousFromObject(pyx, PyArray_DOUBLE, 1, 1);
    y = (PyArrayObject*)PyArray_ContiguousFromObject(pyy, PyArray_DOUBLE, 1, 1);
    if (!x || !y) {
        PyErr_SetString(PyExc_ValueError, "x,y must be 1-D arrays of floats");
        goto fail;
    }

    npoints = PyArray_DIM(x, 0);
    if (PyArray_DIM(y, 0) != npoints) {
        PyErr_SetString(PyExc_ValueError, "x and y must have the same length");
        goto fail;
    }

    mesh = getMesh(npoints, (double*)(x->data), (double*)(y->data));

    if (!mesh) goto fail;

    Py_DECREF(x);
    Py_DECREF(y);

    return mesh;

fail:
    Py_XDECREF(x);
    Py_XDECREF(y);
    return NULL;
}

static PyMethodDef delaunay_methods[] = {
    {"delaunay", (PyCFunction)delaunay_method, METH_VARARGS, 
        "Compute the Delaunay triangulation of a cloud of 2-D points.\n\n"
        "circumcenters, edges, tri_points, tri_neighbors = delaunay(x, y)\n\n"
        "x, y -- shape-(npoints,) arrays of floats giving the X and Y coordinates of the points\n"
        "circumcenters -- shape-(numtri,2) array of floats giving the coordinates of the\n"
        "    circumcenters of each triangle (numtri being the number of triangles)\n"
        "edges -- shape-(nedges,2) array of integers giving the indices into x and y\n"
        "    of each edge in the triangulation\n"
        "tri_points -- shape-(numtri,3) array of integers giving the indices into x and y\n"
        "    of each node in each triangle\n"
        "tri_neighbors -- shape-(numtri,3) array of integers giving the indices into circumcenters\n"
        "    tri_points, and tri_neighbors of the neighbors of each triangle\n"},
    {"compute_planes", (PyCFunction)compute_planes_method, METH_VARARGS,
        ""},
    {"linear_interpolate_grid", (PyCFunction)linear_interpolate_method, METH_VARARGS,
        ""},
    {"nn_interpolate_grid", (PyCFunction)nn_interpolate_method, METH_VARARGS,
        ""},
    {"nn_interpolate_unstructured", (PyCFunction)nn_interpolate_unstructured_method, METH_VARARGS,
        ""},
    {NULL, NULL, 0, NULL}
};


DL_EXPORT(void) init_delaunay()
{
    PyObject* m;
    m = Py_InitModule3("_delaunay", delaunay_methods, 
        "Tools for computing the Delaunay triangulation and some operations on it.\n"
        );
    if (m == NULL)
        return;
    import_array();
}

} // extern "C"
