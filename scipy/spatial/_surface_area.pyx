import numpy as np
cimport numpy as np
cimport cython
from libc.math cimport abs

cdef int vertex_index_strider(int index, int num_vertices):
    cdef int forward_index
    forward_index = index + 1
    if forward_index > (num_vertices - 1):
        forward_index = 0
    return forward_index

@cython.boundscheck(False)
def planar_polygon_area(double[:,::1] vertices):
    cdef int N = vertices.shape[0]
    cdef int i, forward_index, backward_index
    cdef double area = 0
    cdef double delta_x

    for i in range(N):
        forward_index = vertex_index_strider(i, N)
        backward_index = i - 1
        delta_x = (vertices[forward_index][0] -
                   vertices[backward_index][0])
        area += delta_x * vertices[i][1]
    area *= 0.5
    return abs(area)
