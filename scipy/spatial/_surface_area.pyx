import numpy as np
cimport numpy as np
cimport cython

cdef int vertex_index_strider(int index, int num_vertices):
    cdef int forward_index
    forward_index = index + 1
    if forward_index > (num_vertices - 1):
        forward_index = 0
    return forward_index

@cython.boundscheck(False)
@cython.wraparound(False)
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
    return area

def spherical_polygon_area(double[:,:] vertices, double radius):
    cdef double [:] lambda_vals = np.arctan2(vertices[...,1],
                                             vertices[...,0])
    cdef double [:] phi_vals = np.arcsin(np.asarray(vertices[...,2]) / radius)
    cdef int N = vertices.shape[0]
    cdef int i, forward_index, backward_index
    cdef double area = 0
    cdef double delta_lambda

    for i in range(N):
        forward_index = vertex_index_strider(i, N)
        backward_index = i - 1
        delta_lambda = (lambda_vals[forward_index] -
                        lambda_vals[backward_index])
        area += delta_lambda * np.sin(phi_vals[i])
    area = (radius ** 2) * area
    return area
