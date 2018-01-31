import numpy as np
cimport numpy as np
cimport cython
from libc.math cimport (sin, acos, atan2,
                        cos, M_PI)

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

@cython.boundscheck(False)
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
        area += delta_lambda * sin(phi_vals[i])
    area = (radius ** 2) * area
    return area

def _slerp(double[:] start_coord,
           double[:] end_coord,
           int n_pts):
    # spherical linear interpolation between points
    # on great circle arc
    # see: https://en.wikipedia.org/wiki/Slerp#Geometric_Slerp
    # NOTE: could we use scipy.interpolate.RectSphereBivariateSpline instead?
    cdef double[:] t_values
    cdef double omega
    cdef double[:,:] new_pts = np.empty((n_pts, 3), dtype=np.float64)
    cdef double[:] new_pt = np.empty(3)
    cdef int i, j
    cdef double t

    omega = acos(np.dot(start_coord, end_coord))
    t_values = np.linspace(0, 1, n_pts)
    for i in range(n_pts):
        t = t_values[i]
        for j in range(3):
            new_pt[j] = (((sin((1 - t) * omega) / sin(omega)) * start_coord[j]) +
                      ((sin(t * omega) / sin(omega)) * end_coord[j]))
        new_pts[i] = new_pt
    return new_pts

def calc_heading(double[:] lambda_range,
                 double[:] phi_range,
                 int i,
                 int next_index):

    cdef double lambda_1 = lambda_range[i]
    cdef double lambda_2 = lambda_range[next_index]
    cdef double phi_1 = phi_range[i]
    cdef double phi_2 = phi_range[next_index]
    cdef double delta_lambda = lambda_1 - lambda_2
    cdef double term_1 = sin(delta_lambda) * cos(phi_2)
    cdef double term_2 = cos(phi_1) * sin(phi_2)
    cdef double term_3 = sin(phi_1) * cos(phi_2) * cos(delta_lambda)
    cdef double result = atan2(term_1, term_2 - term_3)
    cdef double course_angle = result % (2 * M_PI)
    course_angle = course_angle * 180.0 / M_PI
    return course_angle
