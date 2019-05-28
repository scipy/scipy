import numpy as np
cimport numpy as np
cimport cython
from libc.math cimport (sin, acos, atan2,
                        cos, M_PI, abs)
from libc.stdlib cimport RAND_MAX, rand

@cython.cdivision(True)
@cython.boundscheck(False)
@cython.wraparound(False)
cpdef _geometric_slerp(double[:] start_coord,
                       double[:] end_coord,
                       double[:] t_values):
   cdef:
      double omega = acos(np.dot(start_coord, end_coord))
      double sin_omega = sin(omega)
      int ndims = start_coord.shape[0]
      n_pts = t_values.size
      double[:,:] new_pts = np.empty((n_pts, ndims), dtype=np.float64)
      int i, j
      double factors[2]

   for i in range(n_pts):
      factors[0] = sin((1 - t_values[i]) * omega) / sin_omega
      factors[1] = sin(t_values[i] * omega) / sin_omega
      for j in range(ndims):
          new_pts[i,j] = ((factors[0] * start_coord[j]) +
                          (factors[1] * end_coord[j]))
   return np.asarray(new_pts)
