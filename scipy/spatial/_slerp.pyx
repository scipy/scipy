import numpy as np
cimport numpy as np
cimport cython
from libc.math cimport (sin, acos, atan2,
                        cos, M_PI, abs)
from libc.stdlib cimport RAND_MAX, rand

@cython.cdivision(True)
@cython.boundscheck(False)
@cython.wraparound(False)
cpdef _geometric_slerp(double[:] start,
                       double[:] end,
                       double[:] t):
   cdef:
      double omega = acos(np.dot(start, end))
      double sin_omega = sin(omega)
      int ndims = start.shape[0]
      n_pts = t.size
      double[:,:] new_pts = np.empty((n_pts, ndims), dtype=np.float64)
      int i, j
      double factors[2]

   for i in range(n_pts):
      factors[0] = sin((1 - t[i]) * omega) / sin_omega
      factors[1] = sin(t[i] * omega) / sin_omega
      for j in range(ndims):
          new_pts[i,j] = ((factors[0] * start[j]) +
                          (factors[1] * end[j]))
   return np.asarray(new_pts)
