import numpy as np
cimport numpy as np
cimport cython

__all__ = ['sort_vertices_of_regions']

cdef void remaining_filter(int shape, long[:] remaining, int current_simplex):
   cdef unsigned int i
   for i in range(shape):
     if remaining[i] == current_simplex:
       remaining[i] = -2
    
@cython.boundscheck(False)
def sort_vertices_of_regions(int[:,::1] simplices, regions):
  cdef unsigned int n, k, s, i
  cdef int num_regions = len(regions)
  cdef int current_simplex, current_vertex, remaining_count
  cdef int cs_identified, remaining_size
  cdef long[:] remaining
  cdef long[:] sorted_vertices = np.zeros(max([len(region) for region in
	  regions]), dtype = np.int64)

  for n in range(num_regions):
    remaining_count = 0
    remaining = np.asarray(regions[n][:])
    remaining_size = remaining.shape[0]
    sorted_vertices[:] = -2
    current_simplex = remaining[0]
    for i in range(3):
      k = simplices[current_simplex][i]
      if k != n:
        current_vertex = k
        break
    sorted_vertices[remaining_count] = current_simplex
    remaining_count += 1
    remaining_filter(remaining.shape[0], remaining, current_simplex)
    while remaining_size > remaining_count:
      cs_identified = 0
      for i in range(remaining_size):
        if remaining[i] == -2:
          continue
        s = remaining[i]
        for j in range(3):
          if current_vertex == simplices[s,j]:
            current_simplex = remaining[i]
            cs_identified += 1
            break
        if cs_identified > 0:
          break
      for i in range(3):
        s = simplices[current_simplex, i]
        if s != n and s != current_vertex:
          current_vertex = s
          break
      sorted_vertices[remaining_count] = current_simplex
      remaining_count += 1
      remaining_filter(remaining.shape[0], remaining, current_simplex)
    regions_arr = np.asarray(sorted_vertices)
    regions[n] = regions_arr[regions_arr > -2].tolist()
