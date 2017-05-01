from numpy cimport npy_intp as intp

cdef int _filter1d(double *input_line, intp input_length, double *output_line,
	           intp output_length, void *callback_data)
cdef int _filter2d(double *buffer, intp filter_size, double *res,
	           void *callback_data)
cdef int _transform(intp *output_coordinates, double *input_coordinates,
	            int output_rank, int input_rank, void *callback_data)
