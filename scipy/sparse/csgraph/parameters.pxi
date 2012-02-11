DTYPE = np.float64
ctypedef np.float64_t DTYPE_t

ITYPE = np.int32
ctypedef np.int32_t ITYPE_t

# NULL_IDX is the index used in predecessor matrices to store a non-path
cdef ITYPE_t NULL_IDX = -9999
