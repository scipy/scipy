# -*-cython-*-

cdef extern from "sf_error.h":
     ctypedef enum sf_error_t:
          OK "SF_ERROR_OK"
          SINGULAR "SF_ERROR_SINGULAR"
          UNDERFLOW "SF_ERROR_UNDERFLOW"
          OVERFLOW "SF_ERROR_OVERFLOW"
          SLOW "SF_ERROR_SLOW"
          LOSS "SF_ERROR_LOSS"
          NO_RESULT "SF_ERROR_NO_RESULT"
          DOMAIN "SF_ERROR_DOMAIN"
          ARG "SF_ERROR_ARG"
          OTHER "SF_ERROR_OTHER"

     void error "sf_error" (char *func_name, sf_error_t code, char *fmt, ...) nogil
