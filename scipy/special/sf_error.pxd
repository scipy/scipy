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

     char **sf_error_messages
     void error "sf_error" (char *func_name, sf_error_t code, char *fmt, ...) nogil
     void check_fpe "sf_error_check_fpe" (char *func_name) nogil
     int set_print "sf_error_set_print" (int flag) nogil
     int get_print "sf_error_get_print" () nogil
