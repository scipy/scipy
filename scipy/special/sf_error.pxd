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

    ctypedef enum sf_action_t:
        IGNORE "SF_ERROR_IGNORE"
        WARN "SF_ERROR_WARN"
        RAISE "SF_ERROR_RAISE"

    char **sf_error_messages
    void error "sf_error" (char *func_name, sf_error_t code, char *fmt, ...) nogil
    void check_fpe "sf_error_check_fpe" (char *func_name) nogil
    void set_action "scipy_sf_error_set_action" (sf_error_t code, sf_action_t action) nogil
    sf_action_t get_action "scipy_sf_error_get_action" (sf_error_t code) nogil


cdef inline int _sf_error_test_function(int code) noexcept nogil:
    """Function that can raise every sf_error category for testing
    purposes.

    """
    cdef sf_error_t sf_error_code
    
    if code < 0 or code >= 10:
        sf_error_code = OTHER
    else:
        sf_error_code = <sf_error_t>code
    error('_err_test_function', sf_error_code, NULL)
    return 0
