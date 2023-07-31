# distutils: language=c++
# cython: language_level=3

from .HighsExceptions cimport (
    create_highs_exceptions,
    PresolveException,
    highs_exception_handler,
)

# Create exceptions
create_highs_exceptions()

# Create Python level exception
cdef class PyPresolveException(Exception):
    pass

def get_PresolveException():
    return PyPresolveException

def handle_highs_exceptions():
    highs_exception_handler()
