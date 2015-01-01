from ._fortran import *


# Don't use deprecation Numpy C API.  Define this to a fixed version instead of
# NPY_API_VERSION in order not to break compilation for released Scipy versions
# when Numpy introduces a new deprecation.  Use in setup.py::
#
#   config.add_extension('_name', sources=['source_fname'], **numpy_nodepr_api)
#
numpy_nodepr_api = dict(define_macros=[("NPY_NO_DEPRECATED_API",
                                        "NPY_1_9_API_VERSION")])
