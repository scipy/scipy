""" Converters for all of NumPy's scalar types such as
    int32, float32, complex128, etc.
"""
import numpy
import c_spec

class numpy_complex_scalar_converter(c_spec.complex_converter):
    """ Handles conversion of all the NumPy complex types.
        This uses the same machinery as the standard python
        complex converter.
    """
    def init_info(self):
        # First, set up all the same specifications the normal
        # complex converter uses.
        c_spec.complex_converter.init_info(self)

        # But set this converter up to match the numpy complex
        # types.
        self.matching_types = numpy.sctypes['complex']
