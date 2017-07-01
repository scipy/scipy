from __future__ import division, print_function, absolute_import

import itertools

import numpy as np
from scipy.ndimage._ni_iterators import LineBufferIterator
from numpy.testing import assert_array_equal, run_module_suite


class TestLineBufferIterator:

    valid_types = '?' + np.typecodes['AllInteger'] + 'fd'
    input_array = np.arange(3 * 5 * 7).reshape(3, 5, 7)

    def all_type_pairs(self):
        return itertools.product(self.valid_types, repeat=2)

    def test_all_types(self):
        for input_type, output_type in self.all_type_pairs():
            input_array = self.input_array.astype(input_type)
            output_array = np.empty_like(input_array, dtype=output_type)
            for axis in range(input_array.ndim):
                lbiter = LineBufferIterator(input_array, output_array, axis, 1)
                for input_buffer, output_buffer in lbiter:
                    output_buffer[:] = input_buffer
                if output_type == '?':
                    input_array = input_array.astype('?')
                assert_array_equal(input_array, output_array,
                                   '{} -> {}'.format(input_type, output_type))


if __name__ == '__main__':
    run_module_suite()
