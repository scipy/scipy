from __future__ import division, print_function, absolute_import

from numpy.testing import assert_, assert_array_almost_equal, assert_equal, \
                          assert_almost_equal, assert_array_equal, \
                          assert_raises, run_module_suite, TestCase
import numpy as np

import scipy.ndimage as ndimage

import os.path

types = [np.int8, np.uint8, np.int16,
         np.uint16, np.int32, np.uint32,
         np.int64, np.uint64,
         np.float32, np.float64]


np.mod(1., 1)  # Silence fmod bug on win-amd64. See #1408 and #1238.


class Test_measurements_stats(TestCase):
    """ndimage.measurements._stats() is a utility function used by other functions."""

    def test_a(self):
        x = [0,1,2,6]
        labels = [0,0,1,1]
        index = [0,1]
        for shp in [(4,), (2,2)]:
            x = np.array(x).reshape(shp)
            labels = np.array(labels).reshape(shp)
            counts, sums = ndimage.measurements._stats(x, labels=labels, index=index)
            assert_array_equal(counts, [2, 2])
            assert_array_equal(sums, [1.0, 8.0])

    def test_b(self):
        # Same data as test_a, but different labels.  The label 9 exceeds the
        # length of 'labels', so this test will follow a different code path.
        x = [0,1,2,6]
        labels = [0,0,9,9]
        index = [0,9]
        for shp in [(4,), (2,2)]:
            x = np.array(x).reshape(shp)
            labels = np.array(labels).reshape(shp)
            counts, sums = ndimage.measurements._stats(x, labels=labels, index=index)
            assert_array_equal(counts, [2, 2])
            assert_array_equal(sums, [1.0, 8.0])

    def test_a_centered(self):
        x = [0,1,2,6]
        labels = [0,0,1,1]
        index = [0,1]
        for shp in [(4,), (2,2)]:
            x = np.array(x).reshape(shp)
            labels = np.array(labels).reshape(shp)
            counts, sums, centers = ndimage.measurements._stats(x, labels=labels,
                                        index=index, centered=True)
            assert_array_equal(counts, [2, 2])
            assert_array_equal(sums, [1.0, 8.0])
            assert_array_equal(centers, [0.5, 8.0])

    def test_b_centered(self):
        x = [0,1,2,6]
        labels = [0,0,9,9]
        index = [0,9]
        for shp in [(4,), (2,2)]:
            x = np.array(x).reshape(shp)
            labels = np.array(labels).reshape(shp)
            counts, sums, centers = ndimage.measurements._stats(x, labels=labels,
                                        index=index, centered=True)
            assert_array_equal(counts, [2, 2])
            assert_array_equal(sums, [1.0, 8.0])
            assert_array_equal(centers, [0.5, 8.0])

    def test_nonint_labels(self):
        x = [0,1,2,6]
        labels = [0.0, 0.0, 9.0, 9.0]
        index = [0.0, 9.0]
        for shp in [(4,), (2,2)]:
            x = np.array(x).reshape(shp)
            labels = np.array(labels).reshape(shp)
            counts, sums, centers = ndimage.measurements._stats(x, labels=labels,
                                        index=index, centered=True)
            assert_array_equal(counts, [2, 2])
            assert_array_equal(sums, [1.0, 8.0])
            assert_array_equal(centers, [0.5, 8.0])


class Test_measurements_select(TestCase):
    """ndimage.measurements._select() is a utility function used by other functions."""

    def test_basic(self):
        x = [0,1,6,2]
        cases = [
            ([0,0,1,1], [0,1]),                 # "Small" integer labels
            ([0,0,9,9], [0,9]),                 # A label larger than len(labels)
            ([0.0,0.0,7.0,7.0], [0.0, 7.0]),    # Non-integer labels
        ]
        for labels, index in cases:
            result = ndimage.measurements._select(x, labels=labels, index=index)
            assert_(len(result) == 0)
            result = ndimage.measurements._select(x, labels=labels, index=index, find_max=True)
            assert_(len(result) == 1)
            assert_array_equal(result[0], [1, 6])
            result = ndimage.measurements._select(x, labels=labels, index=index, find_min=True)
            assert_(len(result) == 1)
            assert_array_equal(result[0], [0, 2])
            result = ndimage.measurements._select(x, labels=labels, index=index,
                                find_min=True, find_min_positions=True)
            assert_(len(result) == 2)
            assert_array_equal(result[0], [0, 2])
            assert_array_equal(result[1], [0, 3])
            result = ndimage.measurements._select(x, labels=labels, index=index,
                                find_max=True, find_max_positions=True)
            assert_(len(result) == 2)
            assert_array_equal(result[0], [1, 6])
            assert_array_equal(result[1], [1, 2])


def test_label01():
    "label 1"
    data = np.ones([])
    out, n = ndimage.label(data)
    assert_array_almost_equal(out, 1)
    assert_equal(n, 1)


def test_label02():
    "label 2"
    data = np.zeros([])
    out, n = ndimage.label(data)
    assert_array_almost_equal(out, 0)
    assert_equal(n, 0)


def test_label03():
    "label 3"
    data = np.ones([1])
    out, n = ndimage.label(data)
    assert_array_almost_equal(out, [1])
    assert_equal(n, 1)


def test_label04():
    "label 4"
    data = np.zeros([1])
    out, n = ndimage.label(data)
    assert_array_almost_equal(out, [0])
    assert_equal(n, 0)


def test_label05():
    "label 5"
    data = np.ones([5])
    out, n = ndimage.label(data)
    assert_array_almost_equal(out, [1, 1, 1, 1, 1])
    assert_equal(n, 1)


def test_label06():
    "label 6"
    data = np.array([1, 0, 1, 1, 0, 1])
    out, n = ndimage.label(data)
    assert_array_almost_equal(out, [1, 0, 2, 2, 0, 3])
    assert_equal(n, 3)


def test_label07():
    "label 7"
    data = np.array([[0, 0, 0, 0, 0, 0],
                           [0, 0, 0, 0, 0, 0],
                           [0, 0, 0, 0, 0, 0],
                           [0, 0, 0, 0, 0, 0],
                           [0, 0, 0, 0, 0, 0],
                           [0, 0, 0, 0, 0, 0]])
    out, n = ndimage.label(data)
    assert_array_almost_equal(out, [[0, 0, 0, 0, 0, 0],
                               [0, 0, 0, 0, 0, 0],
                               [0, 0, 0, 0, 0, 0],
                               [0, 0, 0, 0, 0, 0],
                               [0, 0, 0, 0, 0, 0],
                               [0, 0, 0, 0, 0, 0]])
    assert_equal(n, 0)


def test_label08():
    "label 8"
    data = np.array([[1, 0, 0, 0, 0, 0],
                           [0, 0, 1, 1, 0, 0],
                           [0, 0, 1, 1, 1, 0],
                           [1, 1, 0, 0, 0, 0],
                           [1, 1, 0, 0, 0, 0],
                           [0, 0, 0, 1, 1, 0]])
    out, n = ndimage.label(data)
    assert_array_almost_equal(out, [[1, 0, 0, 0, 0, 0],
                               [0, 0, 2, 2, 0, 0],
                               [0, 0, 2, 2, 2, 0],
                               [3, 3, 0, 0, 0, 0],
                               [3, 3, 0, 0, 0, 0],
                               [0, 0, 0, 4, 4, 0]])
    assert_equal(n, 4)


def test_label09():
    "label 9"
    data = np.array([[1, 0, 0, 0, 0, 0],
                           [0, 0, 1, 1, 0, 0],
                           [0, 0, 1, 1, 1, 0],
                           [1, 1, 0, 0, 0, 0],
                           [1, 1, 0, 0, 0, 0],
                           [0, 0, 0, 1, 1, 0]])
    struct = ndimage.generate_binary_structure(2, 2)
    out, n = ndimage.label(data, struct)
    assert_array_almost_equal(out, [[1, 0, 0, 0, 0, 0],
                               [0, 0, 2, 2, 0, 0],
                               [0, 0, 2, 2, 2, 0],
                               [2, 2, 0, 0, 0, 0],
                               [2, 2, 0, 0, 0, 0],
                               [0, 0, 0, 3, 3, 0]])
    assert_equal(n, 3)


def test_label10():
    "label 10"
    data = np.array([[0, 0, 0, 0, 0, 0],
                           [0, 1, 1, 0, 1, 0],
                           [0, 1, 1, 1, 1, 0],
                           [0, 0, 0, 0, 0, 0]])
    struct = ndimage.generate_binary_structure(2, 2)
    out, n = ndimage.label(data, struct)
    assert_array_almost_equal(out, [[0, 0, 0, 0, 0, 0],
                               [0, 1, 1, 0, 1, 0],
                               [0, 1, 1, 1, 1, 0],
                               [0, 0, 0, 0, 0, 0]])
    assert_equal(n, 1)


def test_label11():
    "label 11"
    for type in types:
        data = np.array([[1, 0, 0, 0, 0, 0],
                               [0, 0, 1, 1, 0, 0],
                               [0, 0, 1, 1, 1, 0],
                               [1, 1, 0, 0, 0, 0],
                               [1, 1, 0, 0, 0, 0],
                               [0, 0, 0, 1, 1, 0]], type)
        out, n = ndimage.label(data)
        expected = [[1, 0, 0, 0, 0, 0],
                    [0, 0, 2, 2, 0, 0],
                    [0, 0, 2, 2, 2, 0],
                    [3, 3, 0, 0, 0, 0],
                    [3, 3, 0, 0, 0, 0],
                    [0, 0, 0, 4, 4, 0]]
        assert_array_almost_equal(out, expected)
        assert_equal(n, 4)


def test_label11_inplace():
    "label 11 in place"
    for type in types:
        data = np.array([[1, 0, 0, 0, 0, 0],
                               [0, 0, 1, 1, 0, 0],
                               [0, 0, 1, 1, 1, 0],
                               [1, 1, 0, 0, 0, 0],
                               [1, 1, 0, 0, 0, 0],
                               [0, 0, 0, 1, 1, 0]], type)
        n = ndimage.label(data, output=data)
        expected = [[1, 0, 0, 0, 0, 0],
                    [0, 0, 2, 2, 0, 0],
                    [0, 0, 2, 2, 2, 0],
                    [3, 3, 0, 0, 0, 0],
                    [3, 3, 0, 0, 0, 0],
                    [0, 0, 0, 4, 4, 0]]
        assert_array_almost_equal(data, expected)
        assert_equal(n, 4)


def test_label12():
    "label 12"
    for type in types:
        data = np.array([[0, 0, 0, 0, 1, 1],
                               [0, 0, 0, 0, 0, 1],
                               [0, 0, 1, 0, 1, 1],
                               [0, 0, 1, 1, 1, 1],
                               [0, 0, 0, 1, 1, 0]], type)
        out, n = ndimage.label(data)
        expected = [[0, 0, 0, 0, 1, 1],
                    [0, 0, 0, 0, 0, 1],
                    [0, 0, 1, 0, 1, 1],
                    [0, 0, 1, 1, 1, 1],
                    [0, 0, 0, 1, 1, 0]]
        assert_array_almost_equal(out, expected)
        assert_equal(n, 1)


def test_label13():
    "label 13"
    for type in types:
        data = np.array([[1, 0, 1, 1, 1, 0, 1, 1, 1, 0, 1],
                               [1, 1, 1, 0, 1, 1, 1, 0, 1, 1, 1],
                               [1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1],
                               [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1]],
                              type)
        out, n = ndimage.label(data)
        expected = [[1, 0, 1, 1, 1, 0, 1, 1, 1, 0, 1],
                    [1, 1, 1, 0, 1, 1, 1, 0, 1, 1, 1],
                    [1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1],
                    [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1]]
        assert_array_almost_equal(out, expected)
        assert_equal(n, 1)


def test_label_output_typed():
    "test label with specified output with type"
    data = np.ones([5])
    for t in types:
        output = np.zeros([5], dtype=t)
        n = ndimage.label(data, output=output)
        assert_array_almost_equal(output, 1)
        assert_equal(n, 1)


def test_label_output_dtype():
    "test label with specified output dtype"
    data = np.ones([5])
    for t in types:
        output, n = ndimage.label(data, output=t)
        assert_array_almost_equal(output, 1)
        assert output.dtype == t


def test_label_output_wrong_size():
    "test label with output of wrong size"
    data = np.ones([5])
    for t in types:
        output = np.zeros([10], t)
        assert_raises((RuntimeError, ValueError), ndimage.label, data, output=output)


def test_label_structuring_elements():
    "test label with different structuring element neighborhoods"
    data = np.loadtxt(os.path.join(os.path.dirname(__file__), "data", "label_inputs.txt"))
    strels = np.loadtxt(os.path.join(os.path.dirname(__file__), "data", "label_strels.txt"))
    results = np.loadtxt(os.path.join(os.path.dirname(__file__), "data", "label_results.txt"))
    data = data.reshape((-1, 7, 7))
    strels = strels.reshape((-1, 3, 3))
    results = results.reshape((-1, 7, 7))
    r = 0
    for i in range(data.shape[0]):
        d = data[i, :, :]
        for j in range(strels.shape[0]):
            s = strels[j, :, :]
            assert_equal(ndimage.label(d, s)[0], results[r, :, :])
            r += 1


def test_label_default_dtype():
    test_array = np.random.rand(10, 10)
    label, no_features = ndimage.label(test_array > 0.5)
    assert_(label.dtype in (np.int32, np.int64))
    # Shouldn't raise an exception
    ndimage.find_objects(label)


def test_find_objects01():
    "find_objects 1"
    data = np.ones([], dtype=int)
    out = ndimage.find_objects(data)
    assert_(out == [()])


def test_find_objects02():
    "find_objects 2"
    data = np.zeros([], dtype=int)
    out = ndimage.find_objects(data)
    assert_(out == [])


def test_find_objects03():
    "find_objects 3"
    data = np.ones([1], dtype=int)
    out = ndimage.find_objects(data)
    assert_equal(out, [(slice(0, 1, None),)])


def test_find_objects04():
    "find_objects 4"
    data = np.zeros([1], dtype=int)
    out = ndimage.find_objects(data)
    assert_equal(out, [])


def test_find_objects05():
    "find_objects 5"
    data = np.ones([5], dtype=int)
    out = ndimage.find_objects(data)
    assert_equal(out, [(slice(0, 5, None),)])


def test_find_objects06():
    "find_objects 6"
    data = np.array([1, 0, 2, 2, 0, 3])
    out = ndimage.find_objects(data)
    assert_equal(out, [(slice(0, 1, None),),
                       (slice(2, 4, None),),
                       (slice(5, 6, None),)])


def test_find_objects07():
    "find_objects 7"
    data = np.array([[0, 0, 0, 0, 0, 0],
                        [0, 0, 0, 0, 0, 0],
                        [0, 0, 0, 0, 0, 0],
                        [0, 0, 0, 0, 0, 0],
                        [0, 0, 0, 0, 0, 0],
                        [0, 0, 0, 0, 0, 0]])
    out = ndimage.find_objects(data)
    assert_equal(out, [])


def test_find_objects08():
    "find_objects 8"
    data = np.array([[1, 0, 0, 0, 0, 0],
                           [0, 0, 2, 2, 0, 0],
                           [0, 0, 2, 2, 2, 0],
                           [3, 3, 0, 0, 0, 0],
                           [3, 3, 0, 0, 0, 0],
                           [0, 0, 0, 4, 4, 0]])
    out = ndimage.find_objects(data)
    assert_equal(out, [(slice(0, 1, None), slice(0, 1, None)),
                       (slice(1, 3, None), slice(2, 5, None)),
                       (slice(3, 5, None), slice(0, 2, None)),
                       (slice(5, 6, None), slice(3, 5, None))])


def test_find_objects09():
    "find_objects 9"
    data = np.array([[1, 0, 0, 0, 0, 0],
                           [0, 0, 2, 2, 0, 0],
                           [0, 0, 2, 2, 2, 0],
                           [0, 0, 0, 0, 0, 0],
                           [0, 0, 0, 0, 0, 0],
                           [0, 0, 0, 4, 4, 0]])
    out = ndimage.find_objects(data)
    assert_equal(out, [(slice(0, 1, None), slice(0, 1, None)),
                       (slice(1, 3, None), slice(2, 5, None)),
                       None,
                       (slice(5, 6, None), slice(3, 5, None))])


def test_sum01():
    "sum 1"
    for type in types:
        input = np.array([], type)
        output = ndimage.sum(input)
        assert_equal(output, 0.0)


def test_sum02():
    "sum 2"
    for type in types:
        input = np.zeros([0, 4], type)
        output = ndimage.sum(input)
        assert_equal(output, 0.0)


def test_sum03():
    "sum 3"
    for type in types:
        input = np.ones([], type)
        output = ndimage.sum(input)
        assert_almost_equal(output, 1.0)


def test_sum04():
    "sum 4"
    for type in types:
        input = np.array([1, 2], type)
        output = ndimage.sum(input)
        assert_almost_equal(output, 3.0)


def test_sum05():
    "sum 5"
    for type in types:
        input = np.array([[1, 2], [3, 4]], type)
        output = ndimage.sum(input)
        assert_almost_equal(output, 10.0)


def test_sum06():
    "sum 6"
    labels = np.array([], bool)
    for type in types:
        input = np.array([], type)
        output = ndimage.sum(input, labels=labels)
        assert_equal(output, 0.0)


def test_sum07():
    "sum 7"
    labels = np.ones([0, 4], bool)
    for type in types:
        input = np.zeros([0, 4], type)
        output = ndimage.sum(input, labels=labels)
        assert_equal(output, 0.0)


def test_sum08():
    "sum 8"
    labels = np.array([1, 0], bool)
    for type in types:
        input = np.array([1, 2], type)
        output = ndimage.sum(input, labels=labels)
        assert_equal(output, 1.0)


def test_sum09():
    "sum 9"
    labels = np.array([1, 0], bool)
    for type in types:
        input = np.array([[1, 2], [3, 4]], type)
        output = ndimage.sum(input, labels=labels)
        assert_almost_equal(output, 4.0)


def test_sum10():
    "sum 10"
    labels = np.array([1, 0], bool)
    input = np.array([[1, 2], [3, 4]], bool)
    output = ndimage.sum(input, labels=labels)
    assert_almost_equal(output, 2.0)


def test_sum11():
    "sum 11"
    labels = np.array([1, 2], np.int8)
    for type in types:
        input = np.array([[1, 2], [3, 4]], type)
        output = ndimage.sum(input, labels=labels,
                                       index=2)
        assert_almost_equal(output, 6.0)


def test_sum12():
    "sum 12"
    labels = np.array([[1, 2], [2, 4]], np.int8)
    for type in types:
        input = np.array([[1, 2], [3, 4]], type)
        output = ndimage.sum(input, labels=labels,
                                        index=[4, 8, 2])
        assert_array_almost_equal(output, [4.0, 0.0, 5.0])


def test_mean01():
    "mean 1"
    labels = np.array([1, 0], bool)
    for type in types:
        input = np.array([[1, 2], [3, 4]], type)
        output = ndimage.mean(input, labels=labels)
        assert_almost_equal(output, 2.0)


def test_mean02():
    "mean 2"
    labels = np.array([1, 0], bool)
    input = np.array([[1, 2], [3, 4]], bool)
    output = ndimage.mean(input, labels=labels)
    assert_almost_equal(output, 1.0)


def test_mean03():
    "mean 3"
    labels = np.array([1, 2])
    for type in types:
        input = np.array([[1, 2], [3, 4]], type)
        output = ndimage.mean(input, labels=labels,
                                        index=2)
        assert_almost_equal(output, 3.0)


def test_mean04():
    "mean 4"
    labels = np.array([[1, 2], [2, 4]], np.int8)
    olderr = np.seterr(all='ignore')
    try:
        for type in types:
            input = np.array([[1, 2], [3, 4]], type)
            output = ndimage.mean(input, labels=labels,
                                            index=[4, 8, 2])
            assert_array_almost_equal(output[[0,2]], [4.0, 2.5])
            assert_(np.isnan(output[1]))
    finally:
        np.seterr(**olderr)


def test_minimum01():
    "minimum 1"
    labels = np.array([1, 0], bool)
    for type in types:
        input = np.array([[1, 2], [3, 4]], type)
        output = ndimage.minimum(input, labels=labels)
        assert_almost_equal(output, 1.0)


def test_minimum02():
    "minimum 2"
    labels = np.array([1, 0], bool)
    input = np.array([[2, 2], [2, 4]], bool)
    output = ndimage.minimum(input, labels=labels)
    assert_almost_equal(output, 1.0)


def test_minimum03():
    "minimum 3"
    labels = np.array([1, 2])
    for type in types:
        input = np.array([[1, 2], [3, 4]], type)
        output = ndimage.minimum(input, labels=labels,
                                           index=2)
        assert_almost_equal(output, 2.0)


def test_minimum04():
    "minimum 4"
    labels = np.array([[1, 2], [2, 3]])
    for type in types:
        input = np.array([[1, 2], [3, 4]], type)
        output = ndimage.minimum(input, labels=labels,
                                           index=[2, 3, 8])
        assert_array_almost_equal(output, [2.0, 4.0, 0.0])


def test_maximum01():
    "maximum 1"
    labels = np.array([1, 0], bool)
    for type in types:
        input = np.array([[1, 2], [3, 4]], type)
        output = ndimage.maximum(input, labels=labels)
        assert_almost_equal(output, 3.0)


def test_maximum02():
    "maximum 2"
    labels = np.array([1, 0], bool)
    input = np.array([[2, 2], [2, 4]], bool)
    output = ndimage.maximum(input, labels=labels)
    assert_almost_equal(output, 1.0)


def test_maximum03():
    "maximum 3"
    labels = np.array([1, 2])
    for type in types:
        input = np.array([[1, 2], [3, 4]], type)
        output = ndimage.maximum(input, labels=labels,
                                           index=2)
        assert_almost_equal(output, 4.0)


def test_maximum04():
    "maximum 4"
    labels = np.array([[1, 2], [2, 3]])
    for type in types:
        input = np.array([[1, 2], [3, 4]], type)
        output = ndimage.maximum(input, labels=labels,
                                           index=[2, 3, 8])
        assert_array_almost_equal(output, [3.0, 4.0, 0.0])


def test_maximum05():
    "Ticket #501"
    x = np.array([-3,-2,-1])
    assert_equal(ndimage.maximum(x),-1)


def test_median01():
    "median 1"
    a = np.array([[1, 2, 0, 1],
                  [5, 3, 0, 4],
                  [0, 0, 0, 7],
                  [9, 3, 0, 0]])
    labels = np.array([[1, 1, 0, 2],
                       [1, 1, 0, 2],
                       [0, 0, 0, 2],
                       [3, 3, 0, 0]])
    output = ndimage.median(a, labels=labels, index=[1, 2, 3])
    assert_array_almost_equal(output, [2.5, 4.0, 6.0])


def test_median02():
    "median 2"
    a = np.array([[1, 2, 0, 1],
                  [5, 3, 0, 4],
                  [0, 0, 0, 7],
                  [9, 3, 0, 0]])
    output = ndimage.median(a)
    assert_almost_equal(output, 1.0)


def test_median03():
    "median 3"
    a = np.array([[1, 2, 0, 1],
                  [5, 3, 0, 4],
                  [0, 0, 0, 7],
                  [9, 3, 0, 0]])
    labels = np.array([[1, 1, 0, 2],
                       [1, 1, 0, 2],
                       [0, 0, 0, 2],
                       [3, 3, 0, 0]])
    output = ndimage.median(a, labels=labels)
    assert_almost_equal(output, 3.0)


def test_variance01():
    "variance 1"
    olderr = np.seterr(all='ignore')
    try:
        for type in types:
            input = np.array([], type)
            output = ndimage.variance(input)
            assert_(np.isnan(output))
    finally:
        np.seterr(**olderr)


def test_variance02():
    "variance 2"
    for type in types:
        input = np.array([1], type)
        output = ndimage.variance(input)
        assert_almost_equal(output, 0.0)


def test_variance03():
    "variance 3"
    for type in types:
        input = np.array([1, 3], type)
        output = ndimage.variance(input)
        assert_almost_equal(output, 1.0)


def test_variance04():
    "variance 4"
    input = np.array([1, 0], bool)
    output = ndimage.variance(input)
    assert_almost_equal(output, 0.25)


def test_variance05():
    "variance 5"
    labels = [2, 2, 3]
    for type in types:
        input = np.array([1, 3, 8], type)
        output = ndimage.variance(input, labels, 2)
        assert_almost_equal(output, 1.0)


def test_variance06():
    "variance 6"
    labels = [2, 2, 3, 3, 4]
    olderr = np.seterr(all='ignore')
    try:
        for type in types:
            input = np.array([1, 3, 8, 10, 8], type)
            output = ndimage.variance(input, labels, [2, 3, 4])
            assert_array_almost_equal(output, [1.0, 1.0, 0.0])
    finally:
        np.seterr(**olderr)


def test_standard_deviation01():
    "standard deviation 1"
    olderr = np.seterr(all='ignore')
    try:
        for type in types:
            input = np.array([], type)
            output = ndimage.standard_deviation(input)
            assert_(np.isnan(output))
    finally:
        np.seterr(**olderr)


def test_standard_deviation02():
    "standard deviation 2"
    for type in types:
        input = np.array([1], type)
        output = ndimage.standard_deviation(input)
        assert_almost_equal(output, 0.0)


def test_standard_deviation03():
    "standard deviation 3"
    for type in types:
        input = np.array([1, 3], type)
        output = ndimage.standard_deviation(input)
        assert_almost_equal(output, np.sqrt(1.0))


def test_standard_deviation04():
    "standard deviation 4"
    input = np.array([1, 0], bool)
    output = ndimage.standard_deviation(input)
    assert_almost_equal(output, 0.5)


def test_standard_deviation05():
    "standard deviation 5"
    labels = [2, 2, 3]
    for type in types:
        input = np.array([1, 3, 8], type)
        output = ndimage.standard_deviation(input, labels, 2)
        assert_almost_equal(output, 1.0)


def test_standard_deviation06():
    "standard deviation 6"
    labels = [2, 2, 3, 3, 4]
    olderr = np.seterr(all='ignore')
    try:
        for type in types:
            input = np.array([1, 3, 8, 10, 8], type)
            output = ndimage.standard_deviation(input, labels, [2, 3, 4])
            assert_array_almost_equal(output, [1.0, 1.0, 0.0])
    finally:
        np.seterr(**olderr)


def test_standard_deviation07():
    "standard deviation 7"
    labels = [1]
    olderr = np.seterr(all='ignore')
    try:
        for type in types:
            input = np.array([-0.00619519], type)
            output = ndimage.standard_deviation(input, labels, [1])
            assert_array_almost_equal(output, [0])
    finally:
        np.seterr(**olderr)


def test_minimum_position01():
    "minimum position 1"
    labels = np.array([1, 0], bool)
    for type in types:
        input = np.array([[1, 2], [3, 4]], type)
        output = ndimage.minimum_position(input, labels=labels)
        assert_equal(output, (0, 0))


def test_minimum_position02():
    "minimum position 2"
    for type in types:
        input = np.array([[5, 4, 2, 5],
                                [3, 7, 0, 2],
                                [1, 5, 1, 1]], type)
        output = ndimage.minimum_position(input)
        assert_equal(output, (1, 2))


def test_minimum_position03():
    "minimum position 3"
    input = np.array([[5, 4, 2, 5],
                            [3, 7, 0, 2],
                            [1, 5, 1, 1]], bool)
    output = ndimage.minimum_position(input)
    assert_equal(output, (1, 2))


def test_minimum_position04():
    "minimum position 4"
    input = np.array([[5, 4, 2, 5],
                            [3, 7, 1, 2],
                            [1, 5, 1, 1]], bool)
    output = ndimage.minimum_position(input)
    assert_equal(output, (0, 0))


def test_minimum_position05():
    "minimum position 5"
    labels = [1, 2, 0, 4]
    for type in types:
        input = np.array([[5, 4, 2, 5],
                                [3, 7, 0, 2],
                                [1, 5, 2, 3]], type)
        output = ndimage.minimum_position(input, labels)
        assert_equal(output, (2, 0))


def test_minimum_position06():
    "minimum position 6"
    labels = [1, 2, 3, 4]
    for type in types:
        input = np.array([[5, 4, 2, 5],
                                [3, 7, 0, 2],
                                [1, 5, 1, 1]], type)
        output = ndimage.minimum_position(input, labels, 2)
        assert_equal(output, (0, 1))


def test_minimum_position07():
    "minimum position 7"
    labels = [1, 2, 3, 4]
    for type in types:
        input = np.array([[5, 4, 2, 5],
                                [3, 7, 0, 2],
                                [1, 5, 1, 1]], type)
        output = ndimage.minimum_position(input, labels,
                                                    [2, 3])
        assert_equal(output[0], (0, 1))
        assert_equal(output[1], (1, 2))


def test_maximum_position01():
    "maximum position 1"
    labels = np.array([1, 0], bool)
    for type in types:
        input = np.array([[1, 2], [3, 4]], type)
        output = ndimage.maximum_position(input,
                                                    labels=labels)
        assert_equal(output, (1, 0))


def test_maximum_position02():
    "maximum position 2"
    for type in types:
        input = np.array([[5, 4, 2, 5],
                                [3, 7, 8, 2],
                                [1, 5, 1, 1]], type)
        output = ndimage.maximum_position(input)
        assert_equal(output, (1, 2))


def test_maximum_position03():
    "maximum position 3"
    input = np.array([[5, 4, 2, 5],
                            [3, 7, 8, 2],
                            [1, 5, 1, 1]], bool)
    output = ndimage.maximum_position(input)
    assert_equal(output, (0, 0))


def test_maximum_position04():
    "maximum position 4"
    labels = [1, 2, 0, 4]
    for type in types:
        input = np.array([[5, 4, 2, 5],
                                [3, 7, 8, 2],
                                [1, 5, 1, 1]], type)
        output = ndimage.maximum_position(input, labels)
        assert_equal(output, (1, 1))


def test_maximum_position05():
    "maximum position 5"
    labels = [1, 2, 0, 4]
    for type in types:
        input = np.array([[5, 4, 2, 5],
                                [3, 7, 8, 2],
                                [1, 5, 1, 1]], type)
        output = ndimage.maximum_position(input, labels, 1)
        assert_equal(output, (0, 0))


def test_maximum_position06():
    "maximum position 6"
    labels = [1, 2, 0, 4]
    for type in types:
        input = np.array([[5, 4, 2, 5],
                                [3, 7, 8, 2],
                                [1, 5, 1, 1]], type)
        output = ndimage.maximum_position(input, labels,
                                                    [1, 2])
        assert_equal(output[0], (0, 0))
        assert_equal(output[1], (1, 1))


def test_maximum_position07():
    "maximum position 7 - float labels"
    labels = np.array([1.0, 2.5, 0.0, 4.5])
    for type in types:
        input = np.array([[5, 4, 2, 5],
                          [3, 7, 8, 2],
                          [1, 5, 1, 1]], type)
        output = ndimage.maximum_position(input, labels,
                                          [1.0, 4.5])
        assert_equal(output[0], (0, 0))
        assert_equal(output[1], (0, 3))


def test_extrema01():
    "extrema 1"
    labels = np.array([1, 0], bool)
    for type in types:
        input = np.array([[1, 2], [3, 4]], type)
        output1 = ndimage.extrema(input, labels=labels)
        output2 = ndimage.minimum(input, labels=labels)
        output3 = ndimage.maximum(input, labels=labels)
        output4 = ndimage.minimum_position(input,
                                                     labels=labels)
        output5 = ndimage.maximum_position(input,
                                                     labels=labels)
        assert_equal(output1, (output2, output3, output4, output5))


def test_extrema02():
    "extrema 2"
    labels = np.array([1, 2])
    for type in types:
        input = np.array([[1, 2], [3, 4]], type)
        output1 = ndimage.extrema(input, labels=labels,
                                            index=2)
        output2 = ndimage.minimum(input, labels=labels,
                                            index=2)
        output3 = ndimage.maximum(input, labels=labels,
                                            index=2)
        output4 = ndimage.minimum_position(input,
                                            labels=labels, index=2)
        output5 = ndimage.maximum_position(input,
                                            labels=labels, index=2)
        assert_equal(output1, (output2, output3, output4, output5))


def test_extrema03():
    "extrema 3"
    labels = np.array([[1, 2], [2, 3]])
    for type in types:
        input = np.array([[1, 2], [3, 4]], type)
        output1 = ndimage.extrema(input, labels=labels,
                                            index=[2, 3, 8])
        output2 = ndimage.minimum(input, labels=labels,
                                            index=[2, 3, 8])
        output3 = ndimage.maximum(input, labels=labels,
                                            index=[2, 3, 8])
        output4 = ndimage.minimum_position(input,
                                    labels=labels, index=[2, 3, 8])
        output5 = ndimage.maximum_position(input,
                                    labels=labels, index=[2, 3, 8])
        assert_array_almost_equal(output1[0], output2)
        assert_array_almost_equal(output1[1], output3)
        assert_array_almost_equal(output1[2], output4)
        assert_array_almost_equal(output1[3], output5)


def test_extrema04():
    "extrema 4"
    labels = [1, 2, 0, 4]
    for type in types:
        input = np.array([[5, 4, 2, 5],
                                [3, 7, 8, 2],
                                [1, 5, 1, 1]], type)
        output1 = ndimage.extrema(input, labels, [1, 2])
        output2 = ndimage.minimum(input, labels, [1, 2])
        output3 = ndimage.maximum(input, labels, [1, 2])
        output4 = ndimage.minimum_position(input, labels,
                                                     [1, 2])
        output5 = ndimage.maximum_position(input, labels,
                                                     [1, 2])
        assert_array_almost_equal(output1[0], output2)
        assert_array_almost_equal(output1[1], output3)
        assert_array_almost_equal(output1[2], output4)
        assert_array_almost_equal(output1[3], output5)


def test_center_of_mass01():
    "center of mass 1"
    expected = [0.0, 0.0]
    for type in types:
        input = np.array([[1, 0], [0, 0]], type)
        output = ndimage.center_of_mass(input)
        assert_array_almost_equal(output, expected)


def test_center_of_mass02():
    "center of mass 2"
    expected = [1, 0]
    for type in types:
        input = np.array([[0, 0], [1, 0]], type)
        output = ndimage.center_of_mass(input)
        assert_array_almost_equal(output, expected)


def test_center_of_mass03():
    "center of mass 3"
    expected = [0, 1]
    for type in types:
        input = np.array([[0, 1], [0, 0]], type)
        output = ndimage.center_of_mass(input)
        assert_array_almost_equal(output, expected)


def test_center_of_mass04():
    "center of mass 4"
    expected = [1, 1]
    for type in types:
        input = np.array([[0, 0], [0, 1]], type)
        output = ndimage.center_of_mass(input)
        assert_array_almost_equal(output, expected)


def test_center_of_mass05():
    "center of mass 5"
    expected = [0.5, 0.5]
    for type in types:
        input = np.array([[1, 1], [1, 1]], type)
        output = ndimage.center_of_mass(input)
        assert_array_almost_equal(output, expected)


def test_center_of_mass06():
    "center of mass 6"
    expected = [0.5, 0.5]
    input = np.array([[1, 2], [3, 1]], bool)
    output = ndimage.center_of_mass(input)
    assert_array_almost_equal(output, expected)


def test_center_of_mass07():
    "center of mass 7"
    labels = [1, 0]
    expected = [0.5, 0.0]
    input = np.array([[1, 2], [3, 1]], bool)
    output = ndimage.center_of_mass(input, labels)
    assert_array_almost_equal(output, expected)


def test_center_of_mass08():
    "center of mass 8"
    labels = [1, 2]
    expected = [0.5, 1.0]
    input = np.array([[5, 2], [3, 1]], bool)
    output = ndimage.center_of_mass(input, labels, 2)
    assert_array_almost_equal(output, expected)


def test_center_of_mass09():
    "center of mass 9"
    labels = [1, 2]
    expected = [(0.5, 0.0), (0.5, 1.0)]
    input = np.array([[1, 2], [1, 1]], bool)
    output = ndimage.center_of_mass(input, labels, [1, 2])
    assert_array_almost_equal(output, expected)


def test_histogram01():
    "histogram 1"
    expected = np.ones(10)
    input = np.arange(10)
    output = ndimage.histogram(input, 0, 10, 10)
    assert_array_almost_equal(output, expected)


def test_histogram02():
    "histogram 2"
    labels = [1, 1, 1, 1, 2, 2, 2, 2]
    expected = [0, 2, 0, 1, 1]
    input = np.array([1, 1, 3, 4, 3, 3, 3, 3])
    output = ndimage.histogram(input, 0, 4, 5, labels, 1)
    assert_array_almost_equal(output, expected)


def test_histogram03():
    "histogram 3"
    labels = [1, 0, 1, 1, 2, 2, 2, 2]
    expected1 = [0, 1, 0, 1, 1]
    expected2 = [0, 0, 0, 3, 0]
    input = np.array([1, 1, 3, 4, 3, 5, 3, 3])
    output = ndimage.histogram(input, 0, 4, 5, labels, (1,2))

    assert_array_almost_equal(output[0], expected1)
    assert_array_almost_equal(output[1], expected2)


def test_stat_funcs_2d():
    """Apply the stat funcs to a 2-d array."""
    a = np.array([[5,6,0,0,0], [8,9,0,0,0], [0,0,0,3,5]])
    lbl = np.array([[1,1,0,0,0], [1,1,0,0,0], [0,0,0,2,2]])

    mean = ndimage.mean(a, labels=lbl, index=[1, 2])
    assert_array_equal(mean, [7.0, 4.0])

    var = ndimage.variance(a, labels=lbl, index=[1, 2])
    assert_array_equal(var, [2.5, 1.0])

    std = ndimage.standard_deviation(a, labels=lbl, index=[1, 2])
    assert_array_almost_equal(std, np.sqrt([2.5, 1.0]))

    med = ndimage.median(a, labels=lbl, index=[1, 2])
    assert_array_equal(med, [7.0, 4.0])

    min = ndimage.minimum(a, labels=lbl, index=[1, 2])
    assert_array_equal(min, [5, 3])

    max = ndimage.maximum(a, labels=lbl, index=[1, 2])
    assert_array_equal(max, [9, 5])


if __name__ == "__main__":
    run_module_suite()
