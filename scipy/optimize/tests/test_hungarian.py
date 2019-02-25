import pytest
from scipy.optimize import _hungarian
import numpy as np
from scipy.sparse.sputils import matrix
from numpy.testing import assert_array_equal
from pytest import raises as assert_raises

@pytest.fixture
def matrix_for_tests():
    return np.array([[5, 0, 2, 0], [1, 3, 4, 0], [2, 2, 0, 2]])


"""================= Old Tests ================="""


def test_linear_sum_assignment():
    for cost_matrix, expected_cost in [
        # Square
        ([[400, 150, 400],
          [400, 450, 600],
          [300, 225, 300]],
         [150, 400, 300]
         ),

        # Rectangular variant
        ([[400, 150, 400, 1],
          [400, 450, 600, 2],
          [300, 225, 300, 3]],
         [150, 2, 300]),

        # Square
        ([[10, 10, 8],
          [9, 8, 1],
          [9, 7, 4]],
         [10, 1, 7]),

        # Rectangular variant
        ([[10, 10, 8, 11],
          [9, 8, 1, 1],
          [9, 7, 4, 10]],
         [10, 1, 4]),

        # n == 2, m == 0 matrix
        ([[], []],
         []),
    ]:
        cost_matrix = np.array(cost_matrix)
        row_ind, col_ind = _hungarian.linear_sum_assignment(cost_matrix)
        assert_array_equal(row_ind, np.sort(row_ind))
        assert_array_equal(expected_cost, cost_matrix[row_ind, col_ind])

        cost_matrix = cost_matrix.T
        row_ind, col_ind = _hungarian.linear_sum_assignment(cost_matrix)
        assert_array_equal(row_ind, np.sort(row_ind))
        assert_array_equal(np.sort(expected_cost),
                           np.sort(cost_matrix[row_ind, col_ind]))


def test_linear_sum_assignment_input_validation():
    assert_raises(ValueError, _hungarian.linear_sum_assignment, [1, 2, 3])

    C = [[1, 2, 3], [4, 5, 6]]

    assert_array_equal(linear_sum_assignment(C),
                       linear_sum_assignment(np.asarray(C)))
    assert_array_equal(linear_sum_assignment(C),
                       linear_sum_assignment(matrix(C)))


    I = np.identity(3)
    assert_array_equal(_hungarian.linear_sum_assignment(I.astype(np.bool)),
                       _hungarian.linear_sum_assignment(I))
    assert_raises(ValueError, _hungarian.linear_sum_assignment, I.astype(str))

    I[0][0] = np.nan
    assert_raises(ValueError, _hungarian.linear_sum_assignment, I)

    I = np.identity(3)
    I[1][1] = np.inf
    assert_raises(ValueError, _hungarian.linear_sum_assignment, I)


"""================= Test Internal State of Algorithm ================="""


def test_munkres_init_matrix_negation_minimize_option(matrix_for_tests):
    """Test that on initialization, all entries of cost matrix are negated"""
    munkres = _hungarian.Munkres(-matrix_for_tests)
    matrix = np.array([[-5, 0, -2, 0], [-1, -3, -4, 0], [-2, -2, 0, -2]],
                      dtype=float)
    assert_array_equal(munkres.matrix, matrix)


def test_munkres_init_matrix_negation_maximize_option(matrix_for_tests):
    """Test that on initialization, all entries of cost matrix are negated"""
    munkres = _hungarian.Munkres(matrix_for_tests, maximize=True)
    matrix = np.array([[-5, 0, -2, 0], [-1, -3, -4, 0], [-2, -2, 0, -2]],
                      dtype=float)
    assert_array_equal(munkres.matrix, matrix)


def test_maximal_matching_matrix_adjustment_minimize_option(matrix_for_tests):
    """
    Test that _maximal_matching method correctly subtracts the smallest element
    of each row from every element in the same row
    """
    munkres = _hungarian.Munkres(-matrix_for_tests)
    munkres._maximal_matching()
    matrix = np.array([[0, 5, 3, 5], [3, 1, 0, 4], [0, 0, 2, 0]], dtype=float)
    assert_array_equal(munkres.matrix, matrix)


def test_maximal_matching_matrix_adjustment_maximize_option(matrix_for_tests):
    """
    Test that _maximal_matching method correctly subtracts the smallest element
    of each row from every element in the same row
    """
    munkres = _hungarian.Munkres(matrix_for_tests, maximize=True)
    munkres._maximal_matching()
    matrix = np.array([[0, 5, 3, 5], [3, 1, 0, 4], [0, 0, 2, 0]], dtype=float)
    assert_array_equal(munkres.matrix, matrix)


def test_row_col_saturated_maximal_matching_minimize_option(matrix_for_tests):
    """
    Test that the computation of a maximal matching for the 0-induced graph
    correctly labels the appropriate row and column vertices as saturated
    """
    munkres = _hungarian.Munkres(-matrix_for_tests)
    munkres._maximal_matching()
    assert_array_equal(munkres.row_saturated, np.array([True, True, True],
                                                       dtype=bool))
    assert_array_equal(munkres.col_saturated, np.array([True, True, True,
                                                        False], dtype=bool))


def test_row_col_saturated_maximal_matching_maximize_option(matrix_for_tests):
    """
    Test that the computation of a maximal matching for the 0-induced graph
    correctly labels the appropriate row and column vertices as saturated
    """
    munkres = _hungarian.Munkres(matrix_for_tests, maximize=True)
    munkres._maximal_matching()
    assert_array_equal(munkres.row_saturated, np.array([True, True, True],
                                                       dtype=bool))
    assert_array_equal(munkres.col_saturated, np.array([True, True, True,
                                                        False], dtype=bool))


def test_remove_covers(matrix_for_tests):
    """
    Test that the remove covers function resets all appropriate vectors to have
    all entries False and that marked contain only zeros
    """
    munkres = _hungarian.Munkres(matrix_for_tests)
    munkres.col_saturated += True
    munkres.row_saturated += True
    munkres.marked += True
    munkres.row_marked += True
    munkres.col_marked += True

    assert_array_equal(munkres.marked, np.ones((3, 4), dtype=bool))
    assert_array_equal(munkres.col_saturated, np.array([True, True, True,
                                                        True], dtype=bool))
    assert_array_equal(munkres.row_saturated, np.array([True, True, True],
                                                       dtype=bool))
    assert_array_equal(munkres.col_marked, np.array([True, True, True, True],
                                                    dtype=bool))
    assert_array_equal(munkres.row_marked, np.array([True, True, True],
                                                    dtype=bool))

    munkres._remove_covers()
    assert_array_equal(munkres.marked, np.zeros((3, 4), dtype=bool))
    assert_array_equal(munkres.col_saturated, np.array([False, False, False,
                                                        False], dtype=bool))
    assert_array_equal(munkres.row_saturated, np.array([False, False, False],
                                                       dtype=bool))
    assert_array_equal(munkres.col_marked, np.array([False, False, False,
                                                     False], dtype=bool))
    assert_array_equal(munkres.row_marked, np.array([False, False, False],
                                                    dtype=bool))


"""================ Test Correctness of output of Algorithm ================"""


def test_maximal_matching_marked_minimize_option(matrix_for_tests):
    """
    Test that the matrix encoding the entries of a maximal
    matching of the 0-induced matrix are computed correctly
    """
    # A cost matrix where the ones mark a minimum assignment
    # so should be unchanged by algorithm
    cost_matrix = np.array([[0, -1, 0], [0, 0, -1]])
    munkres = _hungarian.linear_sum_assignment(cost_matrix)
    row_result = np.array([0, 1])
    col_result = np.array([1, 2])
    assert_array_equal(munkres[0], row_result)
    assert_array_equal(munkres[1], col_result)


def test_maximal_matching_marked_maximize_option(matrix_for_tests):
    """
    Test that the matrix encoding the entries of a maximal
    matching of the 0-induced matrix are computed correctly
    """
    # A cost matrix where the ones mark a maximum assignment
    # so should be unchanged by algorithm
    cost_matrix = np.array([[0, 1, 0], [0, 0, 1]])
    munkres = _hungarian.linear_sum_assignment(cost_matrix, maximize=True)
    row_result = np.array([0, 1])
    col_result = np.array([1, 2])
    assert_array_equal(munkres[0], row_result)
    assert_array_equal(munkres[1], col_result)


def test_transposing():
    """
    Test that algorithm can handle both tall and wides matrices and does
    so correctly
    """
    # A cost matrix where the ones mark a maximum assignment
    # so should be unchanged by algorithm
    cost_matrix = np.array([[0, 1, 0], [0, 0, 1]])
    munkres_1 = _hungarian.linear_sum_assignment(cost_matrix,
                                                 maximize=True)
    munkres_2 = _hungarian.linear_sum_assignment(cost_matrix.transpose(),
                                                 maximize=True)
    row_result = np.array([0, 1])
    col_result = np.array([1, 2])
    assert_array_equal(munkres_1[0], row_result)
    assert_array_equal(munkres_1[1], col_result)
    assert_array_equal(munkres_2[0], col_result)
    assert_array_equal(munkres_2[1], row_result)


def test_aug_paths_1():
    """
    Test the algorithm that finds a maximum matching from a maximal matching
    via augmenting paths
    """
    # Original biadjacency matrix where zeros represent an edge and
    # non-zero values represent  non-edges
    munkres = _hungarian.linear_sum_assignment(np.array([[0, 0, 0, 0, 0],
                                                         [0, 0, 0, 0, 1],
                                                         [0, 0, 0, 1, 1],
                                                         [0, 0, 1, 1, 1],
                                                         [0, 1, 1, 1, 1]],
                                                        dtype=float))
    # A maximal matching that is not maximum that will be found by
    # _maximal_matching method
    #  [[1, 0, 0, 0, 0],
    #   [0, 1, 0, 0, 0],
    #   [0, 0, 1, 0, 0],
    #   [0, 0, 0, 0, 0],
    #   [0, 0, 0, 0, 0]]

    # The resulting (unique) maximum matching that should result after
    # calling the _aug_paths method
    marked = np.array([[0, 0, 0, 0, 1],
                       [0, 0, 0, 1, 0],
                       [0, 0, 1, 0, 0],
                       [0, 1, 0, 0, 0],
                       [1, 0, 0, 0, 0]], dtype=bool)
    assert_array_equal(munkres, np.nonzero(marked == 1))


def test_aug_paths_2():
    """
    Test the algorithm that finds a maximum matching from a maximal
    matching via augmenting paths
    """
    # Original biadjacency matrix where zeros represent an edge and non-zero
    # values represent non-edges
    munkres = _hungarian.linear_sum_assignment(np.array([[0, 1, 1, 1, 1, 1],
                                                         [0, 0, 0, 1, 1, 1],
                                                         [0, 1, 1, 0, 0, 1],
                                                         [1, 1, 0, 0, 1, 0],
                                                         [1, 1, 1, 0, 1, 1]],
                                                        dtype=float))
    # A maximal matching that is not maximum that will be found by
    # _maximal_matching method
    #  [[1, 0, 0, 0, 0, 0],
    #   [0, 1, 0, 0, 0, 0],
    #   [0, 0, 0, 1, 0, 0],
    #   [0, 0, 1, 0, 0, 0],
    #   [0, 0, 0, 0, 0, 0]]

    # The resulting maximum matching after _aug_paths method is called
    marked = np.array([[1, 0, 0, 0, 0, 0],
                       [0, 1, 0, 0, 0, 0],
                       [0, 0, 0, 0, 1, 0],
                       [0, 0, 1, 0, 0, 0],
                       [0, 0, 0, 1, 0, 0]], dtype=bool)

    assert_array_equal(munkres, np.nonzero(marked == 1))


def test_min_weight_matching_1(matrix_for_tests):
    """ Test that the correct minimum weight matching is found"""
    # Fully saturated case, wide
    munkres_one = _hungarian.linear_sum_assignment(-matrix_for_tests)
    row_result = np.array([0, 1, 2])
    col_result = np.array([0, 2, 1])
    assert_array_equal(munkres_one[0], row_result)
    assert_array_equal(munkres_one[1], col_result)


def test_min_weight_matching_2():
    """ Test that the correct minimum weight matching is found"""
    # Not saturated case, wide
    munkres_two = _hungarian.linear_sum_assignment(-np.array([[5, 0, 2, 0],
                                                              [1, 3, 4, 0],
                                                              [2, 0, 0, 0]],
                                                             dtype=float))
    row_result = np.array([0, 1, 2])
    col_result = np.array([0, 2, 1])
    assert_array_equal(munkres_two[0], row_result)
    assert_array_equal(munkres_two[1], col_result)


def test_min_weight_matching_3():
    """ Test that the correct minimum weight matching is found"""
    # Not saturated case, tall
    munkres_three = _hungarian.linear_sum_assignment(-np.array([[5, 0, 2, 0],
                                                                [5, 0, 2, 0],
                                                                [5, 0, 2, 0],
                                                                [1, 3, 4, 0],
                                                                [2, 2, 0, 2]],
                                                               dtype=float))
    row_result = np.array([0, 1, 3, 4])
    col_result = np.array([0, 2, 1, 3])
    assert_array_equal(munkres_three[0], row_result)
    assert_array_equal(munkres_three[1], col_result)


def test_min_weight_matching_4():
    """ Test that the correct minimum weight matching is found"""
    # Saturated case tall
    munkres_four = _hungarian.linear_sum_assignment(-np.array([[5, 0, 2, 0],
                                                               [5, 0, 2, 0],
                                                               [5, 0, 2, 0],
                                                               [1, 3, 4, 0],
                                                               [2, 2, 0, 2],
                                                               [2, 2, 0, 2]],
                                                              dtype=float))

    row_result = np.array([0, 3, 4, 5])
    col_result = np.array([0, 2, 3, 1])
    assert_array_equal(munkres_four[0], row_result)
    assert_array_equal(munkres_four[1], col_result)


def test_max_weight_matching_1(matrix_for_tests):
    """ Test that the correct maximum weight matching is found"""
    # Fully saturated case, wide
    munkres_one = _hungarian.linear_sum_assignment(matrix_for_tests, maximize=True)
    row_result = np.array([0, 1, 2])
    col_result = np.array([0, 2, 1])
    assert_array_equal(munkres_one[0], row_result)
    assert_array_equal(munkres_one[1], col_result)


def test_max_weight_matching_2():
    """ Test that the correct maximum weight matching is found"""
    # Not saturated case, wide
    munkres_two = _hungarian.linear_sum_assignment(np.array([[5, 0, 2, 0],
                                                             [1, 3, 4, 0],
                                                             [2, 0, 0, 0]],
                                                            dtype=float),
                                                   maximize=True)
    row_result = np.array([0, 1, 2])
    col_result = np.array([0, 2, 1])
    assert_array_equal(munkres_two[0], row_result)
    assert_array_equal(munkres_two[1], col_result)


def test_max_weight_matching_3():
    """ Test that the correct maximum weight matching is found"""
    # Not saturated case, tall
    munkres_three = _hungarian.linear_sum_assignment(np.array([[5, 0, 2, 0],
                                                               [5, 0, 2, 0],
                                                               [5, 0, 2, 0],
                                                               [1, 3, 4, 0],
                                                               [2, 2, 0, 2]],
                                                              dtype=float),
                                                     maximize=True)
    row_result = np.array([0, 1, 3, 4])
    col_result = np.array([0, 2, 1, 3])
    assert_array_equal(munkres_three[0], row_result)
    assert_array_equal(munkres_three[1], col_result)


def test_max_weight_matching_4():
    """ Test that the correct maximum weight matching is found"""
    # Saturated case tall
    munkres_four = _hungarian.linear_sum_assignment(np.array([[5, 0, 2, 0],
                                                              [5, 0, 2, 0],
                                                              [5, 0, 2, 0],
                                                              [1, 3, 4, 0],
                                                              [2, 2, 0, 2],
                                                              [2, 2, 0, 2]],
                                                             dtype=float),
                                                    maximize=True)

    row_result = np.array([0, 3, 4, 5])
    col_result = np.array([0, 2, 3, 1])
    assert_array_equal(munkres_four[0], row_result)
    assert_array_equal(munkres_four[1], col_result)



