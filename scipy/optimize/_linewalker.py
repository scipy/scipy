"""Module providing a _linewalker function for numpy.minimize_scalar()."""

# To avoid the compilation error when referencing classes that appear later
# in the code, you can use forward references.
from __future__ import annotations
import copy
from dataclasses import dataclass, field
import math
from typing import Any, Callable, Dict, Tuple
import numpy as np
from ._cython_linewalker import CythonLinewalker

@dataclass
class _LinewalkerProgress:
    """ a class to store the linewalker algorithm's progess over the course of major iterations """
    # Function value at the grid index at ix_arg_min_evaluated
    f_min_evaluated: float = 1e15

    # Grid index of the minimum evaluated point
    ix_arg_min_evaluated: int = -1

    # Index of the sample when the ix_arg_min_evaluated was evaluated.
    i_sample_when_minimizer_was_evaluated: int = -1

    # Major iteration (of main "while" loop) when the function was evaluated at ix_arg_min_evaluated
    iter_when_minimizer_was_evaluated: int = 0

    # a positive integer denoting the maximum number of function evaluations (the "budget") that
    # can be made over the entire algorithm
    max_num_function_evaluations: int = 30

    # a positive integer denoting the maximum number of function evaluations to make per major
    # iteration
    max_num_evaluations_per_iteration: int = 1

    # a nonnegative integer denoting the number of major iterations (times entering the main WHILE
    # loop)
    num_major_iterations: int = 0
    previous_evaluated_objval_improvement: float = 0

    # True if we need to save certain metrics in each major iteration; False otherwise
    save_metrics_per_iter: bool = False

    # ixNewNonTabuInIter
    diversify_in_iter: np.ndarray = field(default_factory=lambda: np.array([], dtype=bool))
    num_local_extrema_in_iter: np.ndarray = field(default_factory=lambda: np.array([], dtype=int))
    tabu_tenure_in_iter: np.ndarray = field(default_factory=lambda: np.array([], dtype=int))

    def __init__(self, option: Dict['str', Any]):
        """ Initialize a LinewalkerProgess object. """
        for key, value in option.items():
            if hasattr(self, key):
                setattr(self, key, value)

        # Initialize arrays if saving metrics in each major iteration
        if self.save_metrics_per_iter:
            n: int = self.max_num_function_evaluations - option['initial_number_of_samples']
            self.diversify_in_iter = np.zeros(n, dtype=bool)
            self.num_local_extrema_in_iter = np.zeros(n, dtype=int)
            self.tabu_tenure_in_iter = np.zeros(n, dtype=int)

    def update_minimizer_evaluated(self,
        func_val: float, grid_index: int, i_sample_counter: int) -> None:
        """ Update the best minimizer evaluated and associated information. """
        if func_val >= self.f_min_evaluated:
            return
        self.previous_evaluated_objval_improvement = self.f_min_evaluated - func_val
        self.f_min_evaluated = func_val
        self.ix_arg_min_evaluated = grid_index
        self.i_sample_when_minimizer_was_evaluated = i_sample_counter
        self.iter_when_minimizer_was_evaluated = self.num_major_iterations

    def save_metrics(self,
            diversify: bool, surrogate: _Surrogate, tabu_struct: _TabuStruct) -> None:
        """ Save certain algorithm metrics if desired. """
        if self.save_metrics_per_iter:
            self.diversify_in_iter[self.num_major_iterations] = diversify
            self.num_local_extrema_in_iter[self.num_major_iterations] = \
                len(surrogate.ix_local_maxima) + len(surrogate.ix_local_minima)
            self.tabu_tenure_in_iter[self.num_major_iterations] = \
                tabu_struct.tabu_tenure_num_iterations


@dataclass
class _Surrogate:
    """ a class to store information related to a discrete surrogate. """
    # len(fit)
    grid_size: int = 1000

    # An array of floats (1D array of N elements)
    fit: np.ndarray = field(default_factory=lambda: np.array([], dtype=float))

    # Maximum value over all elements of fit
    max_fit: float = 0

    # Minimum value over all elements of fit
    min_fit: float = 0

    # max_fit-min_fit
    range_fit: float = 0

    # Index of the maximum (predicted) value in fit
    ix_arg_max_fit: int = 0

    # Index of the minimum (predicted) value in fit
    ix_arg_min_fit: int = 0

    # An array of floats (1D array of N elements)
    func: np.ndarray = field(default_factory=lambda: np.array([], dtype=float))

    # An array of distinct nonnegative integers indicating the indices at which
    # function evaluations have been made
    ix: np.ndarray = field(default_factory=lambda: np.array([], dtype=int))

    # ix sorted in ascending order
    ix_sorted: np.ndarray = field(default_factory=lambda: np.array([], dtype=int))

    # Function value at the grid index at ix_arg_min_evaluated
    f_min_evaluated: float = 1e15

    # Index of the minimum evaluated function value\
    ix_arg_min_evaluated: int = None

    # Indices of local maxima
    ix_local_maxima: np.ndarray = field(default_factory=lambda: np.array([], dtype=int))

    # Indices of local minima
    ix_local_minima: np.ndarray = field(default_factory=lambda: np.array([], dtype=int))

    # Indices of new local extrema: (a.k.a. ixNewLocalExtrema =)
    # ix_new = np.union1d(ix_new_local_minima, ix_new_local_maxima)
    ix_new: np.ndarray = field(default_factory=lambda: np.array([], dtype=int))

    # Indices of new local maxima: ix_new_local_maxima = np.setdiff1d(ix_local_maxima, ix)
    ix_new_local_maxima: np.ndarray = field(default_factory=lambda: np.array([], dtype=int))

    # Indices of new local minima: ix_new_local_minima = np.setdiff1d(ix_local_minima, ix)
    ix_new_local_minima: np.ndarray = field(default_factory=lambda: np.array([], dtype=int))

    # boolean array: sp[i]=True implies a function evaluation (sample) will be made at
    # grid point/index i. sp is short for "sample" or "sample point"
    sp: np.ndarray = field(default_factory=lambda: np.array([], dtype=bool))

    # An array of floats containing the x-coordinate associated with each grid index # Formerly Xset
    x_coord: np.ndarray = field(default_factory=lambda: np.array([], dtype=float))

    # First-derivative smoothing parameter
    alpha: float = 0.00

    # Second-derivative smoothing parameter
    mu: float = 0.01

    # a real number >= 0
    max_relative_tolerance: float = 1e-8

    # a real number >= 0
    min_relative_tolerance: float = 1e-8

    a: np.ndarray = field(default_factory=lambda: np.array([], dtype=float))

    # pentadiagonal linear solve parameters

    # Lowest (Leftmost) diagonal of pentadiagonal matrix from compute_fit_vectors5
    d1: np.ndarray = field(default_factory=lambda: np.array([], dtype=float))

    # Lower (second leftmost) diagonal of pentadiagonal matrix from compute_fit_vectors5
    d2: np.ndarray = field(default_factory=lambda: np.array([], dtype=float))

    # Main diagonal of pentadiagonal matrix from compute_fit_vectors5
    d3: np.ndarray = field(default_factory=lambda: np.array([], dtype=float))

    # Higher (second rightmost) diagonal of pentadiagonal matrix from compute_fit_vectors5
    d4: np.ndarray = field(default_factory=lambda: np.array([], dtype=float))

    # Highest (Rightmost) diagonal of pentadiagonal matrix from compute_fit_vectors5
    d5: np.ndarray = field(default_factory=lambda: np.array([], dtype=float))

    # Functions for Acceleration etc.
    penta = None
    compute_fit_vectors5 = None

    def __init__(self, x1: int, x2: int, option: Dict['str', Any]):
        """ Initialize a Surrogate object. """
        for key, value in option.items():
            if hasattr(self, key):
                setattr(self, key, value)

        self.fit = np.zeros(self.grid_size) # Initialize fit (surrogate / function approximation)

        # Define the smoothing matrix or matrices
        self.d1, self.d2, self.d3, self.d4, self.d5 = _compute_fit_vectors5(
            self.grid_size, self.alpha, self.mu)

        # Define set of candidate sample points
        self.x_coord = np.linspace(x1, x2, self.grid_size)

        # Select grid indices of initial samples

        # sp is short for "sample" or "sample point"
        self.select_initial_grid_indices_to_sample(option)

        # vector of indices where samples (function evaluations) have been taken
        # (or soon will be taken)
        self.ix = np.nonzero(self.sp)[0]

        self.ix_sorted = np.sort(self.ix)

    def insert(self, i: int, fnc_val: float) -> None:
        """ Insert the index i and the function value fnc_val into appropriate surrogate members """
        self.func[i] = fnc_val
        self.ix = np.append(self.ix, i)

        # Find the index where the value should be inserted to maintain sorted order
        index = np.searchsorted(self.ix_sorted, i)

        # Insert the value at the found index
        self.ix_sorted = np.insert(self.ix_sorted, index, i)

        self.sp[i] = 1
        if fnc_val < self.f_min_evaluated:
            self.ix_arg_min_evaluated = i
            self.f_min_evaluated = fnc_val

    def update_fit(self) -> None:
        """ Update the fit (the discrete approximation of the true underlying function,
        which we only observe via samples). """
        self.fit = self.penta(
            self.fit.size, self.d1, self.d2,
            self.d3+self.sp, self.d4, self.d5, self.func)


    def find_largest_unexplored_interval(self) -> np.ndarray:
        """Find the largest unexplored interval and bisect it (assumes no prior)

        Returns:
            np.ndarray: a numpy array containing a single element - an int corresponding
            to the grid index at which the largest unexplored interval is bisected
        """
        num_unexplored_points_per_interval_vector = self.ix_sorted[1:] - self.ix_sorted[:-1]
        idx_interval_vector = np.where(
            num_unexplored_points_per_interval_vector >=
            np.max(num_unexplored_points_per_interval_vector) - 1)[0]
        i_min_func_tiebreaker_left_end_point = self.ix_sorted[idx_interval_vector[0]]
        i_min_func_tiebreaker_right_end_point = self.ix_sorted[idx_interval_vector[0] + 1]

        # If there is more than one "max" interval with the same number of unexplored/unsampled
        # points, we break the tie by selecting the interval containing the minimum predicted
        # fit value
        if len(idx_interval_vector) > 1:
            min_func_tiebreaker = 1e15
            for i in idx_interval_vector:
                i_interval_left_end_point = self.ix_sorted[i]
                i_interval_right_end_point = self.ix_sorted[i + 1]
                tmp_min_func = np.min(
                    self.fit[i_interval_left_end_point:i_interval_right_end_point + 1])
                if tmp_min_func < min_func_tiebreaker:
                    min_func_tiebreaker = tmp_min_func
                    i_min_func_tiebreaker_left_end_point = i_interval_left_end_point
                    i_min_func_tiebreaker_right_end_point = i_interval_right_end_point

        ix_new = i_min_func_tiebreaker_left_end_point + \
            math.floor(0.5 * (
                i_min_func_tiebreaker_right_end_point -
                i_min_func_tiebreaker_left_end_point))
        return np.array([ix_new])

    def find_local_extrema(self) -> None:
        """
        Given a surrogate with a 'fit' array of N floats and an int array 'ix' indicating,
        1) find all strictly local (and therefore global) extrema (maxima and minima)
        values and associated indices of fit
        2) find all new strictly local extrema, i.e., indices of local extrema that
        are not already in 'ix'
        """
        fit = self.fit
        previous_range_fit = self.range_fit

        self.ix_arg_max_fit = np.argmax(fit)
        self.max_fit = fit[self.ix_arg_max_fit]
        self.ix_arg_min_fit = np.argmin(fit)
        self.min_fit = fit[self.ix_arg_min_fit]
        self.range_fit = self.max_fit - self.min_fit
        max_absolute_tolerance = self.max_relative_tolerance * previous_range_fit
        min_absolute_tolerance = self.min_relative_tolerance * previous_range_fit

        is_fit_local_maxima = fit[1:-1] > np.maximum(fit[:-2], fit[2:]) + max_absolute_tolerance
        is_fit_local_minima = fit[1:-1] < np.minimum(fit[:-2], fit[2:]) - min_absolute_tolerance

        ix_local_maxima = np.where(is_fit_local_maxima)[0] + 1
        ix_local_minima = np.where(is_fit_local_minima)[0] + 1

        self.ix_local_maxima = np.union1d(ix_local_maxima, self.ix_arg_max_fit)
        self.ix_local_minima = np.union1d(ix_local_minima, self.ix_arg_min_fit)

        self.ix_new_local_maxima = np.setdiff1d(self.ix_local_maxima, self.ix)
        self.ix_new_local_minima = np.setdiff1d(self.ix_local_minima, self.ix)
        self.ix_new = np.union1d(self.ix_new_local_minima, self.ix_new_local_maxima)

    def select_initial_grid_indices_to_sample(self, option: Dict['str', Any]) -> None:
        """
        Select the grid indices that should be sampled at initialization.

        Args:
            option (dict): _description_

        Raises:
            ValueError: If option['initial_sampling_method'] does not have a
            defined/implemented sampling strategy

        Returns:
            bool array: a boolean array sp such that sp[i] = True if index
            i is (or will soon be) sampled; False otherwise
        """
        # binary vector: sp[i]=True implies a function evaluation (sample) will
        # be made at grid point i
        self.sp = np.zeros(self.grid_size, dtype=bool)

        if option['initial_sampling_method'] == 'uniform':
            i_sp = np.round(
                np.linspace(0, self.grid_size-1, option['initial_number_of_samples'])).astype(int)
            for i in i_sp:
                self.sp[i] = True
        elif option['initial_sampling_method'] == 'random':
            # Randomly select k unique positions
            indices = np.random.choice(
                self.grid_size, option['initial_number_of_samples'], replace=False)
            self.sp[indices] = 1
        else:
            raise ValueError(
                'Invalid initial sampling method: ' + option['initial_sampling_method'])

    def determine_if_potential_minimum_is_undersampled(
        self, ix_candidate: int,
        num_immediate_neighbors: int,
        option: Dict[str, Any]) -> bool:
        """ Determine if a potential global minimum is undersampled.
            In early iterations when the fit could be poor/inaccurate, we may not want
            to sample near a potential minimum *too* frequently.
            In later iterations when the fit becomes more accurate, we may want to
            resample near a potential minimum if it does not have *too many* immediate neighbors
        """
        fraction_of_lower_range = option['min_fraction_of_range_improvement']
        maxnum_immediate_neighbors = option['maxnum_immediate_neighbors']
        if len(self.ix) > option['num_samples_before_increasing_maxnum_immediate_neighbors']:
            fraction_of_lower_range = min(
                1.0, option['range_improvement_multiplier'] * fraction_of_lower_range)
            maxnum_immediate_neighbors += 1
        fit_has_potential_better_global_minimum = \
            self.fit[ix_candidate] <= \
            self.func[self.ix_arg_min_evaluated] + fraction_of_lower_range * self.range_fit
        potential_minimum_is_undersampled = \
            fit_has_potential_better_global_minimum and \
            (num_immediate_neighbors <= maxnum_immediate_neighbors)
        return potential_minimum_is_undersampled

    def find_extremum_with_largest_normalized_deviation(self) -> int:
        """
        Find extremum whose interval has the largest normalized deviation from a straight line

        Returns:
            int: Index of the extremum with the largest normalized deviation
        """
        ix_new_largest_normalized_deviation = self.ix_new[0]
        max_normalized_mad = -1

        for ix_candidate in self.ix_new:
            ix_left_ept, ix_right_ept = \
                _find_left_and_right_neighbors(ix_candidate, self.ix_sorted)

            # Compute the straight-line equation (y=mx+b) through the left
            # and right endpoints associated with ix_candidate
            rise = self.func[ix_right_ept] - self.func[ix_left_ept]
            run = ix_right_ept - ix_left_ept
            slope = rise / run
            intercept = -slope * ix_left_ept + self.func[ix_left_ept]
            y_linear = [slope * x + intercept for x in range(ix_left_ept, ix_right_ept + 1)]

            # Compute normalized Mean Absolute Error (mae) between fit and the
            # straightline approximation
            num_unexplored_points = ix_right_ept - ix_left_ept
            mae = sum(
                [abs(y - self.fit[x]) \
                 for x, y in zip(range(ix_left_ept, ix_right_ept + 1), y_linear)])
            normalized_mae = mae / num_unexplored_points

            if normalized_mae > max_normalized_mad:
                max_normalized_mad = normalized_mae
                ix_new_largest_normalized_deviation = ix_candidate

        return ix_new_largest_normalized_deviation


@dataclass
class _TabuStruct:
    """
    A class to store information related tabu (forbidden) grid indices and
    asperation criteria. 
    """

    # An array of distinct nonnegative integers indicating the indices that are
    # currently deemed tabu
    index_list: np.ndarray = field(default_factory=lambda: np.array([], dtype=int))

    # An array of integers: iterAddedList[i] indicates the iteration when indexList[i]
    # was added to the tabu
    iter_added_list: np.ndarray = field(default_factory=lambda: np.array([], dtype=int))

    # Aspiration criterion 2 requires the previously evaluated solution to be
    # option.min_fraction_of_range_improvement*range_fit
    min_fraction_of_range_improvement: float = 0.01

    # a scalar between [0,1]. Usage: Given N grid indices and |ix| samples, we require
    # new sample indices to be at least (roughly) min_grid_points_separation_multiplier*(N/|ix|)
    # indices from an existing sample index. Example: min_grid_points_separation_multiplier=0.1,
    # N=1000, |ix|=17, then ceil(min_grid_points_separation_multiplier*N/|ix|)=6
    min_grid_points_separation_multiplier: float = 0.10

    # The larger the multiplier, the GREATER the number of grid points needed for separation
    max_grid_points_separation_multiplier: float = 0.25

    # Specified as a fraction of the total range of the function
    max_fractional_deviation_from_optimum: float = 0.01

    # grid_size / (2 * max_num_function_evaluations)
    tabu_grid_distance: float = 50

    # a non-negative integer; a value <= 0 turns off all tabu search-related features
    tabu_tenure_num_iterations: int = 5

    # a binary flag: 0 => static; 1 => dynamic
    tabu_tenure_is_dynamic: int = 1

    # a binary flag: 0 => not enforce; 1 => long-term tabu tenure is enforced
    enforce_long_term_tabu_tenure: int = 1

    # Nonnegative integer indicating the number of iterations remaining to override short-
    # and long-term tabu tenures. When the user's goal is to optimize, we may wish to ensure
    # that a sample (function evaluation) is made at a predicted minimizer. This parameter
    # is used in an asperation criterion.  a positive integer K implies that the algorithm
    # may override all tabu tenures when there are K or fewer major iterations remaining.
    # 0 => Linewalker will never force a sample by overriding a tabu tenure
    force_sample_at_predicted_minimizer: int = 0

    def __init__(self, surrogate: _Surrogate, option: Dict['str', Any]):
        """ Initialize a TabuStruct object. """
        # Update attributes based on the option dictionary
        for key, value in option.items():
            if hasattr(self, key):
                setattr(self, key, value)

        self.tabu_grid_distance = option['grid_size'] / (2 * option['max_num_function_evaluations'])

        if self.tabu_tenure_num_iterations > 0:
            self.index_list = copy.deepcopy(surrogate.ix)
            self.iter_added_list = np.zeros(surrogate.ix.size, dtype=int)

    def compute_grid_distance_to_all_samples(
        self, ix_candidate: int,
        surrogate: _Surrogate) -> Tuple[bool, int]:
        """
        Given a candidate index 'ix_candidate' of where to sample, compute the grid distance
        to all existing samples. The grid distance between index i and j is defined as abs(i-j),
        i.e., the number of indices between i and j.
        """
        abs_grid_distance_to_sample = np.abs(ix_candidate - surrogate.ix_sorted)
        if self.enforce_long_term_tabu_tenure == 1:
            # Ratio of (ix_candidate's objective function distance from the min or max of the
            # fit) to (half of the range of the fit). If ix_candidate's objval is near the
            # minimum or the maximum, then this fraction will be very small. If ix_candidate's
            # objval is near the average, then the fraction will be close to 1.
            distance_ratio_to_mid_fit = \
                min(surrogate.max_fit - surrogate.fit[ix_candidate], \
                    surrogate.fit[ix_candidate] - surrogate.min_fit) / \
                (0.5 * surrogate.range_fit)
            tmp_grid_points_separation_multiplier = \
                self.min_grid_points_separation_multiplier + \
                    distance_ratio_to_mid_fit * \
                        (self.max_grid_points_separation_multiplier - \
                         self.min_grid_points_separation_multiplier)
            min_num_grid_points_of_separation_required = \
                tmp_grid_points_separation_multiplier * \
                    surrogate.fit.size / surrogate.ix_sorted.size
            is_min_grid_distance_to_all_samples_satisfied = \
                np.min(abs_grid_distance_to_sample) > \
                min_num_grid_points_of_separation_required
        else:
            is_min_grid_distance_to_all_samples_satisfied = 1

        binary_vector = abs_grid_distance_to_sample <= self.tabu_grid_distance
        num_immediate_neighbors = np.sum(binary_vector)
        return is_min_grid_distance_to_all_samples_satisfied, num_immediate_neighbors

    def is_candidate_sufficiently_far_away(
        self,
        is_min_grid_distance_to_all_samples_satisfied: bool,
        ix_candidate: int) -> bool:
        """ Check if the candidate is sufficiently far away from all tabu points and
        previously-sampled points. """
        abs_grid_distance_to_tabu_sample = np.abs(ix_candidate - self.index_list)
        is_min_grid_distance_to_tabu_samples_satisfied = \
            np.min(abs_grid_distance_to_tabu_sample) > self.tabu_grid_distance
        is_min_grid_distance_condition_true = \
            is_min_grid_distance_to_tabu_samples_satisfied and \
            is_min_grid_distance_to_all_samples_satisfied
        return is_min_grid_distance_condition_true

    def update_short_term_tabu_tenure(self, num_local_extrema: int) -> None:
        """
        Increase the tabu tenure if the current surrogate has at least as many
        local extrema as the current tabu tenure;
        Decrease the tabu tenure if the current surrogate has strictly fewer
        than the number of local extrema;
        Otherwise, do not change the tabu tenure.

        Args:
            num_local_extrema (int): Number of local extrema in the current surrogate
        """
        if self.tabu_tenure_is_dynamic:
            if num_local_extrema >= self.tabu_tenure_num_iterations:
                self.tabu_tenure_num_iterations += 1
            elif num_local_extrema < \
                 self.tabu_tenure_num_iterations - 1 and self.tabu_tenure_num_iterations > 1:
                self.tabu_tenure_num_iterations -= 1

    def update_tabu_lists(self, i: int, _iter: int) -> None:
        """ Update """
        if self.tabu_tenure_num_iterations > 0: # Update tabu lists
            self.index_list = np.append(self.index_list, i)
            self.iter_added_list = np.append(self.iter_added_list, _iter)

    def remove_no_longer_tabu_indices(self, num_major_iterations: int) -> None:
        """
        Remove tabu indices that are no longer tabu (have exceeded their tabu tenure)
        """
        if len(self.index_list) > 0:
            indices_to_remove = \
                np.nonzero((num_major_iterations - self.iter_added_list) > \
                           self.tabu_tenure_num_iterations)[0]
            np.delete(self.index_list, indices_to_remove)
            np.delete(self.iter_added_list, indices_to_remove)

    def manage_tabu_struct(self,
        surrogate: _Surrogate,
        lw_progress: _LinewalkerProgress,
        # num_major_iterations: int,
        # previous_evaluated_objval_improvement: float,
        option: Dict['str', Any]) -> np.ndarray:
        """
        Manage tabu_struct
        """
        if self.tabu_tenure_num_iterations <= 0: # There is no tabu_struct to manage
            return surrogate.ix_new

        num_local_extrema = len(surrogate.ix_local_minima) + len(surrogate.ix_local_maxima)
        self.update_short_term_tabu_tenure(num_local_extrema)
        self.remove_no_longer_tabu_indices(lw_progress.num_major_iterations)

        # Identify which new candidate points are non-tabu
        if self.index_list.size == 0:
            ix_new_non_tabu = surrogate.ix_new
        else:
            # Initialize as an empty list, although this will later be converted to a numpy array
            ix_new_non_tabu = []
            for ix_candidate in surrogate.ix_new:
                ix_left_ept, ix_right_ept = \
                    _find_left_and_right_neighbors(ix_candidate, surrogate.ix_sorted)
                is_min_grid_distance_to_all_samples_satisfied, num_immediate_neighbors = \
                    self.compute_grid_distance_to_all_samples(ix_candidate, surrogate)

                # NonTabu Check 1: Is ix_candidate sufficiently far away from all tabu points
                # AND from all previously-sampled points
                is_min_grid_distance_condition_true = \
                    self.is_candidate_sufficiently_far_away(
                        is_min_grid_distance_to_all_samples_satisfied, ix_candidate)
                if is_min_grid_distance_condition_true:
                    ix_new_non_tabu.append(ix_candidate)
                    continue

                # NonTabu Check 2 - Aspiration criterion 1: Potential minimum is undersampled.
                # Overrides short- and long-term tabu tenures
                # Is ix_candidate a potential minimizer AND is there at most
                # maxnum_immediate_neighbors points in the interval |-----x------| that have
                # been sampled. If so, go back and sample again "around the corner"
                potential_minimum_is_undersampled = \
                    surrogate.determine_if_potential_minimum_is_undersampled(
                        ix_candidate, num_immediate_neighbors, option)
                if potential_minimum_is_undersampled:
                    ix_new_non_tabu.append(ix_candidate)
                    continue

                # NonTabu Check 3 - Aspiration criterion 2: Resample near a newly-found minimum.
                # Overrides short-term tabu tenures only
                # Is ix_candidate near a newly-found minimum that was just evaluated/found in
                # the previous iteration?
                # Note that fit(ix_arg_min_evaluated) may not equal min(fit), meaning that
                # Check 2 may not catch this possibility, so we include an additional catch
                new_minimum_was_just_evaluated = \
                    (ix_left_ept <= surrogate.ix_arg_min_evaluated <= ix_right_ept) \
                    and (lw_progress.previous_evaluated_objval_improvement >= \
                         self.min_fraction_of_range_improvement * surrogate.range_fit) \
                    and is_min_grid_distance_to_all_samples_satisfied
                if new_minimum_was_just_evaluated:
                    ix_new_non_tabu.append(ix_candidate)
                    continue

                # NonTabu Check 4 (Aspiration criterion 3: Force sample at predicted minimizer
                # force_sample_at_predicted_minimizer; overrides short- and long-term tabu tenures.
                # When the user's goal is to optimize, we may wish to ensure that a sample
                # (function evaluation) is made at a predicted minimizer
                if ((lw_progress.max_num_function_evaluations-len(surrogate.ix) <= \
                        self.force_sample_at_predicted_minimizer) and \
                    surrogate.fit[ix_candidate] == surrogate.min_fit):
                    ix_new_non_tabu.append(ix_candidate)

        return np.array(ix_new_non_tabu)


def _compute_fit_vectors5(n: int, alpha: float, mu: float):
    """
    N     ==> Number of grid (line segment) points
    alpha ==> first-derivative smoothing parameter
    mu    ==> second-derivative smoothing parameter

    Source: https://www.math.uakron.edu/~kreider/anpde/_penta.f

    RESULTS:  matrix has 5 bands, EADCF, with D being the main diagonal,
    E and a are the lower diagonals, and C and F are the upper diagonals.

    E is defined for rows i = 2:N-1, but is defined as E[0] to E[N-3]
    a is defined for rows i = 1:N-1, but is defined as a[0] to a[N-2]
    D is defined for rows i = 0:N-1
    C is defined for rows i = 0:N-2, but the last element isn't used
    F is defined for rows i = 0:N-3, but the last 2 elements aren't used
    """

    e = mu * np.ones((n-2), dtype=float)
    f = copy.deepcopy(e)

    a = np.ones((n-1), dtype=float) * -(alpha + 4 * mu)
    a[0]   = -(alpha + 2 * mu)
    a[n-2] = -(alpha + 2 * mu)
    c = copy.deepcopy(a)

    d = (2 * alpha + 6 * mu) * np.ones((n), dtype=float)
    d[0]   = mu
    d[1]   = 2 * alpha + 5 * mu
    d[n-2] = 2 * alpha + 5 * mu
    d[n-1] = mu

    return e, a, d, c, f

def _penta(n, e, _a, _d, _c, f, _b):

    """
    Vectors a,D,C,B are changed in the forward and backward loops below, so
    we need to "deepcopy" them to avoid errors in subsequent calls to _penta.
    Vectors E and F are never changed.
    """
    a = copy.deepcopy(_a)
    d = copy.deepcopy(_d)
    c = copy.deepcopy(_c)
    b = copy.deepcopy(_b)

    # Initialize solution vector X
    x = np.empty(n)

    # Forward elimination
    for i in range(1, n-1):
        xmult = a[i-1] / d[i-1]
        d[i] -= xmult * c[i-1]
        c[i] -= xmult * f[i-1]
        b[i] -= xmult * b[i-1]
        xmult = e[i-1] / d[i-1]
        a[i] -= xmult * c[i-1]
        d[i+1] -= xmult * f[i-1]
        b[i+1] -= xmult * b[i-1]

    # Backward substitution
    xmult = a[n-2] / d[n-2]
    d[n-1] -= xmult * c[n-2]
    x[n-1] = (b[n-1] - xmult * b[n-2]) / d[n-1]
    x[n-2] = (b[n-2] - c[n-2] * x[n-1]) / d[n-2]
    for i in range(n-3, -1, -1):
        x[i] = (b[i] - f[i] * x[i+2] - c[i] * x[i+1]) / d[i]

    return x

def _find_left_and_right_neighbors(ix_candidate: int, ix_sorted: np.ndarray) -> Tuple[int, int]:
    """
    Given an index ix_candidate and a sorted array ix_sorted of distinct int values,
    returns the elements ix_left_ept and ix_right_ept such that ix_sorted[iStar-1] =
    ix_left_ept < ix_candidate < ix_right_ept = ix_sorted[iStar] for some iStar in
    range(1,len(ix_sorted)) Because ix_sorted is sorted in ascending order and contains
    distinct int values, this can be accomplished trivially with binary search
    Example: ix_sorted=[1, 20, 30, 40, 50, 60, 70, 80, 90, 100], ix_candidate=71.
    find_left_and_right_neighbors(ix_candidate, ix_sorted) returns 60,70
    """
    assert ix_sorted[0] < ix_candidate and ix_candidate < ix_sorted[-1]
    i_right = np.searchsorted(ix_sorted, ix_candidate, side="right")
    i_left = i_right-1
    ix_left_ept = ix_sorted[i_left]
    ix_right_ept = ix_sorted[i_right]
    assert ix_left_ept != ix_right_ept
    return ix_left_ept, ix_right_ept

def _initialize_vector_of_function_evaluations(
    func: Callable[[float], float],
    lw_progress: _LinewalkerProgress,
    surrogate: _Surrogate,
    option: Dict['str', Any]) -> None:
    """
    Initialize function evaluations

    Args:
        lw_progress (_type_): _description_
        func (_type_): _description_
        surrogate (Surrogate): _description_
        option (_type_): _description_

    Returns:
        _type_: _description_
    """
    # an integer counter denoting the number of function evaluations (samples) made thus far
    i_sample_counter = 0

    # func stores the true function value at each grid index that is actually evaluated.
    # func would need to be modified if we allowed for noisy function evaluations because
    # we may choose to sample the same grid index multiple times.
    surrogate.func = np.zeros(option['grid_size'], dtype=float)

    # vector of indices where samples (function evaluations) have been taken (or soon will be taken)
    surrogate.ix = np.where(surrogate.sp)[0]

    for grid_index in surrogate.ix:
        i_sample_counter += 1
        # If the function has already been evaluated at the first grid point
        if grid_index == 0 and not np.isnan(option['f1']):
            surrogate.func[grid_index] = option.f1

        # If the function has already been evaluated at the last grid point
        elif grid_index == option['grid_size'] - 1 and not np.isnan(option['f2']):
            surrogate.func[grid_index] = option['f2']

        else:
            foo = func(surrogate.x_coord[grid_index])
            surrogate.func[grid_index] = foo

        if surrogate.func[grid_index] < lw_progress.f_min_evaluated:
            lw_progress.update_minimizer_evaluated(
                surrogate.func[grid_index], grid_index, i_sample_counter)
            surrogate.f_min_evaluated = surrogate.func[grid_index]
            surrogate.ix_arg_min_evaluated = grid_index

def _polish_non_tabu_candidates(
    surrogate: _Surrogate,
    diversify: bool,
    ix_new_to_evaluate: np.ndarray, # ixNewToEvaluate: Any,
    option: Dict[str, Any]) -> np.ndarray:
    """
    Polish non-tabu candidates
    """
    if option['max_fractional_deviation_from_optimum'] <= 0 or diversify:
        return ix_new_to_evaluate

    max_deviation_from_optimum = \
        option['max_fractional_deviation_from_optimum'] * surrogate.range_fit
    num_remaining_evaluations = option['max_num_function_evaluations'] - len(surrogate.ix_sorted)
    ix_new_to_evaluate_updated = np.zeros(len(ix_new_to_evaluate), dtype=int)

    for i, ix_candidate in enumerate(ix_new_to_evaluate):
        # force_sample_at_predicted_minimizer; overrides short- and long-term tabu tenures
        sample_at_predicted_minimizer: bool = \
            (num_remaining_evaluations <= option['force_sample_at_predicted_minimizer'] and \
            surrogate.fit[ix_candidate] == surrogate.min_fit)
        if sample_at_predicted_minimizer:
            ix_new_to_evaluate_updated[i] = ix_candidate
            continue

        ix_left_ept, ix_right_ept = \
            _find_left_and_right_neighbors(ix_candidate, surrogate.ix_sorted)
        ix_bisect_pt = ix_left_ept + math.ceil((ix_right_ept - ix_left_ept) / 2)

        if ix_right_ept - ix_candidate >= ix_candidate - ix_left_ept:
            i_temp = np.where(
                np.abs(surrogate.fit[ix_candidate:ix_bisect_pt+1] - \
                       surrogate.fit[ix_candidate]) <= \
                max_deviation_from_optimum)[0]
            ix_new_to_evaluate_updated[i] = ix_candidate + i_temp[-1] \
                                            if len(i_temp) > 0 \
                                            else ix_candidate
        else:
            i_temp = np.where(
                np.abs(surrogate.fit[ix_bisect_pt:ix_candidate+1] - \
                       surrogate.fit[ix_candidate]) <= \
                max_deviation_from_optimum)[0]
            ix_new_to_evaluate_updated[i] = ix_bisect_pt + i_temp[0] \
                                            if len(i_temp) > 0 \
                                            else ix_candidate

    return ix_new_to_evaluate_updated

def _sort_new_samples_to_evaluate(
    surrogate: _Surrogate,
    ix_new: np.ndarray,
    option: Dict['str', Any]) -> np.ndarray:
    """
    Sort candidate list of new indices to evaluate according to overall 'goal'
    """
    num_remaining_evaluations = \
        min(option['max_num_evaluations_per_iteration'],
            option['max_num_function_evaluations'] - len(surrogate.ix_sorted))
    ix_new_to_evaluate_unpolished = ix_new

    if len(ix_new) > num_remaining_evaluations:
        i_fit_val_sorted = np.argsort(surrogate.fit[ix_new])
        ix_new_to_evaluate_unpolished = ix_new[i_fit_val_sorted[:num_remaining_evaluations]]

    return ix_new_to_evaluate_unpolished

def _linewalker(func, brack, **kwargs):
    """
    Return the minimizer of a function of one variable using the linewalker
    method.

    Given a function of one variable and a bracketing interval,
    return a minimizer of the function or the minimum function value
    found with the function evalution budget max_num_function_evaluations.

    Parameters
    ----------
    func : callable func(x,*args)
        Objective function to minimize.
    brack : tuple
        A pair (xa, xb) such that xa < xb to be used as end points in constructing
        the search grid (see `scipy.optimize.linewalker`).
    kwargs : tuple, optional
        Additional arguments (if present), passed to func.
    grid_size: int
        Number of equally-spaced grid indices (candidate solutions) along the line segment of interest.
        Recommended range: 1000 <= grid_size <= 10000
    max_num_function_evaluations: int
        A positive integer indicating the maximum number of function evaluations that can be taken
        Recommended range: initial_number_of_samples+1 <= max_num_function_evaluations <= 50
    initial_number_of_samples: int
        A positive integer >= 2 denoting the number of initial samples (function evaluations) that must be made.
        Recommended range: initial_number_of_samples >= 10
    force_sample_at_predicted_minimizer: int
        An integer denoting the number of iterations in which to force a sample (function
        evaluation) to be taken at the predicted minimizer. An integer <= 0 means that no
        sample will be forced.  A positive integer K means that, if the grid index
        corresponding to the predicted minimizer has not been sampled, then a sample
        will be taken there in the last K major iterations.
        If K > 'max_num_function_evaluations'-'initial_number_of_samples', then every
        sample after the first 'initial_number_of_samples' initial samples are made will
        be at the index of a predicted global minimizer, assuming that this index has not
        already been sampled. The polishing (a.k.a. "sampling around the bend") is
        disabled if a global minimizer is sampled. It is not recommended to set this
        parameter > 2 since much can be learned from sampling a nonconvex function at
        diverse points.
        Recommended range: 1 <= force_sample_at_predicted_minimizer <= 2

    Returns
    -------
    surrogate : dict
        A dictionary containing various values from the surrogate
            'fun': float
                Minimum evaluated function value
            'x': float
                x coordinate of the minimum evaluated function value
            'nit': int
                Number of major iterations performed
            'nfev': int
                Number of function evaluations taken
            'success': bool
                True if minima was found
            'message': string
                Message
            'f_min_evaluated': float
                Minimum evaluated function value
            'minimizer_evaluated': float
                x coordinate of the minimum evaluated function value
            'num_major_iterations': int
                Number of major iterations performed
            'num_function_evaluations': int
                Number of function evaluations taken
            'f_min_predicted': float
                Minimum predicted function value
            'ix_arg_min_evaluated': int
                Grid index of the minimum evaluated function value
            'grid_size': int
                lenth of fit
            'fit': np.ndarray(float)
                An array of floats (1D array of N elements)
            'max_fit': float
                Maximum value over all elements of fit
            'min_fit': float
                Minimum value over all elements of fit
            'range_fit': float
                max_fit-min_fit
            'ix_arg_max_fit': int
                Index of the maximum (predicted) value in fit
            'ix_arg_min_fit': int
                Index of the minimum (predicted) value in fit
            'ix': np.ndarray(int)
                An array of distinct nonnegative integers indicating the indices at which
                function evaluations have been made
            'ix_sorted': np.ndarray(int)
                ix sorted in ascending order
            'x_coord': np.ndarray(float)
                An array of floats containing the x-coordinate associated with each grid index
            'sp': np.ndarray(bool)
                boolean array: sp[i]=True implies a function evaluation (sample) will be made at
                grid point/index i. sp is short for "sample" or "sample point"

    See also
    --------
    minimize_scalar: Interface to minimization algorithms for scalar
        univariate functions. See the 'Linewalker' `method` in particular.

    Notes
    -----
    Linewalker constructs a smooth surrogate on a set of equally-spaced grid points 
    by evaluating the true function at a sparse set of judiciously chosen grid points.
    At each iteration, the surrogate`s non-tabu local minima and maxima are identified
    as candidates for sampling. Tabu search constructs are also used to promote 
    diversification. If no non-tabu extrema are identified, a simple exploration step
    is taken by sampling the midpoint of the largest unexplored interval. The algorithm
    continues until a user-defined function evaluation limit is reached.

    Linewalker is particularly well-suited for nonconvex (multimodal) functions.
    If the underlying function is known to be convex, use a local solver.

    Examples
    --------
    We illustrate the behaviour of the function when `brack` is of
    size 2 and 3, respectively. In the case where `brack` is of the
    form (xa,xb), we can see for the given values, the output need
    not necessarily lie in the range ``(xa, xb)``.

    >>> import numpy as np
    >>> from scipy.optimize import linewalker, minimize_scalar
    >>> 
    >>> def f(x):
    >>>     return (x-1)**2
    >>> 
    >>> res = linewalker(f, brack=(-1, 2))
    >>> 
    >>> print(f"\
    >>>     {res['fun'] = },\n\
    >>>     {res['x'] = },\n\
    >>>     {res['nit'] = },\n\
    >>>     {res['nfev'] = },\n\
    >>>     {res['success'] = },\n\
    >>>     {res['message'] = }\n\
    >>>     {res['f_min_evaluated'] = },\n\
    >>>     {res['minimizer_evaluated'] = },\n\
    >>>     {res['num_major_iterations'] = },\n\
    >>>     {res['num_function_evaluations'] = },\n\
    >>>     {res['f_min_predicted'] = },\n\
    >>>     {res['ix_arg_min_evaluated'] = },\n\
    >>>     {res['grid_size'] = },\n\
    >>>     res['fit'] len, mean = {len(res['fit'])}{np.mean(res['fit'])},\n\
    >>>     {res['max_fit'] = },\n\
    >>>     {res['min_fit'] = },\n\
    >>>     {res['range_fit'] = },\n\
    >>>     {res['ix_arg_max_fit'] = },\n\
    >>>     {res['ix_arg_min_fit'] = },\n\
    >>>     res['ix'] len, mean = {len(res['ix'])}{np.mean(res['ix'])},\n\
    >>>     res['ix_sorted'] len, mean = {len(res['ix_sorted'])}{np.mean(res['ix_sorted'])},\n\
    >>>     res['x_coord'] len, mean = {len(res['x_coord'])}{np.mean(res['x_coord'])},\n\
    >>>     res['sp'] len, mean = {len(res['sp'])}{np.mean(res['sp'])}\
    >>> ")

    >>>     res['fun'] = np.float64(0.0),
    >>>     res['x'] = np.float64(1.0),
    >>>     res['nit'] = 19,
    >>>     res['nfev'] = 30,
    >>>     res['success'] = True,
    >>>     res['message'] = 'success'
    >>>     res['f_min_evaluated'] = np.float64(0.0),
    >>>     res['minimizer_evaluated'] = np.float64(1.0),
    >>>     res['num_major_iterations'] = 19,
    >>>     res['num_function_evaluations'] = 30,
    >>>     res['f_min_predicted'] = np.float64(-5.3577671966469964e-15),
    >>>     res['ix_arg_min_evaluated'] = np.int64(666),
    >>>     res['grid_size'] = 1000,
    >>>     res['fit'] len, mean = 10001.0015823438786493,
    >>>     res['max_fit'] = np.float64(3.9999999958493935),
    >>>     res['min_fit'] = np.float64(-5.3577671966469964e-15),
    >>>     res['range_fit'] = np.float64(3.999999995849399),
    >>>     res['ix_arg_max_fit'] = np.int64(0),
    >>>     res['ix_arg_min_fit'] = np.int64(666),
    >>>     res['ix'] len, mean = 30570.4666666666667,
    >>>     res['ix_sorted'] len, mean = 30570.4666666666667,
    >>>     res['x_coord'] len, mean = 10000.5,
    >>>     res['sp'] len, mean = 10000.03

    .. versionadded:: 1.15.2

    """

    default_args = {

        # Number of equally-spaced grid indices (candidate solutions) along the line segment of interest
        'grid_size': 1000,

        # First-derivative smoothing parameter for the underlying surrogate
        'alpha': 0.00,

        # Second-derivative smoothing parameter for the underlying surrogate
        'mu': 0.01,

        # A positive integer indicating the maximum number of function evaluations that can be taken
        'max_num_function_evaluations': 30,
        
        # A positive integer indicating the maximum number of function evaluations per major iteration
        'max_num_evaluations_per_iteration': 1,
        
        # A real number >= 0. The predicted value at a local maximum must be at least maximumRelativeTolerance greater than that of its nearest grid indices
        'max_relative_tolerance': 1e-8,
        
        # A real number >= 0. The predicted value at a local minimum must be at least minimumRelativeTolerance less than that of its nearest grid indices
        'min_relative_tolerance': 1e-8,
        
        # 'chaseGlobalOptima','improveOverallFit'
        'goal': 'chaseGlobalOptima',
        
        # Aspiration criterion 2 requires the previously evaluated solution to be min_fraction_of_range_improvement*range_fit, indicating that a significant objective value function improvement was obtained in the previous iteration and that this neighborhood is worth re-exploring
        'min_fraction_of_range_improvement': 0.01,
        
        'min_grid_points_separation_multiplier': 0.10,
        
        'max_grid_points_separation_multiplier': 0.25,
        
        'max_fractional_deviation_from_optimum': 0.01,

        # Used with potentialMinimumIsUndersampled: Recommended value = 2 when goal is
        # 'chaseGlobalOptima'; 1 when goal is 'improveOverallFit'
        'maxnum_immediate_neighbors': 1,
        
        # Used with potentialMinimumIsUndersampled
        'num_samples_before_increasing_maxnum_immediate_neighbors': 30,
        
        # Used with potentialMinimumIsUndersampled
        'range_improvement_multiplier': 10,
        
        # A non-negative integer; a value <= 0 turns off *all* (=short- and long-term) tabu search-related features
        'tabu_tenure_num_iterations': 5,
        
        # A binary flag: 0 ==> static; 1 ==> dynamic. If dynamic, the tabu tenure 'tabu_tenure_num_iterations' changes from major iteration to major iteration
        'tabu_tenure_is_dynamic': 1,
        
        # A binary flag: 0 ==> not enforce; 1 ==> long-term tabu tenure is enforced, i.e., we try to avoid sampling "too close" to an sampled index unless an aspiration criterion is met
        'enforce_long_term_tabu_tenure': 1,
        
        # An integer denoting the number of iterations in which to force a sample (function evaluation) to be taken at the predicted minimizer. An integer <= 0 means that no sample will be forced.  A positive integer K means that, if the grid index corresponding to the predicted minimizer has not been sampled, then a sample will be taken there in the last K major iterations.  If K > 'max_num_function_evaluations'-'initial_number_of_samples', then every sample after the first 'initial_number_of_samples' initial samples are made will be at the index of a predicted global minimizer, assuming that this index has not already been sampled. The polishing (a.k.a. "sampling around the bend") is disabled if a global minimizer is sampled. It is not recommended to set this parameter > 2 since much can be learned from sampling a nonconvex function at diverse points.
        'force_sample_at_predicted_minimizer': 1,
        
        'f1': float('nan'),
        
        'f2': float('nan'),
        
        # A positive integer >= 2 denoting the number of initial samples (function evaluations) that must be made
        'initial_number_of_samples': 10,
        
        'initial_sampling_method': 'uniform',
        
        'save_metrics_per_iter': 0, # False
    }

    default_args.update(kwargs)
    option = default_args

    if option['grid_size'] < 3:
        raise ValueError("grid_size must be greater than or equal to 3.")

    if option['initial_number_of_samples'] < 2:
        raise ValueError("initial_number_of_samples must be greater than or equal to 2.")

    if option['initial_number_of_samples'] > option['max_num_function_evaluations']:
        raise ValueError("initial_number_of_samples must be less than or equal to " \
                         "max_num_function_evaluations.")

    accelerated = option['acceleration'] if 'acceleration' in option else 'CYTHON_PENTA'

    lw_progress = _LinewalkerProgress(option)
    surrogate = _Surrogate(brack[0], brack[1], option)

    # Substitution for Acceleration
    if accelerated == 'NUMPY_ONLY_PENTA':
        surrogate.compute_fit_vectors5 = _compute_fit_vectors5
        surrogate.penta = _penta
    else: # Use CYTHON_PENTA
        adapter = CythonLinewalker(option['grid_size'])
        surrogate.compute_fit_vectors5 = adapter.compute_fit_vectors5
        surrogate.penta = adapter.compute_banded_matrix

    _initialize_vector_of_function_evaluations(func, lw_progress, surrogate, option)

    # ixIterEvaluated = np.zeros(surrogate.ix.size, dtype=int)
    # an integer counter denoting the number of function evaluations (samples) made thus far
    i_sample_counter = surrogate.ix.size

    tabu_struct = _TabuStruct(surrogate,option)

    # Main loop
    lw_progress.num_major_iterations = -1 # Will be immediately updated to 0
    while surrogate.ix.size < lw_progress.max_num_function_evaluations:
        lw_progress.num_major_iterations += 1
        surrogate.update_fit()
        surrogate.find_local_extrema()

        # Find non-tabu candidates
        ix_new_non_tabu = tabu_struct.manage_tabu_struct(surrogate, lw_progress, option)

        # Diversity/Explore
        if len(ix_new_non_tabu)==0:
            # Find the largest unexplored interval and bisect it (assumes no prior)
            ix_new = surrogate.find_largest_unexplored_interval()
            diversify = True
        else:
            ix_new = ix_new_non_tabu
            diversify = False

        # Determine new samples (grid indices) to evaluate
        ix_new_to_evaluate_unpolished = _sort_new_samples_to_evaluate(surrogate, ix_new, option)
        ix_new_to_evaluate_updated = \
            _polish_non_tabu_candidates(surrogate, diversify, ix_new_to_evaluate_unpolished, option)
        ix_new_to_evaluate = ix_new_to_evaluate_updated

        lw_progress.save_metrics(diversify, surrogate, tabu_struct) # Optional

        # Function evaluation Loop: Evaluate each new sample point
        lw_progress.previous_evaluated_objval_improvement = 0
        for i in ix_new_to_evaluate: # Parallelizable
            i_sample_counter += 1
            tmp_fnc_val = func(surrogate.x_coord[i])
            surrogate.insert(i,tmp_fnc_val)
            # ixIterEvaluated = np.append(ixIterEvaluated,num_major_iterations)
            lw_progress.update_minimizer_evaluated(tmp_fnc_val, i, i_sample_counter)
            tabu_struct.update_tabu_lists(i, lw_progress.num_major_iterations)

    surrogate.update_fit()
    surrogate.find_local_extrema()

    # Collect results in a dictionary
    return {
        'fun': lw_progress.f_min_evaluated,
        'x': surrogate.x_coord[surrogate.ix_arg_min_evaluated],
        'nit': lw_progress.num_major_iterations,
        'nfev': surrogate.ix.size,
        'success': True,
        'message': 'success',
        'f_min_evaluated': lw_progress.f_min_evaluated,
        'minimizer_evaluated': surrogate.x_coord[surrogate.ix_arg_min_evaluated],
        'num_major_iterations': lw_progress.num_major_iterations,
        'num_function_evaluations': surrogate.ix.size,
        'f_min_predicted': surrogate.fit[surrogate.ix_arg_min_fit],
        'ix_arg_min_evaluated': lw_progress.ix_arg_min_evaluated,
        'grid_size': surrogate.grid_size,
        'fit': surrogate.fit,
        'max_fit': surrogate.max_fit,
        'min_fit': surrogate.min_fit,
        'range_fit': surrogate.range_fit,
        'ix_arg_max_fit': surrogate.ix_arg_max_fit,
        'ix_arg_min_fit': surrogate.ix_arg_min_fit,
        'ix': surrogate.ix,
        'ix_sorted': surrogate.ix_sorted,
        'x_coord': surrogate.x_coord,
        'sp': surrogate.sp
    }
