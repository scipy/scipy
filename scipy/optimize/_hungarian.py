# Hungarian algorithm (Kuhn-Munkres) for solving the linear sum assignment
# problem. Taken from scikit-learn. Based on original code by Brian Clapper,
# adapted to NumPy by Gael Varoquaux.
# Further improvements by Ben Root, Vlad Niculae and Lars Buitinck.
# Rewritten by Jacob Turner
#
# Copyright (c) 2008 Brian M. Clapper <bmc@clapper.org>, Gael Varoquaux
# Author: Brian M. Clapper, Gael Varoquaux, Jacob Turner
# License: 3-clause BSD

import numpy as np
from collections import deque, namedtuple

Node = namedtuple('Node', 'previous next_row next_col')


def linear_sum_assignment(cost_matrix, maximize=False):

    """
    Algorithm that solves the minimum / maximum weighted assignment problem
    (or linear sum assignment problem) for bipartite graphs whose weighted
    bi-adjacency matrices are not necessarily square.

    Let C be an n x m matrix. Consider a set of pairs


    .. math::
        M = \\{(i, j) \\in [n] \\times [m] \\},

    or equivalently, edges in the bipartite graph.
    We say that row i is assigned to column j if (i, j) in M. M is a matching
    if every row is assigned to at most one column and likewise each column is
    assigned to at most one row.

    The problem is to find a matching M of size min(n, m) such that the sum of
    C[i,j] for (i, j) in M is minimized / maximized.

    Definitions:
        Many times we are performing operations on the graph whose edges are
        designated by zero entries in the assignment matrix and non-zero
        entries are non-edges. We call this graph the **0-induced graph**.

    If a row has been assigned a column, we say that it is **saturated**. In
    the algorithm, this is kept track of via the row_saturated and columns
    saturated vectors.

    Marking is done when we determine a minimum vertex cover of the 0-induced
    graph. We mark all vertices reachable by augmenting paths from free row
    vertices. From the resulting set the column vertices are kept and the
    complement of the row vertices are kept.

    In the algorithm description below, [3] refers to reference 3 below:
    Munkres, J. Algorithms for the Assignment and Transportation Problems.
    *J. SIAM*, 5(1):32-38, March, 1957.

    The algorithm consists of six main parts:
        1. A pre-processing of the input matrix for easy access in the main
           algorithms
        2. Adjust weights so that each row has a zero entry. (Step A in [3])
        3. The computation of a maximal matching, done greedily. If the maximal
           matching is maximum, we terminate (Step B in [3])
        4. Otherwise, we compute a minimum vertex cover of the bipartite graph
           defined by zeros in the weighted bi-adjacency matrix.
           (Step B in [3])
        5. Subtract smallest unmarked element. Subtract from each other
           unmarked element and add to every doubly marked element.  Return to
           step 2 above. (Step C in [3])
        6. Once a result is found, an assignment (i.e. a maximum matching) is
           found from the current known maximal matching by finding augmenting
           paths.
    N.B. in the above the description, we follow [3] at a high level. However,
    the particular algorithm for performing Step B in [3] (steps 3. and 4.
    above) differs from the procedure in the paper.

    Rephrasing Step B in [3] says that we must find a minimum vertex cover
    of the 0-induced graph. Our particular implementation of this is based on
    the proof of Konig's Theorem which provides an algorithm for finding
    a minimum vertex cover from a maximal matching (see e.g. reference [4]
    below).


    Parameters
     ----------
    cost_matrix : array
        The cost matrix of the bipartite graph.
    maximize : bool
        Specifies if the maximum weight matching should be computed. Default is
        False, meaning that minimum weight matching is computed

    Returns
    -------
    row_ind, col_ind : array
        An array of row indices and one of corresponding column indices giving
        the optimal assignment. The cost of the assignment can be computed
        as ``cost_matrix[row_ind, col_ind].sum()``. The row indices will be
        sorted; in the case of a square cost matrix they will be equal to
        ``numpy.arange(cost_matrix.shape[0])``.

    References
    ----------
    1. Harold W. Kuhn. The Hungarian Method for the assignment problem.
       *Naval Research Logistics Quarterly*, 2:83-97, 1955.
    2. Harold W. Kuhn. Variants of the Hungarian method for assignment
       problems. *Naval Research Logistics Quarterly*, 3: 253-258, 1956.
    3. Munkres, J. Algorithms for the Assignment and Transportation Problems.
       *J. SIAM*, 5(1):32-38, March, 1957.
    4. Cook, W. J., Cunningham, W. H., Pulleyblank, W. R., & Schrijver, A.
       (2009). Combinatorial optimization. Springer Oberwolfach Rep.,
       5(4), 2875-2942.
    5. https://en.wikipedia.org/wiki/Hungarian_algorithm

    """

    # Copy matrix to prevent side effects on the input cost matrix
    cost_matrix = np.asarray(cost_matrix).copy()

    # Validation

    if len(cost_matrix.shape) != 2:
        raise ValueError("expected a matrix (2-d array), got a %r array"
                         % (cost_matrix.shape,))

    if not (np.issubdtype(cost_matrix.dtype, np.number) or
            cost_matrix.dtype == np.dtype(np.bool)):
        raise ValueError("expected a matrix containing numerical entries,"
                         " got %s" % (cost_matrix.dtype,))

    if np.any(np.isinf(cost_matrix) | np.isnan(cost_matrix)):
        raise ValueError("matrix contains invalid numeric entries")

    if cost_matrix.dtype == np.dtype(np.bool):
        cost_matrix = cost_matrix.astype(np.int)

    # The internal algorithm assumes that there are at
    # least as many columns as rows
    if cost_matrix.shape[0] <= cost_matrix.shape[1]:
        return Munkres(cost_matrix,
                       transposed=False,
                       maximize=maximize).optimal_weight_matching()
    else:
        return Munkres(cost_matrix.transpose(),
                       transposed=True,
                       maximize=maximize).optimal_weight_matching()


class Munkres(object):
    """
    Class for finding maximum weight matchings and minimum vertex covers
    in bipartite graph
    """

    def __init__(self, matrix, transposed=False, maximize=False):
        if maximize:
            matrix = -matrix
        self.shape = matrix.shape
        self.marked = np.zeros(self.shape, dtype=bool)
        self.row_saturated = np.zeros(self.shape[0], dtype=bool)
        self.col_saturated = np.zeros(self.shape[1], dtype=bool)
        self.row_marked = np.zeros(self.shape[0], dtype=bool)
        self.col_marked = np.zeros(self.shape[1], dtype=bool)
        self.matrix = matrix
        self.transposed = transposed
        self.col_saturated_size = 0

    def _maximal_matching(self):
        """Find a maximal matching greedily"""

        # For each row, find the smallest element in that row and
        # subtract it from each element in its row.
        self.matrix -= self.matrix.min(axis=1)[:, np.newaxis]

        # Iterating through each zero-valued matrix entry, if neither the row
        # or column of the entry has been marked, add entry to the matching
        # and mark the its row and column as being assigned.
        for row, col in zip(*np.nonzero(self.matrix == 0)):
            if not self.row_saturated[row] and not self.col_saturated[col]:
                self.marked[row, col] = self.row_saturated[row]\
                    = self.col_saturated[col] = True
        self.col_saturated_size = self.col_saturated.sum()

    def _remove_covers(self):
        self.row_marked *= False
        self.col_marked *= False
        self.row_saturated *= False
        self.col_saturated *= False
        self.marked *= False
        self.col_saturated_size = 0

    def _min_vertex_cover(self):
        """Find a minimum vertex cover of 0-induced graph"""

        # Start with all unsaturated row vertices
        self.row_marked = self.row_saturated == False

        # We keep trying to find new vertices reachable by augmenting paths.
        while True:
            found = self.col_marked.sum()
            # Saturated column neighbors of rows from previous round
            self.col_marked = np.any(self.matrix[self.row_marked] == 0, axis=0)
            if found == self.col_marked.sum():
                self.col_saturated_size = found
                return
            else:
                # Mark rows that are matched with columns found above
                self.row_marked[np.any(self.marked[:, self.col_marked],
                                       axis=1)] = True

    def _aug_paths(self):
        """Find an augmenting path if one exists from maximal matching."""
        # Check unsaturated row vertices for augmenting paths
        for row in np.nonzero(self.row_saturated == False)[0]:
            path_row, path_col = self._aug_path(row)
            if not path_col:
                continue
            if not len(path_row + path_col) % 2:
                for i in range(len(path_row) - 1):
                    self.marked[path_row[i], path_col[i]] = 1
                    self.marked[path_row[i], path_col[i + 1]] = 0
                self.marked[path_row[-1], path_col[-1]] = 1
                self.row_saturated[path_row[-1]] = \
                    self.col_saturated[path_col[0]] = True

    def _aug_path(self, row):
        """"
        Recursively search for augmenting paths starting at row vertex 'row' in
        the 0-induced graph.
        """

        queue = deque()
        queue.append(Node(previous=None, next_col=None, next_row=row))
        visited_columns = set([])
        last_col = None
        path_row = []
        path_col = []
        found = False

        # We proceed to search for an augmented path
        # via a breadth-first search approach
        while queue and not found:
            current_node = queue.popleft()

            # We now check every column to see if we can extend
            # tentative augmented path with column vertex. We iterate
            # over column vertex neighbors of the 0-induced graph
            for col in np.nonzero(self.matrix[current_node.next_row] == 0)[0]:

                # We do not check vertices already on a tentative path
                if col in visited_columns:
                    continue

                # If col vertex is saturated, it we check to see if
                # it can extend tentative path
                if self.col_saturated[col]:

                    # We find row vertex it is matched with
                    row_index = np.argmax(self.marked[:, col])
                    visited_columns.add(col)

                    # We extend the tentative augmented path via linked list
                    queue.append(Node(previous=current_node,
                                      next_row=row_index,
                                      next_col=col))

                else:

                    # We have found the end of an augmented path.
                    # We add column vertex.
                    last_col = Node(previous=current_node,
                                    next_row=None,
                                    next_col=col)
                    found = True
                    break

        # Recreate augmented path from linked list if possible
        if last_col is not None:

            # add final column vertex
            path_col.append(last_col.next_col)
            last_col = last_col.previous

            # iterate back through linked list
            while last_col.previous is not None:
                path_row.append(last_col.next_row)
                path_col.append(last_col.next_col)
                last_col = last_col.previous

            # add initial row vertex
            path_row.append(last_col.next_row)

        return path_row, path_col

    def optimal_weight_matching(self):
        """Main algorithm. Runs the Hungarian Algorithm."""

        while True:
            # Subtract the smallest element in each row from every other
            # entry in same row and compute maximal matching of resulting
            # 0-inducted graph. (Steps 1 and 2 in Wikipedia)
            self._maximal_matching()

            # Find a minimum vertex cover of the 0-induced graph
            # (Step 3 in Wikipedia)
            self._min_vertex_cover()

            # If all rows are saturated, find the maximum matching and stop
            if self.row_marked.sum() == self.col_saturated_size:
                self._aug_paths()
                break

            # Find minimum value of uncovered edge weights
            minval = np.min(self.matrix[self.row_marked][:, self.col_marked !=
                                                         True])

            # Adjust the matrix weights according to Hungarian algorithm
            # (step 4 Wikipedia)
            self.matrix[self.row_marked] -= minval
            self.matrix[:, self.col_marked] += minval

            # Reset process and run again
            self._remove_covers()

        # If we transposed the input so that it was wide, undo that now
        # before returning answer.
        if self.transposed:
            return np.nonzero(self.marked.transpose() == 1)
        else:
            return np.nonzero(self.marked == 1)
