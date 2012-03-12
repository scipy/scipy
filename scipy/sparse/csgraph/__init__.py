r"""
==============================================================
Compressed Sparse Graph Routines (:mod:`scipy.sparse.csgraph`)
==============================================================

.. currentmodule:: scipy.sparse.csgraph

Fast graph algorithms based on sparse matrix representations.

Contents
========

.. autosummary::
   :toctree: generated/

   connected_components -- determine connected components of a graph
   laplacian -- compute the laplacian of a graph
   shortest_path -- compute the shortest path between points on a positive graph
   dijkstra -- use Dijkstra's algorithm for shortest path
   floyd_warshall -- use the Floyd-Warshall algorithm for shortest path
   bellman_ford -- use the Bellman-Ford algorithm for shortest path
   johnson -- use Johnson's algorithm for shortest path
   breadth_first_order -- compute a breadth-first order of nodes
   depth_first_order -- compute a depth-first order of nodes
   breadth_first_tree -- construct the breadth-first tree from a given node
   depth_first_tree -- construct a depth-first tree from a given node
   minimum_spanning_tree -- construct the minimum spanning tree of a graph

Graph Representations
=====================
This module uses graphs which are stored in a matrix format.  A
graph with N nodes can be represented by an (N x N) adjacency matrix G.
If there is a connection from node i to node j, then G[i, j] = w, where
w is the weight of the connection.  For nodes i and j which are
not connected, G[i, j] = :math:`\infty` (for dense representations) or
is not defined (for sparse representations).

The representation of un-connected nodes (null edges) brings up a point of
potential confusion for converting between sparse and dense graph
representations.  Consider the following example::

    >>> from scipy.sparse import csr_matrix
    >>> X = np.array([[0, 1, 0],
    ...               [1, 0, 2],
    ...               [0, 2, 0]])
    >>> Xcsr = csr_matrix(X)

Here, ``X`` and ``Xcsr`` actually represent different graphs!  They
represent the following undirected graph configurations::

              X                  Xcsr

             (0)                  (0)
            /   \                    \ 
           0     1                    1 
          /       \                    \ 
        (2)---2---(1)        (2)---2---(1)
    
This is because the explicit zero-weights in the dense representation ``X``
are removed in the conversion to sparse format. For this reason, the csgraph
module includes utility routines to convert graphs to the canonical
representation where null edges are indicated by infinite weights, and to
efficiently switch between sparse and dense representations of graphs.
We'll convert our non-canonical graph to canonical dense or sparse form by
indicating that the null value is zero::

    >>> from scipy.sparse import csgraph
    >>> X_canonical = csgraph.canonical_from_dense(X, null_value=0)
    >>> print X_canonical
    [[ inf   1.  inf]
     [  1.  inf   2.]
     [ inf   2.  inf]]
    >>> Xcsr = csgraph.csgraph_from_dense(X, null_value=0)
    >>> print csgraph.csgraph_to_dense(Xcsr)
    [[ inf   1.  inf]
     [  1.  inf   2.]
     [ inf   2.  inf]]

The routines in this module accept as input either sparse (csr or csr) graph
representations, or dense representations in this canonical form.

Directed vs. Undirected
-----------------------
Matrices may represent either directed or undirected graphs.  This is
specified throughout the csgraph module by a boolean keyword.  Graphs are
assumed to be directed by default. In a directed graph, traversal from node
i to node j can be accomplished over the edge G[i, j], but not the edge
G[j, i].  In a non-directed graph, traversal from node i to node j can be
accomplished over either G[i, j] or G[j, i].  If both edges are not null,
and the two have unequal weights, then the smaller of the two is used.
Note that a symmetric matrix will represent an undirected graph, regardless
of whether the 'directed' keyword is set to True or False.
"""

__docformat__ = "restructuredtext en"

__all__ = ['cs_graph_components',
           'connected_components',
           'cs_graph_components',
           'laplacian',
           'shortest_path',
           'floyd_warshall',
           'dijkstra',
           'bellman_ford',
           'johnson',
           'breadth_first_order',
           'depth_first_order',
           'breadth_first_tree',
           'depth_first_tree',
           'minimum_spanning_tree',
           'construct_dist_matrix',
           'reconstruct_path',
           'csgraph_from_dense',
           'csgraph_to_dense',
           'canonical_from_dense',
           'NegativeCycleError']

from _components import cs_graph_components
from _laplacian import laplacian
from _shortest_path import shortest_path, floyd_warshall, dijkstra,\
    bellman_ford, johnson, NegativeCycleError
from _traversal import breadth_first_order, depth_first_order, \
    breadth_first_tree, depth_first_tree, connected_components
from _min_spanning_tree import minimum_spanning_tree
from _tools import construct_dist_matrix, reconstruct_path,\
    csgraph_from_dense, csgraph_to_dense, canonical_from_dense

from numpy import deprecate as _deprecate
cs_graph_components = _deprecate(connected_components,
                                 'cs_graph_components',
                                 'csgraph.connected_components')

from numpy.testing import Tester
test = Tester().test
