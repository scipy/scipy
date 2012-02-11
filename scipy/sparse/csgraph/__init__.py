"""
=====================================================
Compressed Sparse Graph (:mod:`scipy.sparse.csgraph`)
=====================================================

.. currentmodule:: scipy.sparse.csgraph

Fast graph algorithms based on sparse format.

Contents
========

.. autosummary::
   :toctree: generated/

   connected_components - determine connected components of a graph
   laplacian - compute the laplacian of a graph
   shortest_path - compute the shortest path between points on a positive graph
   dijkstra - use Dijkstra's algorithm for shortest path
   floyd_warshall - use the Floyd-Warshall algorithm for shortest path
   breadth_first_order - compute a breadth-first order of nodes
   depth_first_order - compute a depth-first order of nodes
   breadth_first_tree - construct the breadth-first tree from a given node
   depth_first_tree - construct a depth-first tree from a given node
   minimum_spanning_tree - construct the minimum spanning tree of a graph

Format
======
This module uses graphs which are stored in a compressed sparse format.  A
graph with N nodes can be represented by an (N x N) adjacency matrix G.
If there is a connection from node i to node j, then G[i, j] = w, where
w is the (nonzero) weight of the connection.  For nodes i and j which are
not connected, G[i, j] = 0.  Though this convention does not allow for the
relatively uncommon case of graphs with edges of weight zero, it has the
advantage that G can be stored efficiently using scipy's sparse matrix
formats.

Directed vs. Undirected
-----------------------
Matrices may represent either directed or undirected graphs.  This is generally
specified by a boolean keyword.  Unless otherwise specified, graphs are
assumed to be directed.  In a directed graph, traversal from node i to node j
can be accomplished over the edge G[i, j], but not the edge G[j, i].  In a
non-directed graph, traversal from node i to node j can be accomplished over
either G[i, j] or G[j, i].  If both edges are nonzero and have unequal weights,
then the smaller of the two is used.
"""

__docformat__ = "restructuredtext en"

__all__ = ['connected_components',
           'cs_graph_components',
           'laplacian',
           'shortest_path',
           'floyd_warshall',
           'dijkstra',
           'breadth_first_order',
           'depth_first_order',
           'breadth_first_tree',
           'depth_first_tree',
           'minimum_spanning_tree']

from _components import connected_components
from _laplacian import laplacian
from _shortest_path import shortest_path, floyd_warshall, dijkstra
from _traversal import breadth_first_order, depth_first_order, \
    breadth_first_tree, depth_first_tree
from _min_spanning_tree import minimum_spanning_tree
from tools import construct_dist_matrix, reconstruct_path


def cs_graph_components(*args, **kwargs):
    """Deprecated function.  Use csgraph.connected_components"""
    warnings.warn('scipy.sparse.cs_graph_components has been deprecated. '
                  'Use scipy.sparse.csgraph.connected_components instead')
    return connected_components(*args, **kwargs)
