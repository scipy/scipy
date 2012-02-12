"""Sparse Graph Routines"""
import numpy
from info import __doc__

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

cs_graph_components = numpy.deprecate(connected_components,
                                      'cs_graph_components',
                                      'csgraph.connected_components')
