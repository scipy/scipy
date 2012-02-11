"""Compressed Sparse graph algorithms"""

__docformat__ = "restructuredtext en"

__all__ = ['components',
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

from _components import components
from _laplacian import laplacian
from _shortest_path import shortest_path, floyd_warshall, dijkstra
from _traversal import breadth_first_order, depth_first_order, \
    breadth_first_tree, depth_first_tree
from _min_spanning_tree import minimum_spanning_tree
from tools import construct_dist_matrix, reconstruct_path


def cs_graph_components(*args, **kwargs):
    """Deprecated function.  See csgraph.components"""
    warnings.warn('cs_graph_components has been deprecated. '
                  'Use csgraph.components instead')
    return components(*args, **kwargs)
