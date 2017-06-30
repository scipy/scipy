"""

Polygon Surface Area Code

.. versionadded:: 1.0.0

"""

from __future__ import division, print_function, absolute_import
import numpy as np
from scipy.spatial.distance import pdist
from . import _surface_area
from six.moves import xrange

#
# Copyright (C)  James Nichols and Tyler Reddy
# Contributions from Tyler Reddy are owned by LANL
# and respect the scipy distribution license
#
# Distributed under the same BSD license as Scipy.

def _vertex_index_strider(index, num_vertices):
    # handle the wrapping / iteration over
    # polygon vertices in either CW or CCW
    # sort order
    forward_index = index + 1
    backward_index = index - 1
    if forward_index >= num_vertices:
        forward_index = 0
    return forward_index, backward_index

def poly_area(vertices, radius=None, threshold=1e-21,
              cython=None):
    # calculate the surface area of a planar or spherical polygon
    # crude pure Python implementation for handling a single
    # polygon at a time
    # based on JPL Publication 07-3 by Chamberlain and Duquette (2007)
    # for planar polygons we currently still require x,y,z coords
    # can just set i.e., z = 0 for all vertices
    min_vertex_dist = pdist(vertices).min()
    if min_vertex_dist < threshold:
        err_str = '''Duplicate vertices detected based on minimum
                     distance {min_vertex_dist} and threshold value
                     {threshold}.'''.format(min_vertex_dist=min_vertex_dist,
                                            threshold=threshold)
        raise ValueError(err_str)

    num_vertices = vertices.shape[0]
    area_sum = 0

    if radius is not None: # spherical polygons
        if radius <= threshold:
            err_str = 'radius must be > {threshold}'.format(threshold=threshold)
            raise ValueError(err_str)
        if cython is None:

            lambda_vals = np.arctan2(vertices[...,1], vertices[...,0]) # longitudes
            phi_vals = np.arcsin(vertices[...,2] / radius) # latitudes

            for i in xrange(num_vertices):
                forward_index, backward_index = _vertex_index_strider(i, num_vertices)
                delta_lambda = (lambda_vals[forward_index] -
                                lambda_vals[backward_index])
                area_sum += delta_lambda * np.sin(phi_vals[i])
            # the paper divides by 2 here, but my testing
            # suggests we should not do that!
            area = (radius ** 2) * area_sum
        else: # cython code for spherical polygon SA
            area = _surface_area.spherical_polygon_area(vertices, radius)
    else: # planar polygon
        if cython is None:
            for i in xrange(num_vertices):
                forward_index, backward_index = _vertex_index_strider(i, num_vertices)
                delta_x = (vertices[forward_index][0] -
                           vertices[backward_index][0])
                area_sum += delta_x * vertices[i][1]
            area = -0.5 * area_sum
        else: # cython code for planar polygon SA
            area = _surface_area.planar_polygon_area(vertices)
    
    return abs(area)






