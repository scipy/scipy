"""

Polygon Surface Area Code

.. versionadded:: 1.0.0

"""

import numpy as np

#
# Copyright (C)  James Nichols and Tyler Reddy
# Contributions from Tyler Reddy are owned by LANL
# and respect the scipy distribution license
#
# Distributed under the same BSD license as Scipy.

def _wrap_vertices(vertices, num_vertices):
    # wrap vertices so that last point is also first point
    new_vertices = np.zeros((num_vertices + 1, 3))
    new_vertices[:-1,...] = vertices
    new_vertices[-1,...] = vertices[0,...]
    return new_vertices

def poly_area(vertices, radius=None):
    # calculate the surface area of a planar or spherical polygon
    # crude pure Python implementation for handling a single
    # polygon at a time
    # based on JPL Publication 07-3 by Chamberlain and Duquette (2007)
    # for planar polygons we currently still require x,y,z coords
    # can just set i.e., z = 0 for all vertices
    # TODO: check for duplicate vertices in input (and corresponding
    # unit test?)
    num_vertices = vertices.shape[0]
    new_vertices = _wrap_vertices(vertices, num_vertices)
    area_sum = 0

    if radius is not None: # spherical polygons
        if radius <= 0.0:
            raise ValueError('radius must be > 0.0')

        lambda_vals = np.arctan2(new_vertices[...,1], new_vertices[...,0]) # longitudes
        phi_vals = np.arcsin(new_vertices[...,2] / radius) # latitudes

        for i in xrange(0, num_vertices):
            delta_lambda = (lambda_vals[i + 1] -
                            lambda_vals[i - 1])
            area_sum += delta_lambda * np.sin(phi_vals[i])
        
        # the paper divides by 2 here, but my testing
        # suggests we should not do that!
        area = (radius ** 2) * area_sum
    else: # planar polygon
        for i in xrange(0, num_vertices):
            delta_x = (new_vertices[i + 1][0] -
                       new_vertices[i - 1][0])
            area_sum += delta_x * new_vertices[i][1] 
        area = -0.5 * area_sum
    
    return area






