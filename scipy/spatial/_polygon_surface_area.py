"""

Polygon Surface Area Code

.. versionadded:: 1.0.0

"""

from __future__ import division, print_function, absolute_import
import numpy as np
from scipy.spatial.distance import pdist, cdist
from . import _surface_area
from six.moves import xrange
import math

#
# Copyright (C)  James Nichols and Tyler Reddy
#
# Distributed under the same BSD license as Scipy.

def poly_area(vertices,
              radius=None,
              threshold=1e-21,
              discretizations=500,
              n_rot=50,
              n_polygons=1):

    if n_polygons == 1:
        area = _surface_area.poly_area(vertices=vertices,
                                       radius=radius,
                                       threshold=threshold,
                                       discretizations=discretizations,
                                       n_rot=n_rot)
    else:
        area = np.empty(n_polygons)
        for polygon in xrange(n_polygons):
            area[polygon] = _surface_area.poly_area(vertices=vertices[polygon],
                                                    radius=radius,
                                                    threshold=threshold,
                                                    discretizations=discretizations,
                                                    n_rot=n_rot)

    return area
