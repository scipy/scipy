"""
Vector Quantization / Kmeans
============================
Clustering algorithms are useful in information theory, target detection,
communications, compression, and other areas.  The `vq` module only
supports vector quantization and the k-means algorithms.  Development of
self-organizing maps (SOM) and other approaches is underway.

Hierarchical Clustering
=======================
The `hierarchy` module provides functions for hierarchical and
agglomerative clustering.  Its features include generating hierarchical
clusters from distance matrices, computing distance matrices from
observation vectors, calculating statistics on clusters, cutting linkages
to generate flat clusters, and visualizing clusters with dendrograms.

"""
#
# spatial - Distances
#

from info import __doc__

__all__ = ['vq', 'hierarchy']

import vq, hierarchy

from numpy.testing import Tester
test = Tester().test
