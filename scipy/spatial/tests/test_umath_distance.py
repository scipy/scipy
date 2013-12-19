from __future__ import division, absolute_import, print_function

import scipy.spatial.distance as dist
import scipy.spatial.umath_distance as udist
import numpy as np

from numpy.testing import assert_almost_equal, run_module_suite

numerical_distances = [(udist.braycurtis, 'braycurtis', False, None),
                       (udist.chebyshev, 'chebyshev', False, None),
                       (udist.cityblock, 'cityblock', False, None),
                       (udist.correlation, 'correlation', True, None),
                       (udist.cosine, 'cosine', True, None),
                       (udist.euclidean, 'euclidean', False, None),
                       (udist.hamming, 'hamming', False, None),
                       (udist.greatcircle, None, False, 2)]

boolean_distances = [(udist.dice, 'dice', False, None),
                     (udist.hamming, 'hamming', False, None),
                     (udist.jaccard, 'jaccard', False, None),
                     (udist.kulsinski, 'kulsinski', False, None),
                     (udist.matching, 'matching', False, None),
                     (udist.rogerstanimoto, 'rogerstanimoto', False, None),
                     (udist.russellrao, 'russellrao', False, None),
                     (udist.sokalmichener, 'matching', False, None),
                     (udist.sokalsneath, 'sokalsneath', False, None),
                     (udist.yule, 'yule', False, None)]
                     
def check_pdist_vs_scipy(dist_gufunc, sp_func, data):
    assert_almost_equal(dist_gufunc(data),
                        dist.pdist(data, metric=sp_func))

def check_cdist_vs_scipy(dist_gufunc, sp_func, data_a, data_b):
    assert_almost_equal(dist_gufunc(data_a, data_b, cdist=True),
                        dist.cdist(data_a, data_b, metric=sp_func))

def check_pdist_vs_cdist(dist_gufunc, data_a, lowmem):
    n = len(data_a)
    idx = np.triu_indices(n, 1)
    assert_almost_equal(dist_gufunc(data_a),
                        dist_gufunc(data_a, data_a, cdist=True)[idx])
    if lowmem:
        assert_almost_equal(dist_gufunc(data_a, lowmem=True),
                            dist_gufunc(data_a, data_a,
                                        cdist=True, lowmem=False)[idx])
        assert_almost_equal(dist_gufunc(data_a, lowmem=False),
                            dist_gufunc(data_a, data_a,
                                        cdist=True, lowmem=True)[idx])
        assert_almost_equal(dist_gufunc(data_a, lowmem=False),
                            dist_gufunc(data_a, data_a,
                                        cdist=True, lowmem=False)[idx])
def test():
    for gufunc, spname, lowmem, dim in numerical_distances:
        data_a = np.random.rand(100, 100 if dim is None else dim)
        data_b = np.random.rand(200, 100 if dim is None else dim)
        if spname is not None:
            yield check_pdist_vs_scipy, gufunc, spname, data_a
            yield check_cdist_vs_scipy, gufunc, spname, data_a, data_b
        yield check_pdist_vs_cdist, gufunc, data_a, lowmem
    for gufunc, spname, lowmem, dim in boolean_distances:
        data_a = np.random.rand(100, 100 if dim is None else dim) < 0.5
        data_b = np.random.rand(200, 100 if dim is None else dim) < 0.5
        if spname is not None:
            yield check_pdist_vs_scipy, gufunc, spname, data_a
            yield check_cdist_vs_scipy, gufunc, spname, data_a, data_b
        yield check_pdist_vs_cdist, gufunc, data_a, lowmem

if __name__ == "__main__":
    run_module_suite()