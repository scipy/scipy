from __future__ import division, print_function, absolute_import

import numpy as np
from scipy.cluster import vq


def cluster_data(data,cluster_cnt,iter=20,thresh=1e-5):
    """ Group data into a number of common clusters

        data -- 2D array of data points.  Each point is a row in the array.
        cluster_cnt -- The number of clusters to use
        iter -- number of iterations to use for kmeans algorithm
        thresh -- distortion threshold for kmeans algorithm

        return -- list of 2D arrays.  Each array contains the data points
                  that belong to a specific cluster.

        Uses kmeans algorithm to find the clusters.
    """
    wh_data = vq.whiten(data)
    code_book,dist = vq.kmeans(wh_data,cluster_cnt,iter,thresh)
    code_ids, distortion = vq.vq(wh_data,code_book)
    clusters = []
    for i in range(len(code_book)):
        cluster = np.compress(code_ids == i,data,0)
        clusters.append(cluster)
    return clusters

if __name__ == "__main__":

    data = np.array(((400, 79, 5.4),
                     (180, 76, 4.5),
                     (28, 25, 30.),
                     (270, 81, 5.0),
                     (185, 78, 4.6)))

    clusters = cluster_data(data,2)
    for i in range(len(clusters)):
        print('cluster %d:' % i)
        print(clusters[i])
