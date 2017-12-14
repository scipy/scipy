scipy.cluster.vq.kmeans2
========================

Introduction
============

**What is K-means Algorithm in Clustering?**

During `Cluster Analysis <https://en.wikipedia.org/wiki/Cluster_analysis>`_ during `Data Mining <https://docs.oracle.com/cd/B28359_01/datamine.111/b28129/process.htm#CHDFGCIJ>`_, data-points/features are grouped into clusters. The contents of these clusters show similarities and it makes easier to analyse data.

The `K-Means Algorithm <https://en.wikipedia.org/wiki/K-means_clustering>`_ divided N number of data points into K clusters in which each feature belongs to the nearest mean. The results are typically of the form of  `Voronoi Cells <http://mathworld.wolfram.com/VoronoiDiagram.html>`_
.

**How is scipy.cluster.vq.kmeans2 different from scipy.cluster.vq.kmeans?**

*kmeans2* is a different implementation of k-means clustering with more methods for generating initial centroids but without using a distortion change threshold as a stopping criterion as compared to *kmeans*

**Examples:**

Let us consider a sample set of features/data-points as follows

[1.9, 2.3], [1.5, 2.5], [0.8, 0.6], [1.0, 0.8], [1.65, 1.9], [0.9, 1.1]

Plotting the features::

   import matplotlib.pyplot as plt
   plt.plot([1.9,1.65,1.5,0.8,0.9,1.0],[2.3,1.9,2.5,0.6,1.1,0.8],'ro')
   plt.show()

**Syntax for scipy.cluster.vq.kmeans2**::

   scipy.cluster.vq.kmeans2(data, k, iter=10, thresh=1e-05, minit='random', missing='warn', check_finite=True)

Parameters
==========

**data** : *ndarray*

    A ‘M’ by ‘N’ array of ‘M’ observations in ‘N’ dimensions or a length ‘M’ array of ‘M’ one-dimensional observations.

**k** : *int or ndarray*

    The number of clusters to form as well as the number of centroids to generate. If minit initialization string is ‘matrix’, or if a ndarray is given instead, it is interpreted as initial cluster to use instead.

**iter** : *int, optional*

    Number of iterations of the k-means algrithm to run. Note that this differs in meaning from the iters parameter to the kmeans function.

**thresh** : *float, optional*

    (not used yet)

**minit** : *str, optional*

    Method for initialization. Available methods are ‘random’, ‘points’, and ‘matrix’:

    ‘random’: generate k centroids from a Gaussian with mean and variance estimated from the data.

    ‘points’: choose k observations (rows) at random from data for the initial centroids.

    ‘matrix’: interpret the k parameter as a k by M (or length k array for one-dimensional data) array of initial centroids.

**missing** : *str, optional*

    Method to deal with empty clusters. Available methods are ‘warn’ and ‘raise’:

    ‘warn’: give a warning and continue.

    ‘raise’: raise an ClusterError and terminate the algorithm.

**check_finite** : *bool, optional*

    Whether to check that the input matrices contain only finite numbers. Disabling may give a performance gain, but may result in problems (crashes, non-termination) if the inputs do contain infinities or NaNs. Default: True

Results
=======

**centroid** : *ndarray*

    A ‘k’ by ‘N’ array of centroids found at the last iteration of k-means.

**label** : *ndarray*

    label[i] is the code or index of the centroid the i’th observation is closest to.

**data-type** : *data-type*

    Example: int32 i.e a 32-bit integer

**Generating the Centroids to the clusters**::

   >>>from numpy import array
      from scipy.cluster.vq import kmeans2

      features = array([
                 [1.90, 2.30],
                 [1.50, 2.50],
                 [0.80, 0.60],
                 [1.00, 0.80],
                 [1.65, 1.90],
                 [0.90, 1.10]
                 ])

      kmeans2(features,2)

   >>>(array([[ 0.9       ,  0.83333333],
        [ 1.68333333,  2.23333333]]), array([1, 1, 0, 0, 1, 0], dtype=int32))

After importing the *numpy* module we declare a *numpy* array which contains the data-points and store it in the variable *features*. We also import the *kmeans2* function from *scipy.cluster.vq*.

In the above example, we use the *kmeans2* function to act on the data-set *features* containing *N=6 data-points* and divide it into *K=2* clusters.

We get the output as *(centroids, label, data-type)*

*Centroids:*

* [0.9, 0.83333333]
* [1.68333333, 2.23333333]

*Label Array:*

* [1, 1, 0, 0, 1, 0]

*Data-Type:*

* int 32
