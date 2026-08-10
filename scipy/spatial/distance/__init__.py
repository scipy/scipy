"""
Distance computations (:mod:`scipy.spatial.distance`)
=====================================================

.. sectionauthor:: Damian Eads

Function reference
------------------

Distance matrix computation from a collection of raw observation vectors
stored in a rectangular array.

.. autosummary::
   :toctree: generated/

   pdist   -- pairwise distances between observation vectors.
   cdist   -- distances between two collections of observation vectors
   squareform -- convert distance matrix to a condensed one and vice versa
   directed_hausdorff -- directed Hausdorff distance between arrays

Predicates for checking the validity of distance matrices, both
condensed and redundant. Also contained in this module are functions
for computing the number of observations in a distance matrix.

.. autosummary::
   :toctree: generated/

   is_valid_dm -- checks for a valid distance matrix
   is_valid_y  -- checks for a valid condensed distance matrix
   num_obs_dm  -- # of observations in a distance matrix
   num_obs_y   -- # of observations in a condensed distance matrix

Distance functions between two numeric vectors ``u`` and ``v``. Computing
distances over a large collection of vectors is inefficient for these
functions. Use ``pdist`` for this purpose.

.. autosummary::
   :toctree: generated/

   braycurtis       -- the Bray-Curtis distance.
   canberra         -- the Canberra distance.
   chebyshev        -- the Chebyshev distance.
   cityblock        -- the Manhattan distance.
   correlation      -- the Correlation distance.
   cosine           -- the Cosine distance.
   euclidean        -- the Euclidean distance.
   jensenshannon    -- the Jensen-Shannon distance.
   mahalanobis      -- the Mahalanobis distance.
   minkowski        -- the Minkowski distance.
   seuclidean       -- the normalized Euclidean distance.
   sqeuclidean      -- the squared Euclidean distance.

Distance functions between two boolean vectors (representing sets) ``u`` and
``v``.  As in the case of numerical vectors, ``pdist`` is more efficient for
computing the distances between all pairs.

.. autosummary::
   :toctree: generated/

   dice             -- the Dice dissimilarity.
   hamming          -- the Hamming distance.
   jaccard          -- the Jaccard distance.
   rogerstanimoto   -- the Rogers-Tanimoto dissimilarity.
   russellrao       -- the Russell-Rao dissimilarity.
   sokalsneath      -- the Sokal-Sneath dissimilarity.
   yule             -- the Yule dissimilarity.

:func:`hamming` also operates over discrete numerical vectors.
"""

from ._distance import (
    braycurtis,
    canberra,
    cdist,
    chebyshev,
    cityblock,
    correlation,
    cosine,
    dice,
    directed_hausdorff,
    euclidean,
    hamming,
    is_valid_dm,
    is_valid_y,
    jaccard,
    jensenshannon,
    mahalanobis,
    minkowski,
    num_obs_dm,
    num_obs_y,
    pdist,
    rogerstanimoto,
    russellrao,
    seuclidean,
    sokalsneath,
    sqeuclidean,
    squareform,
    yule,
)

__all__ = [
    "braycurtis",
    "canberra",
    "cdist",
    "chebyshev",
    "cityblock",
    "correlation",
    "cosine",
    "dice",
    "directed_hausdorff",
    "euclidean",
    "hamming",
    "is_valid_dm",
    "is_valid_y",
    "jaccard",
    "jensenshannon",
    "mahalanobis",
    "minkowski",
    "num_obs_dm",
    "num_obs_y",
    "pdist",
    "rogerstanimoto",
    "russellrao",
    "seuclidean",
    "sokalsneath",
    "sqeuclidean",
    "squareform",
    "yule",
]
