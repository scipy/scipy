import numpy as np

from .common import Benchmark, LimitedParamBenchmark, safe_import

with safe_import():
    from scipy.spatial import cKDTree, KDTree
with safe_import():
    from scipy.spatial import distance
with safe_import():
    from scipy.spatial import ConvexHull, Voronoi
with safe_import():
    from scipy.spatial import SphericalVoronoi
with safe_import():
    from scipy.spatial import geometric_slerp
with safe_import():
    from scipy.spatial.transform import Rotation


class Build(Benchmark):
    params = [
        [(3,10000,1000), (8,10000,1000), (16,10000,1000)],
        ['KDTree', 'cKDTree'],
    ]
    param_names = ['(m, n, r)', 'class']

    def setup(self, mnr, cls_name):
        self.cls = KDTree if cls_name == 'KDTree' else cKDTree
        m, n, r = mnr

        np.random.seed(1234)
        self.data = np.concatenate((np.random.randn(n//2,m),
                                    np.random.randn(n-n//2,m)+np.ones(m)))

        self.queries = np.concatenate((np.random.randn(r//2,m),
                                       np.random.randn(r-r//2,m)+np.ones(m)))

    def time_build(self, mnr, cls_name):
        """
        Constructing kd-tree
        =======================
        dim | # points |  time
        """
        m, n, r = mnr
        if cls_name == 'cKDTree_flat':
            self.T = self.cls(self.data, leafsize=n)
        else:
            self.cls(self.data)


class PresortedDataSetup(Benchmark):
    params = [
        [(3, 10 ** 4, 1000), (8, 10 ** 4, 1000), (16, 10 ** 4, 1000)],
        [True, False],
        ['random', 'sorted'],
        [0.5]
    ]
    param_names = ['(m, n, r)', 'balanced', 'order', 'radius']

    def setup(self, mnr, balanced, order, radius):
        m, n, r = mnr

        np.random.seed(1234)
        self.data = {
            'random': np.random.uniform(size=(n, m)),
            'sorted': np.repeat(np.arange(n, 0, -1)[:, np.newaxis],
                                m,
                                axis=1) / n
        }

        self.queries = np.random.uniform(size=(r, m))
        self.T = cKDTree(self.data.get(order), balanced_tree=balanced)


class BuildUnbalanced(PresortedDataSetup):
    params = PresortedDataSetup.params[:-1]
    param_names = PresortedDataSetup.param_names[:-1]

    def setup(self, *args):
        super().setup(*args, None)

    def time_build(self, mnr, balanced, order):
        cKDTree(self.data.get(order), balanced_tree=balanced)


class QueryUnbalanced(PresortedDataSetup):
    params = PresortedDataSetup.params[:-1]
    param_names = PresortedDataSetup.param_names[:-1]

    def setup(self, *args):
        super().setup(*args, None)

    def time_query(self, mnr, balanced, order):
        self.T.query(self.queries)


class RadiusUnbalanced(PresortedDataSetup):
    params = PresortedDataSetup.params[:]
    params[0] = [(3, 1000, 30), (8, 1000, 30), (16, 1000, 30)]

    def time_query_pairs(self, mnr, balanced, order, radius):
        self.T.query_pairs(radius)

    def time_query_ball_point(self, mnr, balanced, order, radius):
        self.T.query_ball_point(self.queries, radius)


LEAF_SIZES = [8, 128]
BOX_SIZES = [None, 0.0, 1.0]

class Query(LimitedParamBenchmark):
    params = [
        [(3,10000,1000), (8,10000,1000), (16,10000,1000)],
        [1, 2, np.inf],
        BOX_SIZES, LEAF_SIZES,
    ]
    param_names = ['(m, n, r)', 'p', 'boxsize', 'leafsize']
    num_param_combinations = 21

    @staticmethod
    def do_setup(self, mnr, p, boxsize, leafsize):
        m, n, r = mnr

        np.random.seed(1234)

        self.data = np.random.uniform(size=(n, m))
        self.queries = np.random.uniform(size=(r, m))

        self.T = cKDTree(self.data, leafsize=leafsize, boxsize=boxsize)

    def setup(self, mnr, p, boxsize, leafsize):
        LimitedParamBenchmark.setup(self, mnr, p, boxsize, leafsize)
        Query.do_setup(self, mnr, p, boxsize, leafsize)

    def time_query(self, mnr, p, boxsize, leafsize):
        """
        Querying kd-tree
        dim | # points | # queries |  KDTree  | cKDTree | flat cKDTree
        """
        self.T.query(self.queries, p=p)

    # Retain old benchmark results (remove this if changing the benchmark)
    time_query.version = "327bc0627d5387347e9cdcf4c52a550c813bb80a859eeb0f3e5bfe6650a8a1db"


class Radius(LimitedParamBenchmark):
    params = [
        [(3,10000,1000)],
        [1, 2, np.inf],
        [0.2, 0.5],
        BOX_SIZES, LEAF_SIZES,
    ]
    param_names = ['(m, n, r)', 'p', 'probe radius', 'boxsize', 'leafsize']
    num_param_combinations = 7

    def __init__(self):
        self.time_query_pairs.__func__.params = list(self.params)
        self.time_query_pairs.__func__.params[0] = [(3,1000,30),
                                                    (8,1000,30),
                                                    (16,1000,30)]
        self.time_query_ball_point.__func__.setup = self.setup_query_ball_point
        self.time_query_ball_point_nosort.__func__.setup = self.setup_query_ball_point
        self.time_query_pairs.__func__.setup = self.setup_query_pairs

    def setup(self, *args):
        pass

    def setup_query_ball_point(self, mnr, p, probe_radius, boxsize, leafsize):
        LimitedParamBenchmark.setup(self, mnr, p, probe_radius, boxsize, leafsize,
                                    param_seed=3)
        Query.do_setup(self, mnr, p, boxsize, leafsize)

    def setup_query_pairs(self, mnr, p, probe_radius, boxsize, leafsize):
        # query_pairs is fast enough so we can run all parameter combinations
        Query.do_setup(self, mnr, p, boxsize, leafsize)

    def time_query_ball_point(self, mnr, p, probe_radius, boxsize, leafsize):
        self.T.query_ball_point(self.queries, probe_radius, p=p)

    def time_query_ball_point_nosort(self, mnr, p, probe_radius, boxsize, leafsize):
        self.T.query_ball_point(self.queries, probe_radius, p=p,
                                return_sorted=False)

    def time_query_pairs(self, mnr, p, probe_radius, boxsize, leafsize):
        self.T.query_pairs(probe_radius, p=p)

    # Retain old benchmark results (remove this if changing the benchmark)
    time_query_ball_point.version = "e0c2074b35db7e5fca01a43b0fba8ab33a15ed73d8573871ea6feb57b3df4168"
    time_query_pairs.version = "cf669f7d619e81e4a09b28bb3fceaefbdd316d30faf01524ab33d41661a53f56"


class Neighbors(LimitedParamBenchmark):
    params = [
        [(3,1000,1000),
         (8,1000,1000),
         (16,1000,1000)],
        [1, 2, np.inf],
        [0.2, 0.5],
        BOX_SIZES, LEAF_SIZES,
        ['cKDTree', 'cKDTree_weighted'],
    ]
    param_names = ['(m, n1, n2)', 'p', 'probe radius', 'boxsize', 'leafsize', 'cls']
    num_param_combinations = 17

    def setup(self, mn1n2, p, probe_radius, boxsize, leafsize, cls):
        LimitedParamBenchmark.setup(self, mn1n2, p, probe_radius, boxsize, leafsize, cls)

        m, n1, n2 = mn1n2

        self.data1 = np.random.uniform(size=(n1, m))
        self.data2 = np.random.uniform(size=(n2, m))

        self.w1 = np.ones(n1)
        self.w2 = np.ones(n2)

        self.T1 = cKDTree(self.data1, boxsize=boxsize, leafsize=leafsize)
        self.T2 = cKDTree(self.data2, boxsize=boxsize, leafsize=leafsize)

    def time_sparse_distance_matrix(self, mn1n2, p, probe_radius, boxsize, leafsize, cls):
        self.T1.sparse_distance_matrix(self.T2, probe_radius, p=p)

    def time_count_neighbors(self, mn1n2, p, probe_radius, boxsize, leafsize, cls):
        """
        Count neighbors kd-tree
        dim | # points T1 | # points T2 | p | probe radius |  BoxSize | LeafSize | cls
        """

        if cls != 'cKDTree_weighted':
            self.T1.count_neighbors(self.T2, probe_radius, p=p)
        else:
            self.T1.count_neighbors(self.T2, probe_radius, weights=(self.w1, self.w2), p=p)

    # Retain old benchmark results (remove this if changing the benchmark)
    time_sparse_distance_matrix.version = "9aa921dce6da78394ab29d949be27953484613dcf9c9632c01ae3973d4b29596"
    time_count_neighbors.version = "830287f1cf51fa6ba21854a60b03b2a6c70b2f2485c3cdcfb19a360e0a7e2ca2"


class CNeighbors(Benchmark):
    params = [
        [
          (2,1000,1000),
          (8,1000,1000),
          (16,1000,1000)
        ],
        [2, 10, 100, 400, 1000],
    ]
    param_names = ['(m, n1, n2)', 'Nr']

    def setup(self, mn1n2, Nr):
        m, n1, n2 = mn1n2

        data1 = np.random.uniform(size=(n1, m))
        data2 = np.random.uniform(size=(n2, m))
        self.w1 = np.ones(len(data1))
        self.w2 = np.ones(len(data2))

        self.T1d = cKDTree(data1, leafsize=1)
        self.T2d = cKDTree(data2, leafsize=1)
        self.T1s = cKDTree(data1, leafsize=8)
        self.T2s = cKDTree(data2, leafsize=8)
        self.r = np.linspace(0, 0.5, Nr)

    def time_count_neighbors_deep(self, mn1n2, Nr):
        """
        Count neighbors for a very deep kd-tree
        dim | # points T1 | # points T2 | Nr
        """
        self.T1d.count_neighbors(self.T2d, self.r)

    def time_count_neighbors_shallow(self, mn1n2, Nr):
        """
        Count neighbors for a shallow kd-tree
        dim | # points T1 | # points T2 | Nr
        """
        self.T1s.count_neighbors(self.T2s, self.r)

def generate_spherical_points(num_points):
    # generate uniform points on sphere
    # see: https://stackoverflow.com/a/23785326
    np.random.seed(123)
    points = np.random.normal(size=(num_points, 3))
    points /= np.linalg.norm(points, axis=1)[:, np.newaxis]
    return points

class SphericalVor(Benchmark):
    params = [10, 100, 1000, 5000, 10000]
    param_names = ['num_points']

    def setup(self, num_points):
        self.points = generate_spherical_points(num_points)

    def time_spherical_voronoi_calculation(self, num_points):
        """Perform spherical Voronoi calculation, but not the sorting of
        vertices in the Voronoi polygons.
        """
        SphericalVoronoi(self.points, radius=1, center=np.zeros(3))

class SphericalVorSort(Benchmark):
    params = [10, 100, 1000, 5000, 10000]
    param_names = ['num_points']

    def setup(self, num_points):
        self.points = generate_spherical_points(num_points)
        self.sv = SphericalVoronoi(self.points, radius=1,
                                   center=np.zeros(3))

    def time_spherical_polygon_vertex_sorting(self, num_points):
        """Time the vertex sorting operation in the Spherical Voronoi
        code.
        """
        self.sv.sort_vertices_of_regions()


class SphericalVorAreas(Benchmark):
    params = [10, 100, 1000, 5000, 10000]
    param_names = ['num_points']

    def setup(self, num_points):
        self.points = generate_spherical_points(num_points)
        self.sv = SphericalVoronoi(self.points, radius=1,
                                   center=np.zeros(3))

    def time_spherical_polygon_area_calculation(self, num_points):
        """Time the area calculation in the Spherical Voronoi code."""
        self.sv.calculate_areas()


class Xdist(Benchmark):
    params = ([10, 100, 1000], ['euclidean', 'minkowski', 'cityblock',
    'seuclidean', 'sqeuclidean', 'cosine', 'correlation', 'hamming', 'jaccard',
    'jensenshannon', 'chebyshev', 'canberra', 'braycurtis', 'mahalanobis',
    'yule', 'dice', 'kulsinski', 'rogerstanimoto', 'russellrao',
    'sokalmichener', 'sokalsneath', 'wminkowski', 'minkowski-P3'])
    param_names = ['num_points', 'metric']

    def setup(self, num_points, metric):
        np.random.seed(123)
        self.points = np.random.random_sample((num_points, 3))
        self.metric = metric
        if metric == 'minkowski-P3':
            # p=2 is just the euclidean metric, try another p value as well
            self.kwargs = {'p': 3.0}
            self.metric = 'minkowski'
        elif metric == 'wminkowski':
            # use an equal weight vector since weights are required
            self.kwargs = {'w': np.ones(3)}
        else:
            self.kwargs = {}

    def time_cdist(self, num_points, metric):
        """Time scipy.spatial.distance.cdist over a range of input data
        sizes and metrics.
        """
        distance.cdist(self.points, self.points, self.metric, **self.kwargs)

    def time_pdist(self, num_points, metric):
        """Time scipy.spatial.distance.pdist over a range of input data
        sizes and metrics.
        """
        distance.pdist(self.points, self.metric, **self.kwargs)


class XdistWeighted(Benchmark):
    params = (
        [10, 20, 100],
        ['euclidean', 'minkowski', 'cityblock', 'sqeuclidean', 'cosine',
         'correlation', 'hamming', 'jaccard', 'chebyshev', 'canberra',
         'braycurtis', 'yule', 'dice', 'kulsinski', 'rogerstanimoto',
         'russellrao', 'sokalmichener', 'sokalsneath', 'minkowski-P3'])
    param_names = ['num_points', 'metric']

    def setup(self, num_points, metric):
        np.random.seed(123)
        self.points = np.random.random_sample((num_points, 3))
        self.metric = metric
        if metric == 'minkowski-P3':
            # p=2 is just the euclidean metric, try another p value as well
            self.kwargs = {'p': 3.0}
            self.metric = 'minkowski'
        else:
            self.kwargs = {}
        self.weights = np.ones(3)

    def time_cdist(self, num_points, metric):
        """Time scipy.spatial.distance.cdist for weighted distance metrics."""
        distance.cdist(self.points, self.points, self.metric, w=self.weights,
                       **self.kwargs)

    def time_pdist(self, num_points, metric):
        """Time scipy.spatial.distance.pdist for weighted distance metrics."""
        distance.pdist(self.points, self.metric, w=self.weights, **self.kwargs)


class ConvexHullBench(Benchmark):
    params = ([10, 100, 1000, 5000], [True, False])
    param_names = ['num_points', 'incremental']

    def setup(self, num_points, incremental):
        np.random.seed(123)
        self.points = np.random.random_sample((num_points, 3))

    def time_convex_hull(self, num_points, incremental):
        """Time scipy.spatial.ConvexHull over a range of input data sizes
        and settings.
        """
        ConvexHull(self.points, incremental)


class VoronoiBench(Benchmark):
    params = ([10, 100, 1000, 5000, 10000], [False, True])
    param_names = ['num_points', 'furthest_site']

    def setup(self, num_points, furthest_site):
        np.random.seed(123)
        self.points = np.random.random_sample((num_points, 3))

    def time_voronoi_calculation(self, num_points, furthest_site):
        """Time conventional Voronoi diagram calculation."""
        Voronoi(self.points, furthest_site=furthest_site)

class Hausdorff(Benchmark):
    params = [10, 100, 1000]
    param_names = ['num_points']

    def setup(self, num_points):
        np.random.seed(123)
        self.points1 = np.random.random_sample((num_points, 3))
        np.random.seed(71890)
        self.points2 = np.random.random_sample((num_points, 3))

    def time_directed_hausdorff(self, num_points):
        # time directed_hausdorff code in 3 D
        distance.directed_hausdorff(self.points1, self.points2)

class GeometricSlerpBench(Benchmark):
    params = [10, 1000, 10000]
    param_names = ['num_points']

    def setup(self, num_points):
        points = generate_spherical_points(50)
        # any two points from the random spherical points
        # will suffice for the interpolation bounds:
        self.start = points[0]
        self.end = points[-1]
        self.t = np.linspace(0, 1, num_points)

    def time_geometric_slerp_3d(self, num_points):
        # time geometric_slerp() for 3D interpolation
        geometric_slerp(start=self.start,
                        end=self.end,
                        t=self.t)

class RotationBench(Benchmark):
    params = [1, 10, 1000, 10000]
    param_names = ['num_rotations']

    def setup(self, num_rotations):
        np.random.seed(1234)
        self.rotations = Rotation.random(num_rotations)

    def time_matrix_conversion(self, num_rotations):
        '''Time converting rotation from and to matrices'''
        Rotation.from_matrix(self.rotations.as_matrix())

    def time_euler_conversion(self, num_rotations):
        '''Time converting rotation from and to euler angles'''
        Rotation.from_euler("XYZ", self.rotations.as_euler("XYZ"))

    def time_rotvec_conversion(self, num_rotations):
        '''Time converting rotation from and to rotation vectors'''
        Rotation.from_rotvec(self.rotations.as_rotvec())

    def time_mrp_conversion(self, num_rotations):
        '''Time converting rotation from and to Modified Rodrigues Parameters'''
        Rotation.from_mrp(self.rotations.as_mrp())

    def time_mul_inv(self, num_rotations):
        '''Time multiplication and inverse of rotations'''
        self.rotations * self.rotations.inv()
