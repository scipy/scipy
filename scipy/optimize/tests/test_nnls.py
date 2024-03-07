import numpy as np
from numpy.testing import assert_allclose
from pytest import raises as assert_raises
from scipy.optimize import nnls


class TestNNLS:
    def setup_method(self):
        self.rng = np.random.default_rng(1685225766635251)

    def test_nnls(self):
        a = np.arange(25.0).reshape(-1, 5)
        x = np.arange(5.0)
        y = a @ x
        x, res = nnls(a, y)
        assert res < 1e-7
        assert np.linalg.norm((a @ x) - y) < 1e-7

    def test_nnls_tall(self):
        a = self.rng.uniform(low=-10, high=10, size=[50, 10])
        x = np.abs(self.rng.uniform(low=-2, high=2, size=[10]))
        x[::2] = 0
        b = a @ x
        xact, rnorm = nnls(a, b, atol=500*np.linalg.norm(a, 1)*np.spacing(1.))
        assert_allclose(xact, x, rtol=0., atol=1e-10)
        assert rnorm < 1e-12

    def test_nnls_wide(self):
        # If too wide then problem becomes too ill-conditioned ans starts
        # emitting warnings, hence small m, n difference.
        a = self.rng.uniform(low=-10, high=10, size=[100, 120])
        x = np.abs(self.rng.uniform(low=-2, high=2, size=[120]))
        x[::2] = 0
        b = a @ x
        xact, rnorm = nnls(a, b, atol=500*np.linalg.norm(a, 1)*np.spacing(1.))
        assert_allclose(xact, x, rtol=0., atol=1e-10)
        assert rnorm < 1e-12

    def test_maxiter(self):
        # test that maxiter argument does stop iterations
        a = self.rng.uniform(size=(5, 10))
        b = self.rng.uniform(size=5)
        with assert_raises(RuntimeError):
            nnls(a, b, maxiter=1)


    def test_nnls_max_iterations(self):
        # This is from a real test case that failed due to a problem in the
        # scipy._nnls implementation
        n = np.array(
            [3, 2, 0, 1, 1, 1, 3, 8, 14, 16, 29, 23, 41, 47, 53, 57, 67, 76,
             103, 89, 97, 94, 85, 95, 78, 78, 78, 77, 73, 50, 50, 56, 68, 98,
             95, 112, 134, 145, 158, 172, 213, 234, 222, 215, 216, 216, 206,
             183, 135, 156, 110, 92, 63, 60, 52, 29, 20, 16, 12, 5, 5, 5, 1, 2,
             3, 0, 2])
        k = np.array(
            [0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,
             0., 0., 0., 0.7205812007860187, 0., 1.4411624015720375,
             0.7205812007860187, 2.882324803144075, 5.76464960628815,
             5.76464960628815, 12.249880413362318, 15.132205216506394,
             20.176273622008523, 27.382085629868712, 48.27894045266326,
             47.558359251877235, 68.45521407467177, 97.99904330689854,
             108.0871801179028, 135.46926574777152, 140.51333415327366,
             184.4687874012208, 171.49832578707245, 205.36564222401535,
             244.27702706646033, 214.01261663344755, 228.42424064916793,
             232.02714665309804, 205.36564222401535, 172.9394881886445,
             191.67459940908097, 162.1307701768542, 153.48379576742198,
             110.96950492104689, 103.04311171240067, 86.46974409432225,
             60.528820866025576, 43.234872047161126, 23.779179625938617,
             24.499760826724636, 17.29394881886445, 11.5292992125763,
             5.76464960628815, 5.044068405502131, 3.6029060039300935, 0.,
             2.882324803144075, 0., 0., 0.])
        d = np.array(
            [0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,
             0., 0., 0., 0.003889242101538, 0., 0.007606268390096, 0.,
             0.025457371599973, 0.036952882091577, 0., 0.08518359183449,
             0.048201126400243, 0.196234990022205, 0.144116240157247,
             0.171145134062442, 0., 0., 0.269555036538714, 0., 0., 0.,
             0.010893241091872, 0., 0., 0., 0., 0., 0., 0., 0.,
             0.048167058272886, 0.011238724891049, 0., 0., 0.055162603456078,
             0., 0., 0., 0., 0.027753339088588, 0., 0., 0., 0., 0., 0., 0., 0.,
             0., 0.])
        # The following code sets up a system of equations such that
        # $k_i-p_i*n_i$ is minimized for $p_i$ with weights $n_i$ and
        # monotonicity constraints on $p_i$. This translates to a system of
        # equations of the form $k_i - (d_1 + ... + d_i) * n_i$ and
        # non-negativity constraints on the $d_i$. If $n_i$ is zero the
        # system is modified such that $d_i - d_{i+1}$ is then minimized.
        N = len(n)
        A = np.diag(n) @ np.tril(np.ones((N, N)))
        w = n ** 0.5

        nz = (n == 0).nonzero()[0]
        A[nz, nz] = 1
        A[nz, np.minimum(nz + 1, N - 1)] = -1
        w[nz] = 1
        k[nz] = 0
        W = np.diag(w)

        # Small perturbations can already make the infinite loop go away (just
        # uncomment the next line)
        # k = k + 1e-10 * np.random.normal(size=N)
        dact, _ = nnls(W @ A, W @ k)
        assert_allclose(dact, d, rtol=0., atol=1e-10)

    def test_nnls_max_iterations2(self):
        # This is from a real test case that failed due to a problem in the
        # scipy._nnls implementation

        n = np.array(
            [1, 0, 1, 2, 2, 2, 3, 3, 5, 4, 14, 14, 19, 26, 36, 42, 36, 64, 64,
             64, 81, 85, 85, 95, 95, 95, 75, 76, 69, 81, 62, 59, 68, 64, 71, 67,
             74, 78, 118, 135, 153, 159, 210, 195, 218, 243, 236, 215, 196, 175,
             185, 149, 144, 103, 104, 75, 56, 40, 32, 26, 17, 9, 12, 8, 2, 1, 1,
             1])
        k = np.array(
            [0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,
             0., 0., 0., 0., 0., 0.7064355064917867, 0., 0., 2.11930651947536,
             0.7064355064917867, 0., 3.5321775324589333, 7.064355064917867,
             11.302968103868587, 16.95445215580288, 20.486629688261814,
             20.486629688261814, 37.44108184406469, 55.808405012851146,
             78.41434122058831, 103.13958394780086, 105.965325973768,
             125.74552015553803, 149.057891869767, 176.60887662294667,
             197.09550631120848, 211.930651947536, 204.86629688261814,
             233.8301526487814, 221.1143135319292, 195.6826352982249,
             197.80194181770025, 191.4440222592742, 187.91184472681525,
             144.11284332432447, 131.39700420747232, 116.5618585711448,
             93.24948685691584, 89.01087381796512, 53.68909849337579,
             45.211872415474346, 31.083162285638615, 24.72524272721253,
             16.95445215580288, 9.890097090885014, 9.890097090885014,
             2.8257420259671466, 2.8257420259671466, 1.4128710129835733,
             0.7064355064917867, 1.4128710129835733])
        d = np.array(
            [0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,
             0., 0., 0., 0., 0., 0.0021916146355674473, 0., 0.,
             0.011252740799789484, 0., 0., 0.037746623295934395,
             0.03602328132946222, 0.09509167709829734, 0.10505765870204821,
             0.01391037014274718, 0.0188296228752321, 0.20723559202324254,
             0.3056220879462608, 0.13304643490426477, 0., 0., 0., 0., 0., 0.,
             0., 0., 0., 0., 0., 0.043185876949706214, 0.0037266261379722554,
             0., 0., 0., 0., 0., 0.094797899357143, 0., 0., 0., 0., 0., 0., 0.,
             0., 0.23450935613672663, 0., 0., 0.07064355064917871])
        # The following code sets up a system of equations such that
        # $k_i-p_i*n_i$ is minimized for $p_i$ with weights $n_i$ and
        # monotonicity constraints on $p_i$. This translates to a system of
        # equations of the form $k_i - (d_1 + ... + d_i) * n_i$ and
        # non-negativity constraints on the $d_i$. If $n_i$ is zero the
        # system is modified such that $d_i - d_{i+1}$ is then minimized.
        N = len(n)
        A = np.diag(n) @ np.tril(np.ones((N, N)))
        w = n ** 0.5

        nz = (n == 0).nonzero()[0]
        A[nz, nz] = 1
        A[nz, np.minimum(nz + 1, N - 1)] = -1
        w[nz] = 1
        k[nz] = 0
        W = np.diag(w)

        dact, _ = nnls(W @ A, W @ k, atol=1e-7)

        p = np.cumsum(dact)
        assert np.all(dact >= 0)
        assert np.linalg.norm(k - n * p, ord=np.inf) < 28

        assert_allclose(dact, d, rtol=0., atol=1e-10)
