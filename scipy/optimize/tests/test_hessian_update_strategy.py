import re
from copy import deepcopy

import numpy as np
import pytest
from numpy.linalg import norm
from numpy.testing import (TestCase, assert_array_equal, assert_array_less,
                           assert_allclose)
from scipy.optimize import (BFGS, SR1, LBFGS, minimize)


class Rosenbrock:
    """Rosenbrock function.

    The following optimization problem:
        minimize sum(100.0*(x[1:] - x[:-1]**2.0)**2.0 + (1 - x[:-1])**2.0)
    """

    def __init__(self, n=2, random_state=0):
        rng = np.random.default_rng(random_state)
        self.x0 = rng.uniform(-1, 1, n)
        self.x_opt = np.ones(n)

    def fun(self, x):
        x = np.asarray(x)
        r = np.sum(100.0 * (x[1:] - x[:-1]**2.0)**2.0 + (1 - x[:-1])**2.0,
                   axis=0)
        return r

    def grad(self, x):
        x = np.asarray(x)
        xm = x[1:-1]
        xm_m1 = x[:-2]
        xm_p1 = x[2:]
        der = np.zeros_like(x)
        der[1:-1] = (200 * (xm - xm_m1**2) -
                     400 * (xm_p1 - xm**2) * xm - 2 * (1 - xm))
        der[0] = -400 * x[0] * (x[1] - x[0]**2) - 2 * (1 - x[0])
        der[-1] = 200 * (x[-1] - x[-2]**2)
        return der

    def hess(self, x):
        x = np.atleast_1d(x)
        H = np.diag(-400 * x[:-1], 1) - np.diag(400 * x[:-1], -1)
        diagonal = np.zeros(len(x), dtype=x.dtype)
        diagonal[0] = 1200 * x[0]**2 - 400 * x[1] + 2
        diagonal[-1] = 200
        diagonal[1:-1] = 202 + 1200 * x[1:-1]**2 - 400 * x[2:]
        H = H + np.diag(diagonal)
        return H


class TestHessianUpdateStrategy(TestCase):


    def test_hessian_initialization(self):

        ndims = 5
        symmetric_matrix = np.array([[43, 24, 33, 34, 49],
                                     [24, 36, 44, 15, 44],
                                     [33, 44, 37, 1, 30],
                                     [34, 15, 1, 5, 46],
                                     [49, 44, 30, 46, 22]])
        init_scales = (
            ('auto', np.eye(ndims)),
            (2, np.eye(ndims) * 2),
            (np.arange(1, ndims + 1) * np.eye(ndims),
             np.arange(1, ndims + 1) * np.eye(ndims)),
            (symmetric_matrix, symmetric_matrix),)
        for approx_type in ['hess', 'inv_hess']:
            for init_scale, true_matrix in init_scales:
                # large min_{denominator,curvatur} makes them skip an update,
                # so we can have our initial matrix
                quasi_newton = (BFGS(init_scale=init_scale,
                                     min_curvature=1e50,
                                     exception_strategy='skip_update'),
                                SR1(init_scale=init_scale,
                                    min_denominator=1e50))

                for qn in quasi_newton:
                    qn.initialize(ndims, approx_type)
                    B = qn.get_matrix()

                    assert_array_equal(B, np.eye(ndims))
                    # don't test the auto init scale
                    if isinstance(init_scale, str) and init_scale == 'auto':
                        continue

                    qn.update(np.ones(ndims) * 1e-5, np.arange(ndims) + 0.2)
                    B = qn.get_matrix()
                    assert_array_equal(B, true_matrix)

    # For this list of points, it is known
    # that no exception occur during the
    # Hessian update. Hence no update is
    # skipped or damped.


    def test_initialize_catch_illegal(self):
        ndims = 3
        # no complex allowed
        inits_msg_errtype = ((complex(3.14),
                              r"float\(\) argument must be a string or a "
                              r"(real )?number, not 'complex'",
                              TypeError),

                             (np.array([3.2, 2.3, 1.2]).astype(np.complex128),
                              "init_scale contains complex elements, "
                              "must be real.",
                              TypeError),

                             (np.array([[43, 24, 33],
                                        [24, 36, 44, ],
                                        [33, 44, 37, ]]).astype(np.complex128),
                              "init_scale contains complex elements, "
                              "must be real.",
                              TypeError),

                             # not square
                             (np.array([[43, 55, 66]]),
                              re.escape(
                                  "If init_scale is an array, it must have the "
                                  "dimensions of the hess/inv_hess: (3, 3)."
                                  " Got (1, 3)."),
                              ValueError),

                             # not symmetric
                             (np.array([[43, 24, 33],
                                        [24.1, 36, 44, ],
                                        [33, 44, 37, ]]),
                              re.escape("If init_scale is an array, it must be"
                                        " symmetric (passing scipy.linalg.issymmetric)"
                                        " to be an approximation of a hess/inv_hess."),
                              ValueError),
                             )
        for approx_type in ['hess', 'inv_hess']:
            for init_scale, message, errortype in inits_msg_errtype:
                # large min_{denominator,curvatur} makes it skip an update,
                # so we can retrieve our initial matrix
                quasi_newton = (BFGS(init_scale=init_scale),
                                SR1(init_scale=init_scale))

                for qn in quasi_newton:
                    qn.initialize(ndims, approx_type)
                    with pytest.raises(errortype, match=message):
                        qn.update(np.ones(ndims), np.arange(ndims))

    def test_rosenbrock_with_no_exception(self):
        # Define auxiliary problem
        prob = Rosenbrock(n=5)
        # Define iteration points
        x_list = [[0.0976270, 0.4303787, 0.2055267, 0.0897663, -0.15269040],
                  [0.1847239, 0.0505757, 0.2123832, 0.0255081, 0.00083286],
                  [0.2142498, -0.0188480, 0.0503822, 0.0347033, 0.03323606],
                  [0.2071680, -0.0185071, 0.0341337, -0.0139298, 0.02881750],
                  [0.1533055, -0.0322935, 0.0280418, -0.0083592, 0.01503699],
                  [0.1382378, -0.0276671, 0.0266161, -0.0074060, 0.02801610],
                  [0.1651957, -0.0049124, 0.0269665, -0.0040025, 0.02138184],
                  [0.2354930, 0.0443711, 0.0173959, 0.0041872, 0.00794563],
                  [0.4168118, 0.1433867, 0.0111714, 0.0126265, -0.00658537],
                  [0.4681972, 0.2153273, 0.0225249, 0.0152704, -0.00463809],
                  [0.6023068, 0.3346815, 0.0731108, 0.0186618, -0.00371541],
                  [0.6415743, 0.3985468, 0.1324422, 0.0214160, -0.00062401],
                  [0.7503690, 0.5447616, 0.2804541, 0.0539851, 0.00242230],
                  [0.7452626, 0.5644594, 0.3324679, 0.0865153, 0.00454960],
                  [0.8059782, 0.6586838, 0.4229577, 0.1452990, 0.00976702],
                  [0.8549542, 0.7226562, 0.4991309, 0.2420093, 0.02772661],
                  [0.8571332, 0.7285741, 0.5279076, 0.2824549, 0.06030276],
                  [0.8835633, 0.7727077, 0.5957984, 0.3411303, 0.09652185],
                  [0.9071558, 0.8299587, 0.6771400, 0.4402896, 0.17469338],
                  [0.9190793, 0.8486480, 0.7163332, 0.5083780, 0.26107691],
                  [0.9371223, 0.8762177, 0.7653702, 0.5773109, 0.32181041],
                  [0.9554613, 0.9119893, 0.8282687, 0.6776178, 0.43162744],
                  [0.9545744, 0.9099264, 0.8270244, 0.6822220, 0.45237623],
                  [0.9688112, 0.9351710, 0.8730961, 0.7546601, 0.56622448],
                  [0.9743227, 0.9491953, 0.9005150, 0.8086497, 0.64505437],
                  [0.9807345, 0.9638853, 0.9283012, 0.8631675, 0.73812581],
                  [0.9886746, 0.9777760, 0.9558950, 0.9123417, 0.82726553],
                  [0.9899096, 0.9803828, 0.9615592, 0.9255600, 0.85822149],
                  [0.9969510, 0.9935441, 0.9864657, 0.9726775, 0.94358663],
                  [0.9979533, 0.9960274, 0.9921724, 0.9837415, 0.96626288],
                  [0.9995981, 0.9989171, 0.9974178, 0.9949954, 0.99023356],
                  [1.0002640, 1.0005088, 1.0010594, 1.0021161, 1.00386912],
                  [0.9998903, 0.9998459, 0.9997795, 0.9995484, 0.99916305],
                  [1.0000008, 0.9999905, 0.9999481, 0.9998903, 0.99978047],
                  [1.0000004, 0.9999983, 1.0000001, 1.0000031, 1.00000297],
                  [0.9999995, 1.0000003, 1.0000005, 1.0000001, 1.00000032],
                  [0.9999999, 0.9999997, 0.9999994, 0.9999989, 0.99999786],
                  [0.9999999, 0.9999999, 0.9999999, 0.9999999, 0.99999991]]
        # Get iteration points
        grad_list = [prob.grad(x) for x in x_list]
        delta_x = [np.array(x_list[i+1])-np.array(x_list[i])
                   for i in range(len(x_list)-1)]
        delta_grad = [grad_list[i+1]-grad_list[i]
                      for i in range(len(grad_list)-1)]
        # Check curvature condition
        for s, y in zip(delta_x, delta_grad):
            if np.dot(s, y) <= 0:
                raise ArithmeticError()
        # Define QuasiNewton update
        for quasi_newton in (BFGS(init_scale=1, min_curvature=1e-4),
                             SR1(init_scale=1), LBFGS(memory=40, init_scale=1)):
            hess = deepcopy(quasi_newton)
            inv_hess = deepcopy(quasi_newton)
            hess.initialize(len(x_list[0]), 'hess')
            inv_hess.initialize(len(x_list[0]), 'inv_hess')
            # Compare the hessian and its inverse
            for s, y in zip(delta_x, delta_grad):
                hess.update(s, y)
                inv_hess.update(s, y)
                B = hess.get_matrix()
                H = inv_hess.get_matrix()
                assert_allclose(np.linalg.inv(B), H, rtol=1e-9, atol=1e-10)
            B_true = prob.hess(x_list[len(delta_x)])
            assert_array_less(norm(B - B_true)/norm(B_true), 0.1)

    def test_SR1_skip_update(self):
        # Define auxiliary problem
        prob = Rosenbrock(n=5)
        # Define iteration points
        x_list = [[0.0976270, 0.4303787, 0.2055267, 0.0897663, -0.15269040],
                  [0.1847239, 0.0505757, 0.2123832, 0.0255081, 0.00083286],
                  [0.2142498, -0.0188480, 0.0503822, 0.0347033, 0.03323606],
                  [0.2071680, -0.0185071, 0.0341337, -0.0139298, 0.02881750],
                  [0.1533055, -0.0322935, 0.0280418, -0.0083592, 0.01503699],
                  [0.1382378, -0.0276671, 0.0266161, -0.0074060, 0.02801610],
                  [0.1651957, -0.0049124, 0.0269665, -0.0040025, 0.02138184],
                  [0.2354930, 0.0443711, 0.0173959, 0.0041872, 0.00794563],
                  [0.4168118, 0.1433867, 0.0111714, 0.0126265, -0.00658537],
                  [0.4681972, 0.2153273, 0.0225249, 0.0152704, -0.00463809],
                  [0.6023068, 0.3346815, 0.0731108, 0.0186618, -0.00371541],
                  [0.6415743, 0.3985468, 0.1324422, 0.0214160, -0.00062401],
                  [0.7503690, 0.5447616, 0.2804541, 0.0539851, 0.00242230],
                  [0.7452626, 0.5644594, 0.3324679, 0.0865153, 0.00454960],
                  [0.8059782, 0.6586838, 0.4229577, 0.1452990, 0.00976702],
                  [0.8549542, 0.7226562, 0.4991309, 0.2420093, 0.02772661],
                  [0.8571332, 0.7285741, 0.5279076, 0.2824549, 0.06030276],
                  [0.8835633, 0.7727077, 0.5957984, 0.3411303, 0.09652185],
                  [0.9071558, 0.8299587, 0.6771400, 0.4402896, 0.17469338]]
        # Get iteration points
        grad_list = [prob.grad(x) for x in x_list]
        delta_x = [np.array(x_list[i+1])-np.array(x_list[i])
                   for i in range(len(x_list)-1)]
        delta_grad = [grad_list[i+1]-grad_list[i]
                      for i in range(len(grad_list)-1)]
        hess = SR1(init_scale=1, min_denominator=1e-2)
        hess.initialize(len(x_list[0]), 'hess')
        # Compare the Hessian and its inverse
        for i in range(len(delta_x)-1):
            s = delta_x[i]
            y = delta_grad[i]
            hess.update(s, y)
        # Test skip update
        B = np.copy(hess.get_matrix())
        s = delta_x[17]
        y = delta_grad[17]
        hess.update(s, y)
        B_updated = np.copy(hess.get_matrix())
        assert_array_equal(B, B_updated)

    def test_BFGS_skip_update(self):
        # Define auxiliary problem
        prob = Rosenbrock(n=5)
        # Define iteration points
        x_list = [[0.0976270, 0.4303787, 0.2055267, 0.0897663, -0.15269040],
                  [0.1847239, 0.0505757, 0.2123832, 0.0255081, 0.00083286],
                  [0.2142498, -0.0188480, 0.0503822, 0.0347033, 0.03323606],
                  [0.2071680, -0.0185071, 0.0341337, -0.0139298, 0.02881750],
                  [0.1533055, -0.0322935, 0.0280418, -0.0083592, 0.01503699],
                  [0.1382378, -0.0276671, 0.0266161, -0.0074060, 0.02801610],
                  [0.1651957, -0.0049124, 0.0269665, -0.0040025, 0.02138184]]
        # Get iteration points
        grad_list = [prob.grad(x) for x in x_list]
        delta_x = [np.array(x_list[i+1])-np.array(x_list[i])
                   for i in range(len(x_list)-1)]
        delta_grad = [grad_list[i+1]-grad_list[i]
                      for i in range(len(grad_list)-1)]
        hess = BFGS(init_scale=1, min_curvature=10)
        hess.initialize(len(x_list[0]), 'hess')
        # Compare the Hessian and its inverse
        for i in range(len(delta_x)-1):
            s = delta_x[i]
            y = delta_grad[i]
            hess.update(s, y)
        # Test skip update
        B = np.copy(hess.get_matrix())
        s = delta_x[5]
        y = delta_grad[5]
        hess.update(s, y)
        B_updated = np.copy(hess.get_matrix())
        assert_array_equal(B, B_updated)


@pytest.mark.parametrize('strategy', [BFGS, SR1, LBFGS])
@pytest.mark.parametrize('approx_type', ['hess', 'inv_hess'])
def test_matmul_equals_dot(strategy, approx_type):
    H = strategy(init_scale=1)
    H.initialize(2, approx_type)
    v = np.array([1, 2])
    assert_array_equal(H @ v, H.dot(v))


class TestLBFGS:
    """Tests for the limited-memory BFGS Hessian update strategy."""

    @pytest.mark.parametrize('memory', [1, 3, 10])
    @pytest.mark.parametrize('approx_type', ['hess', 'inv_hess'])
    @pytest.mark.parametrize('init_scale', [2.0, 'auto'])
    def test_dense_chronological_replay(self, memory, approx_type, init_scale):
        rng = np.random.default_rng(42)
        n = 5
        a = self._spd(n, 23)
        qn = LBFGS(memory=memory, init_scale=init_scale)
        qn.initialize(n, approx_type)
        pairs = []
        auto_scale = None
        for _ in range(memory + 5):
            s = rng.normal(size=n)
            y = a @ s
            qn.update(s, y)
            pairs.append((s, y))
            pairs = pairs[-memory:]
            if init_scale == 'auto':
                if auto_scale is None:
                    auto_scale = (y @ y) / (s @ y)
                scale = auto_scale
                if approx_type == 'inv_hess':
                    scale = 1 / scale
            else:
                scale = init_scale
            expected = scale * np.eye(n)
            for step, grad in pairs:
                curvature = step @ grad
                if approx_type == 'hess':
                    bs = expected @ step
                    expected += (np.outer(grad, grad) / curvature
                                 - np.outer(bs, bs) / (step @ bs))
                else:
                    v = np.eye(n) - np.outer(step, grad) / curvature
                    expected = (v @ expected @ v.T
                                + np.outer(step, step) / curvature)
            actual = qn.get_matrix()
            assert_allclose(actual, expected, rtol=2e-11, atol=2e-11)
            assert_allclose(actual, actual.T, rtol=1e-12, atol=1e-12)
            assert np.linalg.eigvalsh(actual).min() > 0
            argument, target = (s, y) if approx_type == 'hess' else (y, s)
            assert_allclose(qn.dot(argument), target, rtol=2e-11, atol=2e-11)
            v = rng.normal(size=n)
            assert_allclose(qn.dot(v), expected @ v, rtol=2e-11, atol=2e-11)

    @pytest.mark.parametrize('memory', [1, 3, 10])
    def test_inverse_consistency_after_eviction(self, memory):
        hess, inv = LBFGS(memory=memory), LBFGS(memory=memory)
        hess.initialize(5, 'hess')
        inv.initialize(5, 'inv_hess')
        rng = np.random.default_rng(16)
        a = self._spd(5, 21)
        for _ in range(memory + 5):
            s = rng.normal(size=5)
            for qn in (hess, inv):
                qn.update(s, a @ s)
            assert_allclose(hess.get_matrix() @ inv.get_matrix(), np.eye(5),
                            rtol=2e-10, atol=2e-10)

    @pytest.mark.parametrize('approx_type', ['hess', 'inv_hess'])
    def test_copy_reinitialize_and_skipped_pairs(self, approx_type):
        qn = LBFGS(memory=3, init_scale=2)
        control = LBFGS(memory=3, init_scale=2)
        for obj in (qn, control):
            obj.initialize(5, approx_type)
        rng = np.random.default_rng(20)
        a = self._spd(5, 12)
        for _ in range(4):
            s = rng.normal(size=5)
            y = a @ s
            control.update(s.copy(), y.copy())
            qn.update(s, y)
            s[:] = 999
            y[:] = -999
        before = qn.get_matrix().copy()
        assert_allclose(before, control.get_matrix(), rtol=1e-12, atol=1e-12)
        s = rng.normal(size=5)
        qn.update(s, -s)
        qn.update(np.zeros(5), s)
        with pytest.warns(UserWarning, match='delta_grad == 0.0'):
            qn.update(s, np.zeros(5))
        assert_array_equal(qn.get_matrix(), before)
        # A subsequent accepted pair detects eviction by a skipped pair.
        for obj in (qn, control):
            obj.update(s, a @ s)
        assert_allclose(qn.get_matrix(), control.get_matrix(),
                        rtol=1e-12, atol=1e-12)
        qn.initialize(3, approx_type)
        fresh = LBFGS(memory=3, init_scale=2)
        fresh.initialize(3, approx_type)
        assert_array_equal(qn.get_matrix(), fresh.get_matrix())
        for obj in (qn, fresh):
            obj.update(np.ones(3), np.arange(1., 4.))
        assert_allclose(qn.get_matrix(), fresh.get_matrix(),
                        rtol=1e-12, atol=1e-12)

    def test_inverse_does_not_compute_w(self, monkeypatch):
        qn = LBFGS(memory=3)
        qn.initialize(5, 'inv_hess')

        def forbidden(*args, **kwargs):
            raise AssertionError('inverse mode computed Hessian products')

        monkeypatch.setattr(qn, '_apply', forbidden)
        rng = np.random.default_rng(3)
        for _ in range(8):
            s = rng.normal(size=5)
            qn.update(s, 2 * s)
            assert np.isfinite(qn.dot(s)).all()
        assert len(qn._w) == 0

    @pytest.mark.parametrize('failure', ['new_factor', 'rebuild'])
    def test_recovery_keeps_current_pair(self, monkeypatch, failure):
        qn = LBFGS(memory=2, init_scale=2)
        qn.initialize(3, 'hess')
        qn.update(np.array([1., 0., 0.]), np.array([3., 1., 0.]))
        if failure == 'rebuild':
            qn.update(np.array([0., 1., 0.]), np.array([1., 4., 0.]))
            monkeypatch.setattr(qn, '_rebuild', lambda: False)
        else:
            apply = qn._apply

            def failed_factor(s, y, w, vector):
                if s:
                    return -vector
                return apply(s, y, w, vector)

            monkeypatch.setattr(qn, '_apply', failed_factor)
        s = np.array([1., 1., 1.])
        y = np.array([3., 4., 5.])
        qn.update(s, y)
        fresh = LBFGS(init_scale=2)
        fresh.initialize(3, 'hess')
        fresh.update(s, y)
        # Compare public results after removing the injected failure.
        monkeypatch.undo()
        assert_allclose(qn.get_matrix(), fresh.get_matrix(),
                        rtol=1e-12, atol=1e-12)

    @pytest.mark.parametrize('approx_type', ['hess', 'inv_hess'])
    def test_auto_scale_is_not_refreshed(self, approx_type):
        qn = LBFGS(memory=1)
        qn.initialize(2, approx_type)
        qn.update(np.array([1., 0.]), np.array([2., 0.]))
        qn.update(np.array([0., 1.]), np.array([0., 100.]))
        expected = 2. if approx_type == 'hess' else 0.5
        assert_allclose(qn.dot(np.array([1., 0.])), [expected, 0.])

    def test_severe_cancellation_keeps_curvature(self):
        # With base I, the pairs below have exact result diag(2, 1).  A
        # dense stored matrix loses the second diagonal entry to
        # cancellation; the factored evaluation over retained pairs must
        # preserve it.
        qn = LBFGS(init_scale=1)
        qn.initialize(2, 'hess')
        qn.update(np.array([1., 0.]), np.array([1., 1e8]))
        qn.update(np.array([1., 0.]), np.array([2., 0.]))
        assert_allclose(qn.get_matrix(), np.diag([2., 1.]),
                        rtol=1e-12, atol=1e-12)
        # A subsequent pair in the weakened direction is accepted and
        # satisfies its secant equation.
        s = np.array([0., 1.])
        y = np.array([0., 3.])
        qn.update(s, y)
        assert_allclose(qn.dot(s), y, rtol=1e-12, atol=1e-12)
        assert np.linalg.eigvalsh(qn.get_matrix()).min() > 0

    @pytest.mark.parametrize('memory', [0, -1, 1.5, True, '3', np.nan, np.inf])
    def test_invalid_memory(self, memory):
        with pytest.raises((TypeError, ValueError)):
            LBFGS(memory=memory)

    @pytest.mark.parametrize('scale', [0, -1, np.nan, np.inf, -np.inf,
                                       1j, 'invalid', np.eye(2)])
    def test_invalid_scale(self, scale):
        with pytest.raises((TypeError, ValueError)):
            qn = LBFGS(init_scale=scale)
            qn.initialize(5, 'hess')
            qn.update(np.ones(5), np.ones(5))

    @pytest.mark.parametrize('n', [0, -1, 1.5, True, '5'])
    def test_invalid_dimension(self, n):
        with pytest.raises((TypeError, ValueError)):
            LBFGS().initialize(n, 'hess')

    def test_invalid_approximation(self):
        with pytest.raises(ValueError):
            LBFGS().initialize(5, 'invalid')

    @pytest.mark.parametrize('bad', [np.ones(4), np.ones((5, 1)), 1.,
                                     np.full(5, np.nan), np.full(5, np.inf),
                                     np.ones(5, dtype=complex) * (1 + 1j)])
    @pytest.mark.parametrize('operation', ['step', 'gradient', 'dot'])
    @pytest.mark.parametrize('approx_type', ['hess', 'inv_hess'])
    def test_invalid_vectors(self, bad, operation, approx_type):
        qn = LBFGS()
        qn.initialize(5, approx_type)
        if operation == 'dot' and np.shape(bad) == (5, 1):
            # SciPy's trust-region solvers also pass column vectors.
            assert qn.dot(bad).shape == (5, 1)
            return
        with pytest.raises((TypeError, ValueError)):
            if operation == 'dot':
                qn.dot(bad)
            elif operation == 'step':
                qn.update(bad, np.ones(5))
            else:
                qn.update(np.ones(5), bad)

    def _spd(self, n, seed):
        rng = np.random.default_rng(seed)
        a = rng.standard_normal((n, n))
        return a @ a.T + n * np.eye(n)

    def test_secant_condition(self):
        # The approximation has to satisfy B s = y for the pair it was
        # updated with.
        n = 12
        hessian = self._spd(n, 3)
        hess = LBFGS()
        hess.initialize(n, 'hess')
        rng = np.random.default_rng(7)
        for _ in range(8):
            s = rng.standard_normal(n)
            y = hessian @ s
            hess.update(s, y)
            assert_allclose(hess.dot(s), y, rtol=1e-11, atol=1e-11)

    def test_symmetric_and_positive_definite(self):
        n = 10
        hessian = self._spd(n, 5)
        hess = LBFGS()
        hess.initialize(n, 'hess')
        rng = np.random.default_rng(11)
        for _ in range(6):
            s = rng.standard_normal(n)
            hess.update(s, hessian @ s)
        a = rng.standard_normal(n)
        b = rng.standard_normal(n)
        assert_allclose(hess.dot(a) @ b, hess.dot(b) @ a,
                        rtol=1e-11, atol=1e-11)
        assert hess.dot(a) @ a > 0

    def test_inverse_hessian(self):
        # The 'inv_hess' product is the inverse of the 'hess' product.
        n = 12
        hessian = self._spd(n, 3)
        rng = np.random.default_rng(7)
        pairs = []
        for _ in range(10):
            s = rng.standard_normal(n)
            pairs.append((s, hessian @ s))
        hess = LBFGS()
        hess.initialize(n, 'hess')
        inv = LBFGS()
        inv.initialize(n, 'inv_hess')
        for s, y in pairs:
            hess.update(s, y)
            inv.update(s, y)
        v = np.random.default_rng(1).standard_normal(n)
        assert_allclose(inv.dot(v), np.linalg.solve(hess.get_matrix(), v),
                        rtol=1e-10, atol=1e-11)

    def test_skip_update(self):
        # A pair without positive curvature, and a zero step, must not
        # change the approximation.
        n = 5
        hess = LBFGS()
        hess.initialize(n, 'hess')
        rng = np.random.default_rng(9)
        s = rng.standard_normal(n)
        hess.update(s, 2 * s)
        B = hess.get_matrix()
        hess.update(s, -s)
        hess.update(np.zeros(n), s)
        assert_array_equal(hess.get_matrix(), B)

    def test_scaled_to_the_problem(self):
        # The initial approximation is scaled to the problem, so that a
        # caller does not have to rescale it by hand.
        n = 6
        scale = 1e4
        hess = LBFGS()
        hess.initialize(n, 'hess')
        rng = np.random.default_rng(21)
        s = rng.standard_normal(n)
        hess.update(s, scale * s)
        p = rng.standard_normal(n)
        p = p - s * (p @ s) / (s @ s)
        ratio = (hess.dot(p) @ p) / (p @ p)
        assert 0.1 * scale <= ratio <= 10 * scale

    @pytest.mark.parametrize('approx_type', ['hess', 'inv_hess'])
    def test_bounded_stored_arrays(self, approx_type):
        # Inspect owned state, not process-global allocation tracking.
        n = 50_000
        memory = 3
        hess = LBFGS(memory=memory)
        rng = np.random.default_rng(0)
        hess.initialize(n, approx_type)

        def arrays(value):
            if isinstance(value, np.ndarray):
                yield value
            elif isinstance(value, dict):
                for item in value.values():
                    yield from arrays(item)
            elif isinstance(value, (list, tuple)):
                for item in value:
                    yield from arrays(item)

        for _ in range(2 * memory + 2):
            s = rng.standard_normal(n)
            hess.update(s, s + 0.01 * rng.standard_normal(n))
            assert np.isfinite(hess.dot(s)).all()
            stored = list(arrays(vars(hess)))
            assert all(a.size <= memory * n for a in stored)
            assert sum(a.size for a in stored) <= (4 * memory + 4) * n

    def test_get_matrix(self):
        # The dense matrix must be the operator that dot() applies.
        n = 8
        hessian = self._spd(n, 17)
        hess = LBFGS()
        hess.initialize(n, 'hess')
        rng = np.random.default_rng(19)
        for _ in range(5):
            s = rng.standard_normal(n)
            hess.update(s, hessian @ s)

        M = hess.get_matrix()
        assert_allclose(M, hess.get_matrix(), rtol=1e-12, atol=1e-12)
        assert M.shape == (n, n)

        v = rng.standard_normal(n)
        assert_allclose(M @ v, hess.dot(v), rtol=1e-11, atol=1e-11)

        inv = LBFGS()
        inv.initialize(n, 'inv_hess')
        for _ in range(5):
            s = rng.standard_normal(n)
            inv.update(s, hessian @ s)

        H = inv.get_matrix()
        assert H.shape == (n, n)
        assert_allclose(H @ v, inv.dot(v), rtol=1e-11, atol=1e-11)

    def test_minimize(self):
        # The strategy can be used as the 'hess' argument of minimize.
        n = 20
        hessian = self._spd(n, 13)
        linear = np.random.default_rng(6).standard_normal(n)
        x0 = np.random.default_rng(8).standard_normal(n)

        def fun(x):
            return 0.5 * x @ hessian @ x + linear @ x

        def jac(x):
            return hessian @ x + linear

        result = minimize(fun, x0, jac=jac, hess=LBFGS(),
                          method='trust-constr',
                          options={'maxiter': 300, 'gtol': 1e-8})
        optimum = np.linalg.solve(hessian, -linear)
        expected = 0.5 * optimum @ hessian @ optimum + linear @ optimum
        assert result.success
        assert abs(result.fun - expected) <= 1e-6 * max(1.0, abs(expected))
