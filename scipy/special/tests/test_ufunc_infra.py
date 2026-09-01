import inspect
import numpy as np
import pickle
import pytest
import warnings

from numpy.testing import assert_equal
from numpy._core._exceptions import UFuncTypeError
from scipy._external.packaging_version import version
# ufunc wrapper with array api dispatch
from scipy.special import mathieu_sem
# ufunc using special._ufunc_tools._with_cache_optimization
from scipy.special._mathieu import mathieu_sem as mathieu_sem_wrapper
# raw ufunc without cache optimization
from scipy.special._ufuncs import _mathieu_sem

from scipy.special._spfun_stats import _poisson_binom_cdf
from scipy.special._ufuncs import betainc

from scipy.special._ufunc_tools import _make_ufunc_wrapper


# Tests that ufunc kwargs still work when _with_cache_optimization is applied
class TestWithCacheOptimization:
    @pytest.mark.parametrize("m_shape,q_shape,x_shape,where_shape", [
        ((5, 1), (1, 1), (1, 10), ()),
        ((2, 3), (2, 3), (2, 3), ()),
        ((1, 5), (1, 5), (10, 5), ()),
        ((5, 1), (1, 1), (1, 10), (1, 10)),
        ((2, 3), (2, 3), (2, 3), (2, 3)),
        ((1, 5), (1, 5), (10, 5), (10, 5)),
        ((5, 1), (1, 1), (1, 10), (5, 1)),
        ((2, 3), (2, 3), (2, 3), (2, 1)),
        ((1, 5), (1, 5), (10, 5), (1, 5)),
        ((5, 1), (5, 1), (5, 1), (3, 1, 1)),
        ((1, 10), (1, 10), (1, 10), (5, 1)),
        ((5,), (5,), (5,), (2, 3, 1, 1)),
        ((5, 1), (5, 1), (1, 1), (1, 10)),
        ((5, 1, 1), (1, 4, 1), (1, 1, 3), (2, 1, 1, 1)),
    ])
    @pytest.mark.parametrize("out_order", ["C", "F"])
    def test_out(self, m_shape, q_shape, x_shape, where_shape, out_order):
        # the call to mathieu_sem in this test raised a segfault on scipy/main prior to
        # https://github.com/scipy/scipy/pull/25127 -- https://github.com/scipy/xsf/pull/146.
        rng = np.random.default_rng(1234)

        m = rng.integers(1, 20, m_shape)
        q = rng.uniform(0, 10, q_shape)
        x = rng.uniform(0, 90, x_shape)

        if where_shape == ():
            batch_shape = np.broadcast_shapes(m_shape, q_shape, x_shape)
            where = True
        else:
            where = rng.choice([True, False], size=where_shape)
            batch_shape = np.broadcast_shapes(m_shape, q_shape, x_shape, where_shape)

        out0 = np.full(batch_shape, np.nan, order=out_order)
        out1 = np.full(batch_shape, np.nan, order=out_order)
        res0, res1 = mathieu_sem(m, q, x, out=(out0, out1), where=where)
        assert res0 is out0 and res1 is out1

        expected0 = np.full(batch_shape, np.nan, order=out_order)
        expected1 = np.full(batch_shape, np.nan, order=out_order)
        _mathieu_sem(m, q, x, out=(expected0, expected1), where=where)
        assert_equal((res0, res1), (expected0, expected1))
        assert res0.flags == expected0.flags
        assert res1.flags == expected1.flags

    def test_out_broadcasts(self):
        m = [[1, 3, 4]]
        q = [[2, 2, 2]]
        x = [[30, 60, 90]]
        out0 = np.full((3, 3), np.nan)
        out1 = np.full((3, 3), np.nan)
        mathieu_sem(m, q, x, out=(out0, out1))
        expected0 = np.full((3, 3), np.nan)
        expected1 = np.full((3, 3), np.nan)
        mathieu_sem(m, q, x, out=(expected0, expected1))
        assert_equal((out0, out1), (expected0, expected1))

    @pytest.mark.parametrize("subok", [True, False])
    def test_subok(self, subok):
        class ArraySubClass(np.ndarray):
            pass
        m = np.asarray([1, 2, 3]).view(ArraySubClass)
        q = np.asarray([2.1, 3.2, 4.3]).view(ArraySubClass)
        x = np.asarray([10, 20, 30]).view(ArraySubClass)
        res0, res1 = mathieu_sem(m, q, x, subok=subok)
        expected0, expected1 = _mathieu_sem(m, q, x, subok=subok)
        assert type(res0) is type(expected0) and type(res1) is type(expected1)

    @pytest.mark.parametrize("order_m", ["C", "F"])
    @pytest.mark.parametrize("order_q", ["C", "F"])
    @pytest.mark.parametrize("order_x", ["C", "F"])
    # order="K" works differently with _with_cache_optimization
    @pytest.mark.parametrize("order", ["C", "F", "A"])
    @pytest.mark.parametrize("m_shape,q_shape,x_shape", [
        ((5, 1), (1, 1), (1, 10)),
        ((2, 3), (2, 3), (2, 3)),
        ((1, 5), (1, 5), (10, 5)),
        ((5, 1), (1, 1), (1, 10)),
        ((2, 3), (2, 3), (2, 3)),
        ((1, 5), (1, 5), (10, 5)),
        ((5, 1), (1, 1), (1, 10)),
        ((2, 3), (2, 3), (2, 3)),
        ((1, 5), (1, 5), (10, 5)),
        ((5, 1), (5, 1), (5, 1)),
        ((1, 10), (1, 10), (1, 10)),
        ((5,), (5,), (5,)),
        ((5, 1), (5, 1), (1, 1)),
        ((5, 1, 1), (1, 4, 1), (1, 1, 3)),
    ])
    def test_order(self,order_m, order_q, order_x, order, m_shape, q_shape, x_shape):
        rng = np.random.default_rng(1234)
        m = np.asarray(rng.integers(1, 20, m_shape), copy=True, order=order_m)
        q = np.asarray(rng.uniform(0, 10, q_shape), copy=True, order=order_q)
        x = np.asarray(rng.uniform(0, 90, x_shape), copy=True, order=order_x)
        res0, res1 = mathieu_sem(m, q, x, order=order)
        expected0, expected1 = _mathieu_sem(m, q, x, order=order)
        assert_equal((res0, res1), (expected0, expected1))

        assert res0.flags == expected0.flags
        assert res1.flags == expected1.flags

    def test_casting(self):
        m = np.float16(2.0)
        q = np.float16(1.0)
        x = np.float16(45.0)
        with pytest.raises(UFuncTypeError):
            res0, res1 = mathieu_sem(m, q, x, casting="no")

        res0, res1 = mathieu_sem(m, q, x, casting="safe")
        expected0, expected1 = _mathieu_sem(m, q, x, casting="safe")
        assert_equal((res0, res1), (expected0, expected1))

    @pytest.mark.parametrize("dtype_kwarg", [np.float32, np.float64])
    @pytest.mark.parametrize("input_dtype", [np.float16, np.float32, np.float64])
    def test_dtype(self, dtype_kwarg, input_dtype):
        m = np.asarray([1, 2, 3], dtype=input_dtype)
        q = np.asarray([2, 3, 4], dtype=input_dtype)
        x = np.asarray([30, 60, 90], dtype=input_dtype)

        res0, res1 = mathieu_sem(m, q, x, dtype=dtype_kwarg)
        expected0, expected1 = _mathieu_sem(m, q, x, dtype=dtype_kwarg)
        assert_equal((res0, res1), (expected0, expected1))
        assert res0.dtype == res1.dtype == expected0.dtype == expected1.dtype

    @pytest.mark.parametrize(
        "signature", [
            (np.float32, np.float32, np.float32, np.float32, np.float32),
            (np.float64, np.float64, np.float64, np.float64, np.float64),
        ],
    )
    def test_signature(self, signature):
        m = [1, 2, 3]
        q = [2, 3, 4]
        x = [30, 60, 90]
        res0, res1 = mathieu_sem(m, q, x, signature=signature)
        expected0, expected1 = _mathieu_sem(m, q, x, signature=signature)
        assert_equal((res0, res1), (expected0, expected1))
        assert res0.dtype == res1.dtype == expected0.dtype == expected1.dtype

    @pytest.mark.parametrize("m", [2.0, [2.0, 2.0, 2.0]])
    @pytest.mark.parametrize("where", [[True, False, True], False])
    def test_where_without_out(self, m, where):
        # ``where`` must be honored whether or not the arguments participating
        # in the cache have batch dimensions, i.e. on both the fast path (m
        # scalar) and the transposed path (m batched). gh-25127
        q = 3.0
        x = [10.0, 20.0, 30.0]
        if version.parse(np.__version__) >= version.parse("2.4.0"):
            with pytest.warns(UserWarning, match="'where' used without 'out'"):
                res0, res1 = mathieu_sem(m, q, x, where=where)
        else:
            res0, res1 = mathieu_sem(m, q, x, where=where)

        mask = np.broadcast_to(where, (3,))
        expected0, expected1 = mathieu_sem(m, q, x)
        assert_equal(res0[mask], expected0[mask])
        assert_equal(res1[mask], expected1[mask])

    @pytest.mark.parametrize("m", [2.0, [2.0, 2.0, 2.0]])
    @pytest.mark.parametrize("where", [[True, False, True], False])
    def test_where_with_explicit_out_none(self, m, where):
        q = 3.0
        x = [10.0, 20.0, 30.0]

        with warnings.catch_warnings():
            warnings.filterwarnings(
                "error",
                message="'where' used without 'out'",
                category=UserWarning,
            )
            res0, res1 = mathieu_sem(m, q, x, where=where, out=(None, None))

        mask = np.broadcast_to(where, (3,))
        expected0, expected1 = mathieu_sem(m, q, x)
        assert_equal(res0[mask], expected0[mask])
        assert_equal(res1[mask], expected1[mask])

    def test_where_broadcasts_like_ufunc(self):
        rng = np.random.default_rng(1234)
        m = rng.integers(1, 20, (5, 1))
        q = rng.uniform(0, 10, (1, 4))
        x = rng.uniform(0, 90, (5, 4))
        where = rng.choice([True, False], size=(5, 4))

        out = [np.full((5, 4), np.nan) for _ in range(2)]
        expected = [np.full((5, 4), np.nan) for _ in range(2)]
        mathieu_sem(m, q, x, out=tuple(out), where=where)
        _mathieu_sem(m, q, x, out=tuple(expected), where=where)
        assert_equal(out, expected)

    def test_out_with_none_entry(self):
        # A None entry means "allocate this output for me", as for a plain
        # ufunc.
        m = [[1, 3, 4]]
        q = [[2, 2, 2]]
        x = [[30, 60, 90]]
        out0 = np.full((3, 3), np.nan)
        expected0 = np.full((3, 3), np.nan)
        res0, res1 = mathieu_sem(m, q, x, out=(out0, None))
        exp0, exp1 = _mathieu_sem(m, q, x, out=(expected0, None))
        assert res0 is out0
        assert_equal((res0, res1), (exp0, exp1))

    def test_out_invalid(self):
        m = [[1, 3, 4]]
        q = [[2, 2, 2]]
        x = [[30, 60, 90]]
        out0 = np.full((3, 3), np.nan)
        # mathieu_sem has two outputs, so a bare array is not acceptable.
        with pytest.raises(TypeError, match="must be a tuple of arrays"):
            mathieu_sem(m, q, x, out=out0)
        with pytest.raises(TypeError, match="must be a tuple of arrays"):
            mathieu_sem(m, q, x, out=[out0, out0.copy()])
        with pytest.raises(ValueError, match="exactly 2 entries"):
            mathieu_sem(m, q, x, out=(out0,))

    def test_unexpected_kwarg_names_wrapper(self):
        with pytest.raises(TypeError, match="mathieu_sem"):
            mathieu_sem(1, 2, 3, not_a_ufunc_kwarg=True)

    def test_wrapper_metadata(self):
        for attr in ["nin", "nout", "nargs", "ntypes", "types", "signature"]:
            assert getattr(mathieu_sem, attr) == getattr(_mathieu_sem, attr)

    @pytest.mark.parametrize("func", [mathieu_sem, mathieu_sem_wrapper])
    def test_pickle(self, func):
        assert pickle.loads(pickle.dumps(func)) is func


def _make_passthrough(ufunc):
    def wrapper(*args, **kwargs):
        return ufunc(*args, **kwargs)
    return wrapper

_betainc_wrapper = _make_ufunc_wrapper(
            _make_passthrough(betainc),
            betainc,
            "betainc",
            ["a", "b", "x"],
            "Wrapper for betainc.",
)

_poisson_binom_cdf_wrapper = _make_ufunc_wrapper(
            _make_passthrough(_poisson_binom_cdf),
            _poisson_binom_cdf,
            "poisson_binom_cdf",
            ["k", "p"],
            "Wrapper for _poisson_binom_cdf.",
)

_mathieu_sem_wrapper_wrapper = _make_ufunc_wrapper(
    _make_passthrough(mathieu_sem_wrapper),
    mathieu_sem_wrapper,
    "mathieu_sem",
    ["m", "q", "x"],
    "Wrapper for wrapper of mathieu_sem",
)

_vecdot_wrapper = _make_ufunc_wrapper(
    _make_passthrough(np.vecdot),
    np.vecdot,
    "vecdot",
    ["x1", "x2"],
    "Wrapper for vecdot."
)


def _assert_same_result(actual, desired):
    if isinstance(desired, tuple):
        assert isinstance(actual, tuple)
        assert len(actual) == len(desired)
        for actual_i, desired_i in zip(actual, desired):
            _assert_same_result(actual_i, desired_i)
        return

    np.testing.assert_equal(actual, desired)

    if isinstance(desired, np.ndarray):
        assert type(actual) is type(desired)
        assert actual.dtype == desired.dtype
        assert actual.shape == desired.shape
        assert actual.strides == desired.strides


class TestMakeUFuncWrapper:
    @pytest.mark.parametrize(
        "func,func_wrapper",
        [
            [betainc, _betainc_wrapper],
            [mathieu_sem_wrapper, _mathieu_sem_wrapper_wrapper],
            [_poisson_binom_cdf, _poisson_binom_cdf_wrapper],
            [np.vecdot, _vecdot_wrapper],
        ]
    )
    def test_attrs(self, func, func_wrapper):
        for attr in ["nin", "nout", "nargs", "ntypes", "types", "signature"]:
            assert getattr(func_wrapper, attr) == getattr(func, attr)

        # check that mutating the returned types list does not mutate the
        # wrapper state.
        types = func_wrapper.types
        types[0] = None
        assert func_wrapper.types == func.types

    def test_simple_elementwise_values(self):
        a = np.linspace(1, 100, 9)
        b = np.linspace(1, 100, 8)
        x = np.linspace(0, 1, 10)
        desired = betainc(a[:, None], b[None, :], x[:, None, None])
        actual = _betainc_wrapper(a[:, None], b[None, :], x[:, None, None])
        np.testing.assert_equal(actual, desired)

    def test_simple_elementwise_values_out(self):
        a = np.linspace(1, 100, 9)
        b = np.linspace(1, 100, 8)
        x = np.linspace(0, 1, 10)

        out_desired = np.empty((10, 9, 8))
        betainc(a[:, None], b[None, :], x[:, None, None], out=out_desired)
        out_actual = np.empty((10, 9, 8))
        actual = _betainc_wrapper(
            a[:, None], b[None, :], x[:, None, None], out=out_actual
        )
        assert out_actual is actual
        np.testing.assert_equal(out_actual, out_desired)

    def test_multiple_outputs_elementwise_values(self):
        m = np.asarray([1, 4])
        q = np.linspace(1, 50, 10)
        x = np.linspace(0, 360, 10)
        desired0, desired1 = mathieu_sem_wrapper(
            m[:, None], q[None, :], x[:, None, None]
        )
        actual0, actual1 = _mathieu_sem_wrapper_wrapper(
            m[:, None], q[None, :], x[:, None, None]
        )
        np.testing.assert_equal(actual0, desired0)
        np.testing.assert_equal(actual1, desired1)

    def test_multiple_outputs_values_out(self):
        m = np.asarray([1, 4])
        q = np.linspace(1, 50, 10)
        x = np.linspace(0, 360, 10)
        out_desired0, out_desired1 = np.empty((10, 2, 10)), np.empty((10, 2, 10))
        desired0, desired1 = mathieu_sem_wrapper(
            m[:, None], q[None, :], x[:, None, None], out=(out_desired0, out_desired1)
        )
        out_actual0, out_actual1 = np.empty((10, 2, 10)), np.empty((10, 2, 10))
        actual0, actual1 = _mathieu_sem_wrapper_wrapper(
            m[:, None], q[None, :], x[:, None, None], out=(out_actual0, out_actual1)
        )
        assert out_actual0 is actual0
        assert out_actual1 is actual1
        np.testing.assert_equal(out_actual0, out_desired0)
        np.testing.assert_equal(out_actual1, out_desired1)

    def test_gufunc_values_out(self):
        p = np.asarray([
            [0.1, 0.2, 0.3],
            [0.4, 0.5, 0.6],
        ])
        k = np.asarray([
            [0, 1, 2],
            [2, 0, 1],
        ])

        out_desired = np.empty((2, 3, 2))
        _poisson_binom_cdf(k[:, :, None], p, out=out_desired)
        out_actual = np.empty((2, 3, 2))
        actual = _poisson_binom_cdf_wrapper(k[:, :, None], p, out=out_actual)
        assert out_actual is actual
        np.testing.assert_equal(out_actual, out_desired)

    @pytest.mark.parametrize(
        "a_dtype", [np.dtype("int32"), np.dtype("float32"), np.dtype("float64")]
    )
    @pytest.mark.parametrize(
        "b_dtype", [np.dtype("int32"), np.dtype("float32"), np.dtype("float64")]
    )
    @pytest.mark.parametrize("x_dtype", [np.dtype("float32"), np.dtype("float64")])
    @pytest.mark.parametrize(
        "out_dtype", [None, np.dtype("float32"), np.dtype("float64")]
    )
    def test_resolve_dtypes(self, a_dtype, b_dtype, x_dtype, out_dtype):
        desired = betainc.resolve_dtypes((a_dtype, b_dtype, x_dtype, out_dtype))
        actual = _betainc_wrapper.resolve_dtypes((a_dtype, b_dtype, x_dtype, out_dtype))
        assert actual == desired

    def test_elementwise_where(self):
        a = np.array([1., 2., 3.])
        b = np.array([2., 3., 4.])
        x = np.array([0.2, 0.4, 0.6])
        where = np.array([True, False, True])

        actual = np.full(3, -1.)
        desired = actual.copy()

        _betainc_wrapper(a, b, x, out=actual, where=where)
        betainc(a, b, x, out=desired, where=where)

        np.testing.assert_equal(actual, desired)

    @pytest.mark.parametrize(
        "kwargs",
        [
            {"casting": "same_kind"},
            {"casting": "unsafe"},
            {"order": "C"},
            {"order": "F"},
            {"order": "A"},
            {"order": "K"},
            {"dtype": np.float32},
            {"dtype": np.float64},
            {"signature": betainc.types[0]},
            {"signature": betainc.types[-1]},
        ],
    )
    def test_kwargs(self, kwargs):
        a = np.linspace(1, 100, 9)[:, None, None]
        b = np.linspace(1, 100, 8)[None, :, None]
        x = np.linspace(0, 1, 10)[None, None, :]

        desired = betainc(a, b, x, **kwargs)
        actual = _betainc_wrapper(a, b, x, **kwargs)

        _assert_same_result(actual, desired)


    @pytest.mark.parametrize("subok", [True, False])
    def test_subok(self, subok):
        class ArraySubclass(np.ndarray):
            pass

        a = np.linspace(1, 100, 9).view(ArraySubclass)
        b = np.linspace(1, 100, 9).view(ArraySubclass)
        x = np.linspace(0, 1, 9).view(ArraySubclass)

        desired = betainc(a, b, x, subok=subok)
        actual = _betainc_wrapper(a, b, x, subok=subok)

        _assert_same_result(actual, desired)

    @pytest.mark.parametrize(
        "kwargs",
        [
            {"dtype": np.float32, "casting": "no"},
            {"signature": "fff->f", "casting": "no"},
        ],
    )
    def test_kwargs_errors(self, kwargs):
        a = np.array([1., 2., 3.], dtype=np.float64)
        b = np.array([2., 3., 4.], dtype=np.float64)
        x = np.array([0.2, 0.4, 0.6], dtype=np.float64)

        with pytest.raises(TypeError):
            betainc(a, b, x, **kwargs)

        with pytest.raises(TypeError):
            _betainc_wrapper(a, b, x, **kwargs)

    @pytest.mark.parametrize("axis", [0, 1, 2, -1])
    def test_gufunc_axis(self, axis):
        x1 = np.arange(24.0).reshape(2, 3, 4)
        x2 = np.arange(24.0).reshape(2, 3, 4)

        desired = np.vecdot(x1, x2, axis=axis)
        actual = _vecdot_wrapper(x1, x2, axis=axis)

        _assert_same_result(actual, desired)

    @pytest.mark.parametrize(
        "axes",
        [
            [(0,), (0,), ()],
            [(1,), (1,), ()],
            [(2,), (2,), ()],
        ],
    )
    def test_gufunc_axes(self, axes):
        x1 = np.arange(24.0).reshape(2, 3, 4)
        x2 = np.arange(24.0).reshape(2, 3, 4)

        desired = np.vecdot(x1, x2, axes=axes)
        actual = _vecdot_wrapper(x1, x2, axes=axes)

        _assert_same_result(actual, desired)

    @pytest.mark.parametrize("keepdims", [True, False])
    def test_gufunc_keepdims(self, keepdims):
        x1 = np.arange(24.0).reshape(2, 3, 4)
        x2 = np.arange(24.0).reshape(2, 3, 4)

        desired = np.vecdot(x1, x2, keepdims=keepdims)
        actual = _vecdot_wrapper(x1, x2, keepdims=keepdims)

        _assert_same_result(actual, desired)

    def test_gufunc_axes_keepdims(self):
        x1 = np.arange(24.0).reshape(2, 3, 4)
        x2 = np.arange(24.0).reshape(2, 3, 4)

        axes = [(1,), (1,), (1,)]

        desired = np.vecdot(x1, x2, axes=axes, keepdims=True)
        actual = _vecdot_wrapper(x1, x2, axes=axes, keepdims=True)

        _assert_same_result(actual, desired)

    @pytest.mark.parametrize(
        "func_wrapper, expected",
        [
            (_betainc_wrapper, "(a, b, x, /, out=None, **kwargs)"),
            (_mathieu_sem_wrapper_wrapper, "(m, q, x, /, out=None, **kwargs)"),
            (_poisson_binom_cdf_wrapper, "(k, p, /, out=None, **kwargs)"),
            (_vecdot_wrapper, "(x1, x2, /, out=None, **kwargs)"),
        ],
    )
    def test_call_signature(self, func_wrapper, expected):
        assert str(inspect.signature(func_wrapper)) == expected

    def test_out_none_vs_not_supplied(self):
        def func(*args, **kwargs):
            return kwargs

        wrapper = _make_ufunc_wrapper(
            func,
            betainc,
            "test_func",
            ["a", "b", "x"],
            "Test wrapper.",
        )

        kwargs = wrapper(1.0, 2.0, 0.5)
        assert "out" not in kwargs

        kwargs = wrapper(1.0, 2.0, 0.5, out=None)
        assert "out" in kwargs
        assert kwargs["out"] is None

        kwargs = wrapper(1.0, 2.0, 0.5, None)
        assert "out" in kwargs
        assert kwargs["out"] is None

    @pytest.mark.parametrize("kwarg", ["signature", "casting"])
    def test_resolve_dtypes_explicit_none(self, kwarg):
        dtypes = (
            np.dtype("float64"),
            np.dtype("float64"),
            np.dtype("float64"),
            None,
        )
        with pytest.raises(TypeError):
            _betainc_wrapper.resolve_dtypes(dtypes, **{kwarg: None})
