import warnings

import numpy as np
import pytest

from scipy.interpolate import LSQMultivariateSpline


def test_exact_quadratic_1d():
    x = np.linspace(-1.0, 1.0, 25)
    y = 1.5 * x**2 - 0.25 * x + 2.0

    spline = LSQMultivariateSpline(x=x, y=y, t=1, k=2)
    x_eval = np.linspace(-0.8, 0.8, 9)

    np.testing.assert_allclose(
        spline(x_eval),
        1.5 * x_eval**2 - 0.25 * x_eval + 2.0,
        atol=1e-12,
    )


def test_explicit_full_knot_construction():
    x = np.linspace(0.0, 1.0, 20)
    y = x.copy()
    spline = LSQMultivariateSpline(x=x, y=y, t=[0.25, 0.5, 0.75], k=2)

    np.testing.assert_allclose(
        spline.get_knots()[0],
        [0.0, 0.0, 0.0, 0.25, 0.5, 0.75, 1.0, 1.0, 1.0],
    )


def test_scalar_eval_returns_scalar_for_1d():
    x = np.linspace(0.0, 1.0, 20)
    y = 2.0 * x + 1.0
    spline = LSQMultivariateSpline(x=x, y=y, t=1, k=1)

    value = spline(0.25)

    assert np.ndim(value) == 0
    np.testing.assert_allclose(value, 1.5)


def test_coordinate_list_input_form():
    x0 = np.array([-1.0, -1.0, 1.0, 1.0])
    x1 = np.array([-1.0, 1.0, -1.0, 1.0])
    x = np.vstack((x0, x1))
    y = 1.0 + x0 - 2.0 * x1

    spline = LSQMultivariateSpline(x=x, y=y, t=[1, 1], k=1)

    np.testing.assert_allclose(spline([[0.5, -0.5]]), [2.5], atol=1e-12)


def test_exact_plane_2d_and_derivatives():
    x0_axis = np.linspace(-1.0, 1.0, 7)
    x1_axis = np.linspace(-1.5, 1.5, 8)
    x0_mesh, x1_mesh = np.meshgrid(x0_axis, x1_axis, indexing="ij")
    x = np.column_stack((x0_mesh.ravel(), x1_mesh.ravel()))
    y = 2.0 * x[:, 0] - 3.0 * x[:, 1] + 1.0

    spline = LSQMultivariateSpline(x=x, y=y, t=[1, 1], k=1)
    x_eval = np.array(
        [
            [-0.5, -0.25],
            [0.0, 0.0],
            [0.75, 1.0],
        ]
    )

    expected = 2.0 * x_eval[:, 0] - 3.0 * x_eval[:, 1] + 1.0
    np.testing.assert_allclose(spline(x_eval), expected, atol=1e-12)
    np.testing.assert_allclose(spline(x_eval, nu=[1, 0]), 2.0, atol=1e-12)
    np.testing.assert_allclose(spline(x_eval, nu=[0, 1]), -3.0, atol=1e-12)


def test_2d_derivative_matches_finite_difference():
    x0 = np.linspace(-1.0, 1.0, 12)
    x1 = np.linspace(-1.0, 1.0, 13)
    x0_mesh, x1_mesh = np.meshgrid(x0, x1, indexing="ij")
    y = np.sin(x0_mesh) + np.cos(x1_mesh)
    spline = LSQMultivariateSpline.from_grid((x0, x1), y, t=[4, 4], k=3)

    point = np.array([[0.2, -0.3]])
    step = 1e-6
    forward = spline(point + [[step, 0.0]])
    backward = spline(point - [[step, 0.0]])
    finite_difference = (forward - backward) / (2.0 * step)

    np.testing.assert_allclose(
        spline(point, nu=[1, 0]),
        finite_difference,
        rtol=1e-5,
        atol=1e-7,
    )


def test_derivative_order_above_degree_is_zero():
    x = np.linspace(0.0, 1.0, 20)
    y = 2.0 * x + 1.0
    spline = LSQMultivariateSpline(x=x, y=y, t=1, k=1)

    np.testing.assert_allclose(spline([0.25, 0.5], nu=2), 0.0)


def test_exact_linear_time_3d():
    x0_axis = np.linspace(-1.0, 1.0, 5)
    x1_axis = np.linspace(-1.0, 1.0, 6)
    time_axis = np.linspace(0.0, 4.0, 5)
    x0_mesh, x1_mesh, time_mesh = np.meshgrid(
        x0_axis,
        x1_axis,
        time_axis,
        indexing="ij",
    )
    y = 1.0 + x0_mesh + 2.0 * x1_mesh - 0.5 * time_mesh

    spline = LSQMultivariateSpline.from_grid(
        (x0_axis, x1_axis, time_axis),
        y,
        t=[1, 1, 1],
        k=1,
    )
    x_eval = np.array(
        [
            [0.2, -0.4, 1.5],
            [-0.8, 0.5, 3.0],
        ]
    )
    expected = 1.0 + x_eval[:, 0] + 2.0 * x_eval[:, 1] - 0.5 * x_eval[:, 2]

    np.testing.assert_allclose(spline(x_eval), expected, atol=1e-12)


def test_weights_are_applied_to_residual():
    x = np.linspace(0.0, 1.0, 8)
    y = np.array([1.0, 1.2, 1.6, 2.1, 2.4, 2.8, 3.3, 3.7])
    w = np.linspace(1.0, 2.0, x.size)

    spline = LSQMultivariateSpline(x=x, y=y, t=1, k=1, w=w)
    residual = np.sum((w * (spline(x) - y)) ** 2)

    np.testing.assert_allclose(spline.get_residual(), residual)


def test_explicit_bbox_flat_form():
    x = np.linspace(0.2, 0.8, 20)
    y = x.copy()

    spline = LSQMultivariateSpline(x=x, y=y, t=[0.5], bbox=[0.0, 1.0], k=1)

    np.testing.assert_allclose(spline.get_knots()[0], [0.0, 0.0, 0.5, 1.0, 1.0])


def test_automatic_knot_count_for_uniform_data():
    x = np.linspace(0.0, 10.0, 101)
    y = np.sin(x)

    spline = LSQMultivariateSpline(x=x, y=y, t=4, k=3)
    full_knots = spline.get_knots()[0]

    np.testing.assert_allclose(full_knots[4:-4], [2.5, 5.0, 7.5])


def test_matches_scipy_lsq_univariate_spline():
    interpolate = pytest.importorskip("scipy.interpolate")

    x = np.linspace(0.0, 2.0 * np.pi, 60)
    y = np.sin(x) + 0.25 * np.cos(2.0 * x)
    interior_knots = np.linspace(x.min(), x.max(), 7)[1:-1]

    ours = LSQMultivariateSpline(x=x, y=y, t=interior_knots, k=3)
    scipy_spline = interpolate.LSQUnivariateSpline(
        x,
        y,
        interior_knots,
        k=3,
    )
    x_eval = np.linspace(x.min(), x.max(), 25)

    np.testing.assert_allclose(ours(x_eval), scipy_spline(x_eval), atol=1e-10)


def test_exact_bilinear_2d():
    x0_axis = np.linspace(-1.0, 1.0, 8)
    x1_axis = np.linspace(-1.0, 1.0, 9)
    x0_mesh, x1_mesh = np.meshgrid(x0_axis, x1_axis, indexing="ij")
    z = 1.0 + 2.0 * x0_mesh - 0.5 * x1_mesh + 0.25 * x0_mesh * x1_mesh
    x = np.column_stack((x0_mesh.ravel(), x1_mesh.ravel()))
    x_eval = np.array([[-0.6, -0.4], [0.0, 0.0], [0.7, 0.5]])
    expected = (
        1.0
        + 2.0 * x_eval[:, 0]
        - 0.5 * x_eval[:, 1]
        + 0.25 * x_eval[:, 0] * x_eval[:, 1]
    )

    spline = LSQMultivariateSpline(x=x, y=z.ravel(), t=[[0.0], [0.0]], k=1)

    np.testing.assert_allclose(spline(x_eval), expected, atol=1e-12)


def test_too_few_points_for_unregularized_fit():
    x = np.array([0.0, 1.0])
    y = np.array([0.0, 1.0])

    with pytest.raises(ValueError, match="not enough data points"):
        LSQMultivariateSpline(x=x, y=y, t=[0.5], k=1)


def test_rank_deficient_design_is_rejected():
    x0 = np.linspace(0.0, 1.0, 5)
    x1 = x0.copy()
    x = np.column_stack((x0, x1))
    y = np.zeros(x0.shape)

    with pytest.raises(ValueError, match="rank deficient"):
        LSQMultivariateSpline(
            x=x,
            y=y,
            t=[1, 1],
            bbox=[[0.0, 1.0], [0.0, 1.0]],
            k=1,
        )


def test_unsupported_basis_is_rejected():
    x = np.linspace(0.0, 0.4, 6)
    y = x.copy()

    with pytest.raises(ValueError, match="support every tensor-product basis"):
        LSQMultivariateSpline(x=x, y=y, t=[0.5], bbox=[0.0, 1.0], k=1)


def test_unsupported_basis_can_warn_and_continue():
    x = np.linspace(0.0, 0.4, 6)
    y = x.copy()

    with pytest.warns(RuntimeWarning) as warning_record:
        spline = LSQMultivariateSpline(
            x=x,
            y=y,
            t=[0.5],
            bbox=[0.0, 1.0],
            k=1,
            unsupported="warn",
        )

    messages = [str(warning.message) for warning in warning_record]
    assert any("do not support" in message for message in messages)
    assert spline.get_coeffs().shape == (3,)
    np.testing.assert_allclose(spline([0.0, 0.4]), [0.0, 0.4], atol=1e-14)


def test_unsupported_basis_can_ignore_and_continue():
    x = np.linspace(0.0, 0.4, 6)
    y = x.copy()

    with warnings.catch_warnings(record=True) as warning_record:
        warnings.simplefilter("always")
        spline = LSQMultivariateSpline(
            x=x,
            y=y,
            t=[0.5],
            bbox=[0.0, 1.0],
            k=1,
            unsupported="ignore",
        )

    assert warning_record == []
    assert spline.get_coeffs().shape == (3,)


def test_high_dimensional_unregularized_fit_is_rejected():
    x = np.linspace(0.0, 1.0, 12)[:, None]
    x = np.repeat(x, 4, axis=1)
    y = np.zeros(x.shape[0])

    with pytest.raises(ValueError, match="not enough data points"):
        LSQMultivariateSpline(x=x, y=y, t=[4, 4, 4, 4], k=3)


@pytest.mark.parametrize(
    "kwargs, message",
    [
        ({"k": 0}, "spline degrees must be positive"),
        ({"t": 0}, "automatic knot counts"),
        ({"w": np.ones(4)}, "same length as y"),
        ({"bbox": [1.0, 0.0]}, "lower bound"),
        ({"smoothing": -1.0}, "nonnegative"),
        ({"penalty_order": 0}, "positive integer"),
        ({"sparse": "yes"}, "sparse"),
        ({"unsupported": "continue"}, "unsupported"),
    ],
)
def test_invalid_constructor_inputs(kwargs, message):
    x = np.linspace(0.0, 1.0, 5)
    y = x.copy()
    params = {"x": x, "y": y, "t": 1, "k": 1}
    params.update(kwargs)

    with pytest.raises(ValueError, match=message):
        LSQMultivariateSpline(**params)


@pytest.mark.parametrize(
    "kwargs, message",
    [
        ({"y": np.ones((5, 1, 1))}, "1D or 2D"),
        ({"t": [0.25, 0.25]}, "strictly increasing"),
        ({"t": [-0.1]}, "strictly inside bbox"),
        ({"bbox": [0.2, 0.8]}, "inside bbox"),
    ],
)
def test_invalid_shape_and_domain_inputs(kwargs, message):
    x = np.linspace(0.0, 1.0, 5)
    y = x.copy()
    params = {"x": x, "y": y, "t": [0.5], "k": 1}
    params.update(kwargs)

    with pytest.raises(ValueError, match=message):
        LSQMultivariateSpline(**params)


def test_check_finite_rejects_nan_values():
    x = np.linspace(0.0, 1.0, 5)
    y = x.copy()
    y[2] = np.nan

    with pytest.raises(ValueError, match="finite"):
        LSQMultivariateSpline(x=x, y=y, t=1, k=1, check_finite=True)


def test_invalid_eval_shape_and_derivative_order():
    x0 = np.linspace(0.0, 1.0, 4)
    x1 = np.linspace(0.0, 1.0, 5)
    x0_mesh, x1_mesh = np.meshgrid(x0, x1, indexing="ij")
    y = x0_mesh + x1_mesh
    spline = LSQMultivariateSpline.from_grid((x0, x1), y, t=[1, 1], k=1)

    with pytest.raises(ValueError, match="shape"):
        spline([0.1, 0.2, 0.3])

    with pytest.raises(ValueError, match="nonnegative"):
        spline([[0.1, 0.2]], nu=[-1, 0])

    with pytest.raises(ValueError, match="scalar evaluation"):
        spline(0.1)


def test_evaluation_outside_bbox_is_rejected():
    x = np.linspace(0.0, 1.0, 10)
    y = x.copy()
    spline = LSQMultivariateSpline(x=x, y=y, t=1, k=1)

    with pytest.raises(ValueError, match="inside bbox"):
        spline([1.1])


def test_automatic_knots_and_derivative_1d():
    x = np.linspace(0.0, 1.0, 30)
    y = 2.0 * x + 1.0

    spline = LSQMultivariateSpline(x=x, y=y, t=1, k=1)

    x_eval = np.array([0.2, 0.5, 0.8])
    np.testing.assert_allclose(spline(x_eval), 2.0 * x_eval + 1.0)
    np.testing.assert_allclose(spline(x_eval, nu=1), 2.0)


def test_vector_valued_output():
    x = np.linspace(0.0, 1.0, 30)
    y = np.column_stack((2.0 * x + 1.0, -3.0 * x + 4.0))

    spline = LSQMultivariateSpline(x=x, y=y, t=1, k=1)
    values = spline([0.25, 0.75])

    assert values.shape == (2, 2)
    expected = np.column_stack(
        (
            2.0 * np.array([0.25, 0.75]) + 1.0,
            -3.0 * np.array([0.25, 0.75]) + 4.0,
        )
    )
    np.testing.assert_allclose(values, expected)


def test_tensor_shaped_coefficients_for_vector_output():
    x0 = np.linspace(-1.0, 1.0, 5)
    x1 = np.linspace(-1.0, 1.0, 6)
    x0_mesh, x1_mesh = np.meshgrid(x0, x1, indexing="ij")
    y0 = x0_mesh + x1_mesh
    y1 = x0_mesh - x1_mesh
    y = np.stack((y0, y1), axis=-1)

    spline = LSQMultivariateSpline.from_grid((x0, x1), y, t=[1, 1], k=1)
    coeffs = spline.get_coeffs()
    coeffs_tensor = spline.get_coeffs_tensor()

    assert coeffs.shape == (4, 2)
    assert coeffs_tensor.shape == (2, 2, 2)
    np.testing.assert_allclose(coeffs_tensor.reshape(4, 2), coeffs)


def test_get_coeffs_returns_copy():
    x = np.linspace(0.0, 1.0, 20)
    y = 2.0 * x + 1.0
    spline = LSQMultivariateSpline(x=x, y=y, t=1, k=1)

    coeffs = spline.get_coeffs()
    coeffs[...] = 0.0

    np.testing.assert_allclose(spline([0.25]), [1.5])


def test_from_grid_scalar_values():
    x0 = np.linspace(-1.0, 1.0, 6)
    x1 = np.linspace(-1.0, 1.0, 7)
    x0_mesh, x1_mesh = np.meshgrid(x0, x1, indexing="ij")
    y = 0.5 * x0_mesh - 0.25 * x1_mesh + 2.0

    spline = LSQMultivariateSpline.from_grid((x0, x1), y, t=[1, 1], k=1)
    value = spline([0.5, -0.5])

    np.testing.assert_allclose(value, 2.375)


def test_from_grid_vector_values_output_first_axis():
    x0 = np.linspace(-1.0, 1.0, 5)
    x1 = np.linspace(-1.0, 1.0, 6)
    x0_mesh, x1_mesh = np.meshgrid(x0, x1, indexing="ij")
    y0 = 1.0 + x0_mesh
    y1 = 2.0 - x1_mesh
    y = np.stack((y0, y1), axis=0)

    spline = LSQMultivariateSpline.from_grid((x0, x1), y, t=[1, 1], k=1)
    value = spline([[0.25, -0.5]])

    np.testing.assert_allclose(value, [[1.25, 2.5]], atol=1e-12)


def test_from_grid_weights_are_applied():
    x = np.linspace(0.0, 1.0, 8)
    y = np.array([1.0, 1.1, 1.5, 2.0, 2.6, 3.1, 3.6, 4.2])
    w = np.linspace(1.0, 2.0, x.size)

    spline = LSQMultivariateSpline.from_grid((x,), y, w=w, t=1, k=1)
    residual = np.sum((w * (spline(x) - y)) ** 2)

    np.testing.assert_allclose(spline.get_residual(), residual)


def test_sparse_fit_matches_dense_fit():
    pytest.importorskip("scipy")

    x = np.linspace(0.0, 1.0, 40)
    y = 2.0 * x + 1.0

    dense = LSQMultivariateSpline(x=x, y=y, t=1, k=1)
    sparse = LSQMultivariateSpline(x=x, y=y, t=1, k=1, sparse=True)
    x_eval = np.linspace(0.1, 0.9, 5)

    np.testing.assert_allclose(sparse(x_eval), dense(x_eval), atol=1e-10)


def test_sparse_fit_matches_dense_fit_2d():
    pytest.importorskip("scipy")

    x0 = np.linspace(-1.0, 1.0, 7)
    x1 = np.linspace(-1.0, 1.0, 8)
    x0_mesh, x1_mesh = np.meshgrid(x0, x1, indexing="ij")
    x = np.column_stack((x0_mesh.ravel(), x1_mesh.ravel()))
    y = 1.0 + 2.0 * x[:, 0] - 0.5 * x[:, 1]

    dense = LSQMultivariateSpline(x=x, y=y, t=[1, 1], k=1)
    sparse = LSQMultivariateSpline(x=x, y=y, t=[1, 1], k=1, sparse=True)
    x_eval = np.array([[-0.25, 0.1], [0.5, -0.8]])

    np.testing.assert_allclose(sparse(x_eval), dense(x_eval), atol=1e-10)


def test_smoothing_reduces_coefficient_roughness():
    pytest.importorskip("scipy")

    rng = np.random.default_rng(0)
    x = np.linspace(0.0, 1.0, 80)
    y = np.sin(12.0 * x) + 0.2 * rng.normal(size=x.shape)

    rough = LSQMultivariateSpline(x=x, y=y, t=14, k=3)
    smooth = LSQMultivariateSpline(x=x, y=y, t=14, k=3, smoothing=10.0)

    roughness_plain = np.sum(np.diff(rough.get_coeffs_tensor(), n=2) ** 2)
    roughness_smooth = np.sum(np.diff(smooth.get_coeffs_tensor(), n=2) ** 2)

    assert roughness_smooth < roughness_plain
    assert smooth.get_residual() > rough.get_residual()


def test_zero_smoothing_matches_plain_least_squares():
    x = np.linspace(0.0, 1.0, 40)
    y = np.sin(4.0 * x)

    plain = LSQMultivariateSpline(x=x, y=y, t=5, k=3)
    smooth_zero = LSQMultivariateSpline(x=x, y=y, t=5, k=3, smoothing=0.0)
    x_eval = np.linspace(0.0, 1.0, 12)

    np.testing.assert_allclose(smooth_zero(x_eval), plain(x_eval))
