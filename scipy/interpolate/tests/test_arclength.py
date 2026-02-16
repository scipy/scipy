"""
Tests for arc-length parameterization of splines.
"""
import numpy as np
from numpy.testing import assert_allclose, assert_array_less
import pytest
from scipy.interpolate import splprep, splev, splev_arclength, spline_arclength


class TestArcLength:
    """Tests for arc-length spline functionality."""

    def test_circle_arclength(self):
        """Test that a circle has the correct arc length."""
        # Create a unit circle
        theta = np.linspace(0, 2*np.pi, 50, endpoint=False)
        x = np.cos(theta)
        y = np.sin(theta)
        
        # Fit a periodic spline
        tck, u = splprep([x, y], s=0, per=True)
        
        # Compute arc length
        length = spline_arclength(tck)
        
        # Should be close to 2*pi
        assert_allclose(length, 2*np.pi, rtol=1e-3)

    def test_equal_spacing_circle(self):
        """Test that points sampled with equal arc length are equally spaced."""
        # Create a circle
        theta = np.linspace(0, 2*np.pi, 50, endpoint=False)
        x = np.cos(theta)
        y = np.sin(theta)
        
        # Fit a periodic spline
        tck, u = splprep([x, y], s=0, per=True)
        
        # Sample with equal arc-length spacing
        n_samples = 100
        points, u_values = splev_arclength(tck, n_samples)
        
        # Compute distances between consecutive points
        diffs = np.diff(points, axis=0)
        distances = np.sqrt(np.sum(diffs**2, axis=1))
        
        # All distances should be approximately equal
        mean_dist = np.mean(distances)
        # Allow 5% tolerance
        assert_allclose(distances, mean_dist, rtol=0.05)

    def test_line_segment(self):
        """Test arc-length sampling on a straight line."""
        # Create a line from (0,0) to (3,4), length = 5
        t_param = np.linspace(0, 1, 10)
        x = 3 * t_param
        y = 4 * t_param
        
        # Fit a spline
        tck, u = splprep([x, y], s=0)
        
        # Compute arc length
        length = spline_arclength(tck)
        assert_allclose(length, 5.0, rtol=1e-3)
        
        # Sample with equal arc-length spacing
        points, u_values = splev_arclength(tck, 11)
        
        # Compute distances
        diffs = np.diff(points, axis=0)
        distances = np.sqrt(np.sum(diffs**2, axis=1))
        
        # Should all be 0.5
        assert_allclose(distances, 0.5, rtol=1e-2)

    def test_1d_spline(self):
        """Test that 1D splines work correctly."""
        from scipy.interpolate import splrep
        
        # Simple 1D function
        x = np.linspace(0, 2*np.pi, 20)
        y = np.sin(x)
        
        # Fit spline
        tck = splrep(x, y, s=0)
        
        # Compute arc length
        length = spline_arclength(tck)
        
        # Arc length of sin from 0 to 2pi is approximately 7.64
        assert length > 7.5
        assert length < 8.0

    def test_return_shapes(self):
        """Test that return values have correct shapes."""
        # 2D parametric curve
        theta = np.linspace(0, np.pi, 20)
        x = np.cos(theta)
        y = np.sin(theta)
        
        tck, u = splprep([x, y], s=0)
        
        n_samples = 50
        points, u_values = splev_arclength(tck, n_samples)
        
        # Check shapes
        assert points.shape == (n_samples, 2)
        assert u_values.shape == (n_samples,)
        
        # Check parameter ordering
        assert np.all(np.diff(u_values) >= 0)  # Should be monotonically increasing

    def test_partial_range(self):
        """Test arc-length computation on a partial parameter range."""
        # Create a circle
        theta = np.linspace(0, 2*np.pi, 50, endpoint=False)
        x = np.cos(theta)
        y = np.sin(theta)
        
        tck, u = splprep([x, y], s=0, per=True)
        
        # Compute arc length for half the curve
        u_min, u_max = 0.0, 0.5
        half_length = spline_arclength(tck, u_range=(u_min, u_max))
        
        # Should be approximately pi (half of 2*pi)
        full_length = spline_arclength(tck)
        assert_allclose(half_length, full_length / 2, rtol=0.05)

    def test_3d_curve(self):
        """Test that 3D curves work correctly."""
        # Create a helix
        t = np.linspace(0, 4*np.pi, 100)
        x = np.cos(t)
        y = np.sin(t)
        z = t
        
        # Fit spline
        tck, u = splprep([x, y, z], s=0)
        
        # Sample with equal arc-length spacing
        n_samples = 200
        points, u_values = splev_arclength(tck, n_samples)
        
        # Check shape
        assert points.shape == (n_samples, 3)
        
        # Compute distances
        diffs = np.diff(points, axis=0)
        distances = np.sqrt(np.sum(diffs**2, axis=1))
        
        # All distances should be approximately equal
        mean_dist = np.mean(distances)
        # Allow 10% tolerance for helix (more complex curve)
        assert_allclose(distances, mean_dist, rtol=0.1)

    def test_degenerate_curve(self):
        """Test handling of degenerate (zero-length) curves."""
        # All points the same
        x = np.ones(10)
        y = np.ones(10)
        
        tck, u = splprep([x, y], s=0)
        
        # Arc length should be zero
        length = spline_arclength(tck)
        assert_allclose(length, 0.0, atol=1e-10)
        
        # Sampling should not crash
        points, u_values = splev_arclength(tck, 10)
        assert points.shape == (10, 2)

    def test_comparison_with_splev(self):
        """Test that endpoints match regular splev evaluation."""
        theta = np.linspace(0, np.pi, 20)
        x = np.cos(theta)
        y = np.sin(theta)
        
        tck, u = splprep([x, y], s=0)
        
        # Sample with arc-length
        points, u_values = splev_arclength(tck, 50)
        
        # First and last points should match curve endpoints
        t_knots = tck[0]
        k = tck[2]
        u_min = t_knots[k]
        u_max = t_knots[-k-1]
        
        start_point = np.array(splev(u_min, tck))
        end_point = np.array(splev(u_max, tck))
        
        assert_allclose(points[0], start_point, rtol=1e-5)
        assert_allclose(points[-1], end_point, rtol=1e-5)


class TestArcLengthEdgeCases:
    """Test edge cases and error handling."""

    def test_small_n_samples(self):
        """Test with very small number of samples."""
        theta = np.linspace(0, np.pi, 20)
        x = np.cos(theta)
        y = np.sin(theta)
        
        tck, u = splprep([x, y], s=0)
        
        # Should work with n_samples=2
        points, u_values = splev_arclength(tck, 2)
        assert points.shape == (2, 2)
        
        # Should work with n_samples=1
        points, u_values = splev_arclength(tck, 1)
        assert points.shape == (1, 2)

    def test_high_curvature(self):
        """Test with high curvature regions."""
        # Create a sharp curve
        t = np.linspace(0, 1, 100)
        x = t
        y = np.where(t < 0.5, 0, 10 * (t - 0.5)**2)
        
        tck, u = splprep([x, y], s=0)
        
        # Should handle high curvature
        points, u_values = splev_arclength(tck, 100)
        
        # Check that spacing is more uniform than uniform parameter spacing
        diffs = np.diff(points, axis=0)
        arc_distances = np.sqrt(np.sum(diffs**2, axis=1))
        
        # Compare with uniform parameter spacing
        u_uniform = np.linspace(u[0], u[-1], 100)
        points_uniform = np.column_stack(splev(u_uniform, tck))
        diffs_uniform = np.diff(points_uniform, axis=0)
        uniform_distances = np.sqrt(np.sum(diffs_uniform**2, axis=1))
        
        # Arc-length spacing should be more uniform (lower std dev)
        assert np.std(arc_distances) < np.std(uniform_distances)
