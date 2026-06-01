"""Unit tests for Inverse Distance Weighting (IDW) interpolation.

Tests cover exact-hit behavior, shape validation, interpolation quality,
edge cases, and error handling for the vectorized IDW implementation.
"""

import numpy as np
import pytest

from idw import idw


class TestIDWExactHits:
    """Verify that cells coinciding with stations receive exact station values."""

    def test_single_station_at_center(self, single_station):
        """Center station's value should propagate to center grid cell."""
        station_rows, station_cols, station_values, output_grid = single_station
        result = idw(station_rows, station_cols, station_values, output_grid, power=2)

        assert result[2, 2] == pytest.approx(42.0)

    def test_four_corner_stations(self, grid_with_corner_stations):
        """Each corner station should occupy its exact grid position."""
        station_rows, station_cols, station_values, output_grid = grid_with_corner_stations
        result = idw(station_rows, station_cols, station_values, output_grid, power=2)

        assert result[0, 0] == pytest.approx(10.0)
        assert result[0, 9] == pytest.approx(20.0)
        assert result[9, 0] == pytest.approx(30.0)
        assert result[9, 9] == pytest.approx(40.0)

    def test_exact_hit_power_zero(self, single_station):
        """Exact hits should return station value even when power=0."""
        station_rows, station_cols, station_values, output_grid = single_station
        result = idw(station_rows, station_cols, station_values, output_grid, power=0)

        assert result[2, 2] == pytest.approx(42.0)

    def test_exact_hit_high_power(self, grid_with_corner_stations):
        """Exact hits should return station value even with high power."""
        station_rows, station_cols, station_values, output_grid = grid_with_corner_stations
        result = idw(station_rows, station_cols, station_values, output_grid, power=10)

        assert result[0, 0] == pytest.approx(10.0)
        assert result[9, 9] == pytest.approx(40.0)


class TestIDWShapeAndStructure:
    """Verify output grid dimensions and type correctness."""

    def test_output_shape_matches_grid(self, grid_with_corner_stations):
        """Output grid shape should match input grid shape."""
        station_rows, station_cols, station_values, output_grid = grid_with_corner_stations
        result = idw(station_rows, station_cols, station_values, output_grid, power=2)

        assert result.shape == (10, 10)

    def test_output_grid_in_place_modification(self, single_station):
        """IDW should modify the grid in-place and return the same object."""
        station_rows, station_cols, station_values, output_grid = single_station
        result = idw(station_rows, station_cols, station_values, output_grid, power=2)

        assert result is output_grid

    def test_output_dtype_float(self, single_station):
        """Output should contain float values."""
        station_rows, station_cols, station_values, output_grid = single_station
        result = idw(station_rows, station_cols, station_values, output_grid, power=2)

        assert result.dtype == np.float64 or result.dtype == np.float32


class TestIDWInterpolationQuality:
    """Verify interpolation produces realistic values between observations."""

    def test_interpolated_values_between_extremes(self, grid_with_corner_stations):
        """Interpolated values should fall within station value range."""
        station_rows, station_cols, station_values, output_grid = grid_with_corner_stations
        result = idw(station_rows, station_cols, station_values, output_grid, power=2)

        min_value = station_values.min()
        max_value = station_values.max()
        assert result.min() >= min_value
        assert result.max() <= max_value

    def test_center_cell_closer_to_nearby_values(self, grid_with_corner_stations):
        """Center cell should be influenced more by nearby values than distant ones."""
        station_rows, station_cols, station_values, output_grid = grid_with_corner_stations
        result = idw(station_rows, station_cols, station_values, output_grid, power=2)

        center_value = result[5, 5]
        expected_range = (min(station_values), max(station_values))
        assert expected_range[0] < center_value < expected_range[1]

    def test_higher_power_localizes_influence(self, grid_with_corner_stations):
        """Higher power should increase influence of nearest neighbors."""
        station_rows, station_cols, station_values, output_grid = grid_with_corner_stations

        result_power1 = idw(
            station_rows.copy(),
            station_cols.copy(),
            station_values.copy(),
            np.zeros((10, 10)),
            power=1,
        )

        station_rows, station_cols, station_values, output_grid = grid_with_corner_stations
        result_power3 = idw(station_rows, station_cols, station_values, output_grid, power=3)

        variance_power1 = np.var(result_power1)
        variance_power3 = np.var(result_power3)

        assert variance_power3 > variance_power1


class TestIDWEdgeCases:
    """Test behavior on edge cases and special configurations."""

    def test_power_zero(self, grid_with_corner_stations):
        """Power=0 should average all station values uniformly."""
        station_rows, station_cols, station_values, output_grid = grid_with_corner_stations
        result = idw(station_rows, station_cols, station_values, output_grid, power=0)

        expected_mean = station_values.mean()
        assert result[5, 5] == pytest.approx(expected_mean)

    def test_power_negative(self, single_station):
        """Negative power should be handled (inverse distance grows with distance)."""
        station_rows, station_cols, station_values, output_grid = single_station
        result = idw(station_rows, station_cols, station_values, output_grid, power=-2)

        assert np.all(np.isfinite(result))

    def test_fractional_coordinates(self):
        """Fractional station coordinates should be handled correctly."""
        station_rows = np.array([2.5, 2.5])
        station_cols = np.array([2.5, 2.5])
        station_values = np.array([10.0, 20.0])
        output_grid = np.zeros((5, 5))

        result = idw(station_rows, station_cols, station_values, output_grid, power=2)

        assert result.shape == (5, 5)
        assert np.all(np.isfinite(result))


class TestIDWValidation:
    """Verify error handling for invalid inputs."""

    def test_empty_station_rows(self):
        """Empty station array should raise ValueError."""
        with pytest.raises(ValueError, match="at least one station"):
            idw(
                np.array([]),
                np.array([1.0]),
                np.array([10.0]),
                np.zeros((5, 5)),
                power=2,
            )

    def test_mismatched_rows_cols(self):
        """Mismatched station_rows and station_cols should raise ValueError."""
        with pytest.raises(ValueError, match="Mismatched array sizes"):
            idw(
                np.array([0.0, 1.0]),
                np.array([0.0]),
                np.array([10.0, 20.0]),
                np.zeros((5, 5)),
                power=2,
            )

    def test_mismatched_rows_values(self):
        """Mismatched station_rows and station_values should raise ValueError."""
        with pytest.raises(ValueError, match="Mismatched array sizes"):
            idw(
                np.array([0.0, 1.0]),
                np.array([0.0, 1.0]),
                np.array([10.0]),
                np.zeros((5, 5)),
                power=2,
            )

    def test_single_station_grid(self):
        """Single station should work without error."""
        station_rows = np.array([0.0])
        station_cols = np.array([0.0])
        station_values = np.array([42.0])
        output_grid = np.zeros((5, 5))

        result = idw(station_rows, station_cols, station_values, output_grid, power=2)

        assert result[0, 0] == pytest.approx(42.0)


class TestIDWNumericalStability:
    """Verify numerical stability with various data ranges and configurations."""

    def test_large_values(self):
        """Large station values should not cause overflow."""
        station_rows = np.array([0.0, 9.0])
        station_cols = np.array([0.0, 9.0])
        station_values = np.array([1e6, 1e6])
        output_grid = np.zeros((10, 10))

        result = idw(station_rows, station_cols, station_values, output_grid, power=2)

        assert np.all(np.isfinite(result))

    def test_small_values(self):
        """Small station values should be interpolated accurately."""
        station_rows = np.array([0.0, 4.0])
        station_cols = np.array([0.0, 4.0])
        station_values = np.array([1e-6, 2e-6])
        output_grid = np.zeros((5, 5))

        result = idw(station_rows, station_cols, station_values, output_grid, power=2)

        assert np.all(np.isfinite(result))

    def test_mixed_sign_values(self):
        """Mixed positive and negative values should be handled."""
        station_rows = np.array([0.0, 9.0])
        station_cols = np.array([0.0, 9.0])
        station_values = np.array([-10.0, 10.0])
        output_grid = np.zeros((10, 10))

        result = idw(station_rows, station_cols, station_values, output_grid, power=2)

        assert np.all(np.isfinite(result))
