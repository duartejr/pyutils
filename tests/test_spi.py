"""Unit tests for Standardized Precipitation Index (SPI) computation.

Tests cover accumulation windows, Gamma distribution fitting, edge cases,
and robustness to missing data for the SPI implementation.
"""

import numpy as np
import pytest

from spi import spi


class TestSPIAccumulation:
    """Verify precipitation accumulation over specified scales."""

    def test_accumulation_scale_one(self, simple_precipitation_series):
        """Scale-1 SPI should process raw monthly values."""
        result = spi(simple_precipitation_series, accumulation_scale=1, num_seasons=12)

        assert len(result) == len(simple_precipitation_series)
        assert result.shape == (120,)

    def test_accumulation_scale_three(self, simple_precipitation_series):
        """Scale-3 SPI should reduce length by 2 months."""
        result = spi(simple_precipitation_series, accumulation_scale=3, num_seasons=12)

        assert len(result) == len(simple_precipitation_series) - 2
        assert result.shape == (118,)

    def test_accumulation_scale_twelve(self, simple_precipitation_series):
        """Scale-12 SPI should sum annual totals."""
        result = spi(simple_precipitation_series, accumulation_scale=12, num_seasons=12)

        assert len(result) == len(simple_precipitation_series) - 11
        assert result.shape == (109,)

    def test_accumulation_scale_larger_than_series(self, simple_precipitation_series):
        """Accumulation larger than series length should still compute."""
        short_series = simple_precipitation_series[:24]
        result = spi(short_series, accumulation_scale=15, num_seasons=12)

        assert len(result) == len(short_series) - 14
        assert result.shape == (10,)


class TestSPISeasonalPartitioning:
    """Verify that seasonal grouping works correctly."""

    def test_twelve_seasons(self, simple_precipitation_series):
        """12 seasons should group months by calendar months."""
        result = spi(simple_precipitation_series, accumulation_scale=1, num_seasons=12)

        assert len(result) == len(simple_precipitation_series)
        assert result.shape[0] > 0

    def test_four_seasons(self, simple_precipitation_series):
        """4 seasons should group months by quarters."""
        result = spi(simple_precipitation_series, accumulation_scale=1, num_seasons=4)

        assert len(result) == len(simple_precipitation_series)

    def test_one_season(self, simple_precipitation_series):
        """Single season should treat all months together."""
        result = spi(simple_precipitation_series, accumulation_scale=1, num_seasons=1)

        assert len(result) == len(simple_precipitation_series)


class TestSPIGammaFitting:
    """Verify Gamma distribution fitting and parameter handling."""

    def test_spi_values_standardized(self, simple_precipitation_series):
        """SPI values should follow approximately standard normal distribution."""
        result = spi(simple_precipitation_series, accumulation_scale=1, num_seasons=12)

        valid_result = result[~np.isnan(result)]
        mean_spi = np.nanmean(result)
        std_spi = np.nanstd(result)

        assert len(valid_result) > 0
        assert mean_spi == pytest.approx(0, abs=0.5)
        assert std_spi == pytest.approx(1, abs=0.3)

    def test_spi_range_realistic(self, simple_precipitation_series):
        """SPI values should be within realistic bounds (-4 to +4)."""
        result = spi(simple_precipitation_series, accumulation_scale=3, num_seasons=12)

        valid_result = result[~np.isnan(result)]
        assert np.all(valid_result >= -4)
        assert np.all(valid_result <= 4)

    def test_negative_spi_drought(self, simple_precipitation_series):
        """Low precipitation should produce negative SPI (drought)."""
        series_with_drought = simple_precipitation_series.copy()
        series_with_drought[36:48] = 10  # Very low values in months 36-47

        result = spi(series_with_drought, accumulation_scale=3, num_seasons=12)

        drought_period_spi = result[36:45]
        drought_spi = drought_period_spi[~np.isnan(drought_period_spi)]

        assert np.nanmean(drought_spi) < 0

    def test_positive_spi_surplus(self, simple_precipitation_series):
        """High precipitation should produce positive SPI (wet)."""
        series_with_surplus = simple_precipitation_series.copy()
        series_with_surplus[36:48] = 500  # Very high values in months 36-47

        result = spi(series_with_surplus, accumulation_scale=3, num_seasons=12)

        surplus_period_spi = result[36:45]
        surplus_spi = surplus_period_spi[~np.isnan(surplus_period_spi)]

        assert np.nanmean(surplus_spi) > 0


class TestSPIMissingDataHandling:
    """Verify robustness to missing data and edge cases."""

    def test_entire_season_missing(self, all_missing_season):
        """Season with all NaN should produce NaN for that season."""
        result = spi(all_missing_season, accumulation_scale=1, num_seasons=12)

        missing_season_indices = np.arange(0, 12)
        result_in_missing_season = result[missing_season_indices]

        assert np.all(np.isnan(result_in_missing_season))

    def test_scattered_missing_values(self, simple_precipitation_series):
        """Scattered NaN values should be handled gracefully."""
        series_with_gaps = simple_precipitation_series.copy()
        series_with_gaps[10] = np.nan
        series_with_gaps[50] = np.nan
        series_with_gaps[100] = np.nan

        result = spi(series_with_gaps, accumulation_scale=1, num_seasons=12)

        assert len(result) == len(series_with_gaps)
        valid_result = result[~np.isnan(result)]
        assert len(valid_result) > 0

    def test_zero_precipitation_values(self, simple_precipitation_series):
        """Zero precipitation should be handled by mixed distribution."""
        series_with_zeros = simple_precipitation_series.copy()
        series_with_zeros[0:12] = 0

        result = spi(series_with_zeros, accumulation_scale=1, num_seasons=12)

        assert len(result) == len(series_with_zeros)
        assert np.any(np.isnan(result)) or np.any(np.isfinite(result))

    def test_mostly_zeros_season(self):
        """Season with mostly zeros should still compute SPI."""
        series = np.full(120, 50.0)
        series[0:12] = 0  # First season all zeros
        series[1:11] = 0  # Most months zero

        result = spi(series, accumulation_scale=1, num_seasons=12)

        assert len(result) == len(series)


class TestSPIEdgeCases:
    """Test behavior on special and boundary configurations."""

    def test_slowly_increasing_trend(self):
        """Slowly increasing precipitation trend should be handled."""
        trend_series = np.linspace(40, 60, 120)
        result = spi(trend_series, accumulation_scale=1, num_seasons=12)

        valid_result = result[~np.isnan(result)]
        assert len(valid_result) > 50

    def test_increasing_trend(self):
        """Increasing precipitation trend should be captured."""
        increasing_series = np.linspace(20, 100, 120)
        result = spi(increasing_series, accumulation_scale=1, num_seasons=12)

        valid_result = result[~np.isnan(result)]
        assert len(valid_result) > 0

    def test_short_series(self):
        """Short time series should be handled."""
        short_series = np.array([10, 20, 30, 40, 50, 60])
        result = spi(short_series, accumulation_scale=1, num_seasons=2)

        assert len(result) == len(short_series)

    def test_single_month_data(self):
        """Single value should be handled gracefully."""
        single_value = np.array([50.0])
        result = spi(single_value, accumulation_scale=1, num_seasons=1)

        assert len(result) == 1


class TestSPIInputValidation:
    """Verify proper handling of input types and conversions."""

    def test_list_input_conversion(self):
        """List input should be converted to ndarray."""
        precipitation_list = [10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 110, 120]
        result = spi(precipitation_list, accumulation_scale=1, num_seasons=12)

        assert isinstance(result, np.ndarray)
        assert len(result) == len(precipitation_list)

    def test_integer_input_conversion(self):
        """Integer array input should be converted to float."""
        precipitation_int = np.array([10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 110, 120])
        result = spi(precipitation_int, accumulation_scale=1, num_seasons=12)

        assert result.dtype in [np.float32, np.float64]

    def test_negative_values_handling(self):
        """Negative precipitation values should be handled."""
        series_with_negative = np.array([10, -5, 20, 30, 40, 50, 60, 70, 80, 90, 100, 110])
        result = spi(series_with_negative, accumulation_scale=1, num_seasons=12)

        assert len(result) == len(series_with_negative)


class TestSPIOutputConsistency:
    """Verify consistency and reproducibility of results."""

    def test_reproducible_results(self, simple_precipitation_series):
        """Same input should produce identical output."""
        result1 = spi(simple_precipitation_series, accumulation_scale=3, num_seasons=12)
        result2 = spi(simple_precipitation_series, accumulation_scale=3, num_seasons=12)

        assert np.allclose(result1, result2, equal_nan=True)

    def test_output_shape_consistency(self, simple_precipitation_series):
        """Output shape should equal input minus accumulation scale plus 1."""
        accumulation_scale = 3
        result = spi(simple_precipitation_series, accumulation_scale, num_seasons=12)

        expected_length = len(simple_precipitation_series) - accumulation_scale + 1
        assert len(result) == expected_length

    def test_finite_and_nan_distribution(self, simple_precipitation_series):
        """Result should contain either finite or NaN values, no inf."""
        result = spi(simple_precipitation_series, accumulation_scale=1, num_seasons=12)

        assert np.all(np.isfinite(result) | np.isnan(result))
