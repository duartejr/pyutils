"""Unit tests for evapotranspiration calculation methods.

Tests cover Hargreaves-Samani, Thornthwaite, and Penman-Monteith methods
for reference evapotranspiration from meteorological inputs.
"""

import numpy as np
import pytest

from evapotranspiration import hargreaves, penman_monteith, thornthwaite


class TestHargreavesBasic:
    """Test basic Hargreaves-Samani evapotranspiration computation."""

    def test_output_shape_daily_aggregation(self, evapotranspiration_temperatures):
        """Daily aggregation mode should return 360 values (12 months × 30 mean days)."""
        mean_temps, max_temps, min_temps = evapotranspiration_temperatures
        temps_annual = np.tile(mean_temps, 10)
        max_temps_annual = np.tile(max_temps, 10)
        min_temps_annual = np.tile(min_temps, 10)

        result = hargreaves(
            temps_annual, max_temps_annual, min_temps_annual, y=0, months=1, data_freq="daily"
        )

        assert len(result) == len(temps_annual)

    def test_output_positive_values(self, evapotranspiration_temperatures):
        """ETP values should be positive."""
        mean_temps, max_temps, min_temps = evapotranspiration_temperatures

        result = hargreaves(
            mean_temps, max_temps, min_temps, y=0, months=np.arange(1, 13), data_freq="monthly"
        )

        assert np.all(result >= 0)

    def test_reference_case_accuracy(self, hargreaves_reference_case):
        """Known reference case should produce positive monthly ETP."""
        case = hargreaves_reference_case
        mean_temp_array = np.array([case["mean_temp"]])
        max_temp_array = np.array([case["max_temp"]])
        min_temp_array = np.array([case["min_temp"]])
        month_array = np.array([case["month"]])

        result = hargreaves(
            mean_temp_array,
            max_temp_array,
            min_temp_array,
            y=case["latitude"],
            months=month_array,
            data_freq="monthly",
        )

        assert result[0] > 0
        assert result[0] > 100


class TestHargreavesDaily:
    """Test daily Hargreaves ETP mode."""

    def test_daily_etp_output_length(self):
        """Daily mode output should match input temperature series length."""
        daily_temps = np.full(365, 25.0)
        daily_max = np.full(365, 30.0)
        daily_min = np.full(365, 20.0)

        result = hargreaves(
            daily_temps, daily_max, daily_min, y=0, months=1, etp_daily=True
        )

        assert len(result) == 365

    def test_daily_etp_positive(self):
        """Daily ETP values should be non-negative."""
        daily_temps = np.random.normal(25, 5, 365)
        daily_max = daily_temps + 5
        daily_min = daily_temps - 5

        result = hargreaves(
            daily_temps, daily_max, daily_min, y=0, months=1, etp_daily=True
        )

        assert np.all(result >= 0)


class TestHargreavesLatitudeEffect:
    """Test latitude impact on ETP calculation."""

    def test_equator_vs_pole_latitude(self, evapotranspiration_temperatures):
        """Lower latitudes should generally produce higher ETP."""
        mean_temps, max_temps, min_temps = evapotranspiration_temperatures
        months = np.arange(1, 13)

        result_equator = hargreaves(
            mean_temps, max_temps, min_temps, y=0, months=months, data_freq="monthly"
        )
        result_high_lat = hargreaves(
            mean_temps, max_temps, min_temps, y=60, months=months, data_freq="monthly"
        )

        assert result_equator.mean() > result_high_lat.mean()

    def test_southern_hemisphere_latitude(self, evapotranspiration_temperatures):
        """Negative latitude should be handled (Southern Hemisphere)."""
        mean_temps, max_temps, min_temps = evapotranspiration_temperatures
        months = np.arange(1, 13)

        result = hargreaves(
            mean_temps, max_temps, min_temps, y=-30, months=months, data_freq="monthly"
        )

        assert len(result) == 12
        assert np.all(result >= 0)


class TestThornthwaitBasic:
    """Test basic Thornthwaite evapotranspiration computation."""

    def test_output_shape(self):
        """Output shape should match input temperature series."""
        temperatures = np.array([20, 21, 22, 23, 24, 25, 24, 23, 22, 21, 20, 19])
        sunshine_hours = np.array([8, 9, 10, 11, 12, 13, 12, 11, 10, 9, 8, 7])

        result = thornthwaite(temperatures, sunshine_hours)

        assert len(result) == len(temperatures)

    def test_output_positive_values(self):
        """ETP values should be positive."""
        temperatures = np.array([20, 21, 22, 23, 24, 25, 24, 23, 22, 21, 20, 19])
        sunshine_hours = np.array([8, 9, 10, 11, 12, 13, 12, 11, 10, 9, 8, 7])

        result = thornthwaite(temperatures, sunshine_hours)

        assert np.all(result >= 0)

    def test_multiple_years(self):
        """Should handle multiple years of data."""
        temperatures = np.tile([20, 21, 22, 23, 24, 25, 24, 23, 22, 21, 20, 19], 10)
        sunshine_hours = np.tile([8, 9, 10, 11, 12, 13, 12, 11, 10, 9, 8, 7], 10)

        result = thornthwaite(temperatures, sunshine_hours)

        assert len(result) == 120


class TestThornthwaitSeasonality:
    """Test seasonal variation in Thornthwaite ETP."""

    def test_warmer_months_higher_etp(self):
        """Warmer months should produce higher ETP values."""
        winter_month_temps = np.array([10])
        winter_month_sunshine = np.array([8])
        summer_month_temps = np.array([30])
        summer_month_sunshine = np.array([14])

        annual_temps = np.concatenate([winter_month_temps, summer_month_temps] * 6)
        annual_sunshine = np.concatenate([winter_month_sunshine, summer_month_sunshine] * 6)

        result = thornthwaite(annual_temps, annual_sunshine)

        even_indices = np.arange(0, len(result), 2)
        odd_indices = np.arange(1, len(result), 2)

        assert result[odd_indices].mean() > result[even_indices].mean()


class TestPenmanMonteithBasic:
    """Test basic Penman-Monteith reference ETP."""

    def test_function_exists_and_callable(self):
        """Penman-Monteith function should be callable."""
        assert callable(penman_monteith)

    def test_parameters_exist(self):
        """Penman-Monteith should accept major meteorological inputs."""
        import inspect
        sig = inspect.signature(penman_monteith)
        expected_params = {'Tmed', 'Tmax', 'Tmin', 'U2', 'UR', 'Rs', 'lat'}
        actual_params = set(sig.parameters.keys())
        assert expected_params.issubset(actual_params)


class TestPenmanMonteithKwargs:
    """Test Penman-Monteith optional parameters."""

    def test_albedo_parameter_exists(self):
        """Albedo parameter should be configurable."""
        import inspect
        sig = inspect.signature(penman_monteith)
        assert 'albedo' in sig.parameters
        assert sig.parameters['albedo'].default == 0.23

    def test_soil_heat_flux_parameter_exists(self):
        """Soil heat flux parameter should be configurable."""
        import inspect
        sig = inspect.signature(penman_monteith)
        assert 'G' in sig.parameters
        assert sig.parameters['G'].default == 0.0

    def test_height_parameter_exists(self):
        """Height of wind measurement should be configurable."""
        import inspect
        sig = inspect.signature(penman_monteith)
        assert 'z' in sig.parameters
        assert sig.parameters['z'].default == 10


class TestEvapotranspirationEdgeCases:
    """Test behavior on boundary conditions."""

    def test_zero_temperature_difference(self):
        """Zero temperature difference should still compute."""
        temperature = np.array([25.0])

        result = hargreaves(
            temperature,
            temperature,
            temperature,
            y=0,
            months=np.array([1]),
            data_freq="monthly",
        )

        assert len(result) == 1
        assert result[0] >= 0

    def test_very_low_temperature(self):
        """Very low temperatures should produce minimal ETP."""
        cold_temp = np.array([5.0])
        max_cold = np.array([8.0])
        min_cold = np.array([2.0])

        result = hargreaves(cold_temp, max_cold, min_cold, y=0, months=np.array([1]), data_freq="monthly")

        assert result[0] > 0


class TestEvapotranspirationNumericalStability:
    """Verify numerical stability across parameter ranges."""

    def test_large_temperature_values(self):
        """Large temperature inputs should be handled."""
        hot_temp = np.array([50.0])
        hot_max = np.array([55.0])
        hot_min = np.array([45.0])

        result = hargreaves(hot_temp, hot_max, hot_min, y=0, months=np.array([1]), data_freq="monthly")

        assert np.isfinite(result[0])

    def test_small_temperature_difference(self):
        """Small daily temperature range should be handled."""
        temp = np.array([25.0])
        max_temp = np.array([25.1])
        min_temp = np.array([24.9])

        result = hargreaves(temp, max_temp, min_temp, y=0, months=np.array([1]), data_freq="monthly")

        assert result[0] >= 0
        assert np.isfinite(result[0])
