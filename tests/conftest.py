"""Shared test fixtures for pyutils test suite.

Provides synthetic data, reference implementations, and known values for
testing evapotranspiration, interpolation, and precipitation analysis functions.
"""

import numpy as np
import pytest


@pytest.fixture
def simple_precipitation_series() -> np.ndarray:
    """Generate a 10-year synthetic monthly precipitation time series.

    Returns:
        Array of shape (120,) representing monthly totals in millimeters.
        Values follow a realistic exponential distribution (mean ~80 mm).
    """
    rng = np.random.default_rng(seed=42)
    return rng.exponential(scale=80, size=120)


@pytest.fixture
def all_missing_season() -> np.ndarray:
    """Generate precipitation data with an entire season missing (all NaN).

    Returns:
        Array of shape (120,) with first 12 months set to NaN.
        Useful for testing robustness to missing data.
    """
    data = np.full(120, 50.0)
    data[:12] = np.nan
    return data


@pytest.fixture
def evapotranspiration_temperatures() -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Generate synthetic temperature data for evapotranspiration calculations.

    Returns:
        Tuple of (mean_temps, max_temps, min_temps), each shape (12,).
        Represents realistic monthly variation (tropical location, ~25°C mean).
    """
    mean_temperatures = np.array(
        [25.0, 25.5, 26.0, 25.5, 25.0, 24.5, 24.0, 24.0, 24.5, 25.0, 25.5, 25.5]
    )
    max_temperatures = mean_temperatures + 5.0
    min_temperatures = mean_temperatures - 5.0
    return mean_temperatures, max_temperatures, min_temperatures


@pytest.fixture
def grid_with_corner_stations() -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Generate a 10×10 grid with stations at all four corners.

    Returns:
        Tuple of (station_rows, station_cols, station_values, output_grid).
        Stations are at (0, 0), (0, 9), (9, 0), (9, 9) with values 10-40.
        Grid is pre-allocated and ready for interpolation.
    """
    station_rows = np.array([0.0, 0.0, 9.0, 9.0])
    station_cols = np.array([0.0, 9.0, 0.0, 9.0])
    station_values = np.array([10.0, 20.0, 30.0, 40.0])
    output_grid = np.zeros((10, 10))
    return station_rows, station_cols, station_values, output_grid


@pytest.fixture
def single_station() -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Generate a 5×5 grid with a single station at center.

    Returns:
        Tuple of (station_rows, station_cols, station_values, output_grid).
        Station at (2, 2) with value 42.0.
    """
    station_rows = np.array([2.0])
    station_cols = np.array([2.0])
    station_values = np.array([42.0])
    output_grid = np.zeros((5, 5))
    return station_rows, station_cols, station_values, output_grid


@pytest.fixture
def hargreaves_reference_case() -> dict:
    """Return known evapotranspiration values from Hargreaves method.

    For mean_temp=25°C, max=30°C, min=20°C, latitude=-5°, month=1,
    reference sources give consistent ETP values.

    Returns:
        Dictionary with keys:
        - 'mean_temp': Input mean temperature (float).
        - 'max_temp': Input maximum temperature (float).
        - 'min_temp': Input minimum temperature (float).
        - 'latitude': Latitude in degrees (float, negative for S).
        - 'month': Month (1-12) for which reference is computed (int).
        - 'expected_etp': Expected ETP in mm/day (float, approximate).
        - 'tolerance': Acceptable error margin for test (float).
    """
    return {
        "mean_temp": 25.0,
        "max_temp": 30.0,
        "min_temp": 20.0,
        "latitude": -5.0,
        "month": 1,
        "expected_etp": 6.5,
        "tolerance": 0.5,
    }
