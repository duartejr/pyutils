"""Standardized Precipitation Index (SPI) computation.

The SPI quantifies precipitation anomalies relative to a long-term baseline.
Each month's accumulated precipitation is fitted to a Gamma distribution per
season, then transformed to a standard normal variate. Negative SPI values
indicate drought; positive values indicate wet conditions.

References:
    McKee, T. B., Doesken, N. J., & Kleist, J. (1993). The relationship of
    drought frequency and duration to time scales. Proceedings of the 8th
    Conference on Applied Climatology, 179-184.

    Adapted from Taesam Lee (Dec. 2009), INRS-ETE, Quebec, Canada.
"""

__author__ = ["Paulo Jarbas Camurca"]
__credits__ = ["Fco Vasconcelos", "Marcelo Rodrigues", "Daniel Pinto"]
__license__ = "GPL"
__version__ = "1.0"
__email__ = "pjarbas312@gmail.com"

import numpy as np
from scipy.stats import gamma, norm


def _accumulate_precipitation(
    monthly_values: np.ndarray,
    accumulation_scale: int,
) -> np.ndarray:
    """Slide a window of length *accumulation_scale* over monthly precipitation.

    Args:
        monthly_values: Monthly precipitation time series. Shape ``(n_months,)``.
        accumulation_scale: Number of months to accumulate (e.g. 1, 3, 6, 12).

    Returns:
        Accumulated precipitation for each valid time step.
        Shape ``(n_months - accumulation_scale + 1,)``.
    """
    window_slices = [
        monthly_values[offset : len(monthly_values) - accumulation_scale + offset + 1]
        for offset in range(accumulation_scale)
    ]
    stacked_windows = np.array(window_slices).T  # (n_valid_steps, accumulation_scale)

    if stacked_windows.ndim == 2:
        return stacked_windows.sum(axis=1)
    return stacked_windows


def _remove_zeros_and_nans(values: np.ndarray) -> np.ndarray:
    """Return *values* with zeros and NaNs removed.

    Args:
        values: Input array, may contain zeros and NaNs.

    Returns:
        Filtered array containing only finite, positive values.
    """
    nonzero_values = values[np.nonzero(values)]
    return nonzero_values[~np.isnan(nonzero_values)]


def _compute_zero_probability(season_values: np.ndarray) -> float:
    """Estimate the probability of a zero-precipitation month in a season.

    Args:
        season_values: Precipitation values for one season across all years.

    Returns:
        Proportion of time steps with zero precipitation (excluding NaNs).
        Returns NaN if no valid (non-NaN) time steps exist.
    """
    zero_occurrences = season_values[season_values == 0]
    nan_count = int(np.sum(np.isnan(season_values)))
    valid_count = len(season_values) - nan_count

    if valid_count == 0:
        return np.nan
    return len(zero_occurrences) / float(valid_count)


def _gamma_to_standard_normal(
    season_values: np.ndarray,
    zero_probability: float,
    shape_parameter: float,
    scale_parameter: float,
) -> np.ndarray:
    """Transform Gamma-distributed values to standard normal variates.

    Applies the mixed distribution transformation:
        H(x) = zero_probability + (1 - zero_probability) * Gamma_CDF(x)
        SPI   = Phi^-1(H(x))

    where Phi^-1 is the inverse standard normal CDF.

    Args:
        season_values: Accumulated precipitation for one season.
        zero_probability: Probability of zero precipitation.
        shape_parameter: Fitted Gamma shape parameter (alpha).
        scale_parameter: Fitted Gamma scale parameter (beta).

    Returns:
        Standard normal variates (SPI values) for each time step.
        Returns NaN for any time step when parameters are NaN.
    """
    if np.isnan(zero_probability) or np.isnan(shape_parameter) or np.isnan(scale_parameter):
        return np.full_like(season_values, np.nan)

    gamma_cdf_values = (
        zero_probability
        + (1 - zero_probability)
        * gamma.cdf(season_values, shape_parameter, scale=scale_parameter)
    )
    return norm.ppf(gamma_cdf_values)


def spi(
    data: np.ndarray | None = None,
    scale: int | None = None,
    nseas: int | None = None,
    *,
    monthly_values: np.ndarray | None = None,
    accumulation_scale: int | None = None,
    num_seasons: int | None = None,
) -> np.ndarray:
    """Compute the Standardized Precipitation Index for a monthly time series.

    Steps:
        1. Accumulate precipitation over *accumulation_scale* months.
        2. For each season, fit a Gamma distribution to non-zero values.
        3. Apply the mixed-distribution CDF and transform to standard normal.

    **Parameter names:** The function accepts both old (data, scale, nseas) and
    new (monthly_values, accumulation_scale, num_seasons) parameter names for
    backward compatibility. If both are provided, the new names take precedence.

    Args:
        data: Monthly precipitation time series (deprecated, use monthly_values).
            Shape ``(n_months,)``.
        scale: Aggregation window in months (deprecated, use accumulation_scale).
            Common values: 1, 3, 6, 12. A scale of 1 uses raw monthly totals;
            larger scales capture medium- to long-term drought.
        nseas: Number of seasons per year (deprecated, use num_seasons).
            Typically 12 for monthly, 4 for quarterly. Determines which months
            are grouped together when fitting the Gamma distribution.
        monthly_values: Monthly precipitation time series. Shape ``(n_months,)``.
            Preferred over *data*.
        accumulation_scale: Aggregation window in months. Preferred over *scale*.
        num_seasons: Number of seasons per year. Preferred over *nseas*.

    Returns:
        SPI time series of the same length as the accumulated series.
        Negative values indicate drought; positive values indicate surplus.
        NaN values appear when a season has no valid (non-NaN) data.
        Shape ``(n_months - accumulation_scale + 1,)``.

    Example:
        >>> import numpy as np
        >>> from spi import spi
        >>> rng = np.random.default_rng(0)
        >>> precipitation = rng.exponential(scale=80, size=120)  # 10 years
        >>> standardized_index = spi(precipitation, scale=3, nseas=12)
        >>> standardized_index.shape
        (118,)
    """
    # Resolve parameter names: new names override old names for backward compat.
    monthly_values_arr = monthly_values if monthly_values is not None else data
    accumulation_scale_val = accumulation_scale if accumulation_scale is not None else scale
    num_seasons_val = num_seasons if num_seasons is not None else nseas

    if monthly_values_arr is None:
        raise ValueError("Must provide monthly_values (or deprecated data)")
    if accumulation_scale_val is None:
        raise ValueError("Must provide accumulation_scale (or deprecated scale)")
    if num_seasons_val is None:
        raise ValueError("Must provide num_seasons (or deprecated nseas)")

    monthly_values_arr = np.asarray(monthly_values_arr, dtype=float)

    accumulated_values = _accumulate_precipitation(monthly_values_arr, accumulation_scale_val)
    standardized_index = np.zeros(accumulated_values.shape)

    for season in range(num_seasons_val):
        season_indices = np.arange(season, len(accumulated_values), num_seasons_val)
        season_values = accumulated_values[season_indices]

        zero_probability = _compute_zero_probability(season_values)
        clean_values = _remove_zeros_and_nans(season_values)

        # Handle all-missing seasons gracefully: fill SPI with NaNs.
        if clean_values.size == 0 or np.isnan(zero_probability):
            standardized_index[season_indices] = np.nan
            continue

        # Fit Gamma to non-zero, non-NaN values.
        try:
            fitted_params = gamma.fit(clean_values, floc=0)
            shape_parameter = fitted_params[0]
            scale_parameter = fitted_params[2]
        except (ValueError, RuntimeError):
            # gamma.fit can fail on edge cases (e.g., all identical values).
            # Fill season with NaN and move on.
            standardized_index[season_indices] = np.nan
            continue

        standardized_index[season_indices] = _gamma_to_standard_normal(
            season_values, zero_probability, shape_parameter, scale_parameter
        )

    return standardized_index
