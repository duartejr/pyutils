"""Hydrological indices and drought indicators."""

from typing import Tuple, Optional
import numpy as np
from scipy.stats import gamma
from pyutils.core import ValidationError


class StandardizedPrecipitationIndex:
    """Standardized Precipitation Index (SPI) for drought assessment.

    Quantifies precipitation deficit on multiple timescales, enabling
    comparison across regions with different precipitation regimes.
    Consolidates previous spi.py implementation with improved numerical stability.

    Reference
    ---------
    McKee, T.B., N.J. Doesken, and J. Kleist. 1993. The relationship of
    drought frequency and duration to time scales. Proceedings of the 8th
    Conference on Applied Climatology.
    """

    def __init__(self, timescale: int = 3, distribution: str = "gamma"):
        """Initialize SPI calculator.

        Parameters
        ----------
        timescale : int
            Accumulation period in months (1, 3, 6, 12, etc.).
        distribution : str, optional
            Probability distribution ('gamma' or 'lognormal'). Default is gamma.

        Raises
        ------
        ValidationError
            If parameters are invalid.
        """
        if timescale < 1:
            raise ValidationError("Timescale must be >= 1")
        if distribution not in {"gamma", "lognormal"}:
            raise ValidationError(f"Unknown distribution: {distribution}")

        self.timescale = timescale
        self.distribution = distribution

    def compute(
        self, precipitation: np.ndarray, climatology_period: Optional[Tuple[int, int]] = None
    ) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
        """Compute SPI values.

        Parameters
        ----------
        precipitation : np.ndarray
            Monthly precipitation series (mm). Should be long enough to
            establish robust statistics (minimum 30 years recommended).
        climatology_period : tuple of int, optional
            (start_index, end_index) for climatology period. If None, uses entire series.

        Returns
        -------
        tuple
            - spi_values: SPI values (dimensionless)
            - accumulated: Accumulated precipitation (mm)
            - probabilities: Empirical cumulative probabilities

        Raises
        ------
        ValidationError
            If input data is invalid.
        """
        precipitation = np.asarray(precipitation, dtype=float)

        if len(precipitation) < self.timescale:
            raise ValidationError(
                f"Data length ({len(precipitation)}) must be >= timescale ({self.timescale})"
            )

        if not np.all(precipitation >= 0):
            raise ValidationError("Precipitation must be non-negative")

        # Accumulate precipitation
        accumulated = self._accumulate_precipitation(precipitation)

        # Establish climatology
        if climatology_period is None:
            climatology_period = (0, len(accumulated))

        start_idx, end_idx = climatology_period
        climatology_data = accumulated[start_idx:end_idx]

        # Fit distribution
        if self.distribution == "gamma":
            shape, scale = self._fit_gamma(climatology_data)
            probabilities = self._gamma_cdf(accumulated, shape, scale)
        else:  # lognormal
            mu, sigma = self._fit_lognormal(climatology_data)
            probabilities = self._lognormal_cdf(accumulated, mu, sigma)

        # Convert probabilities to SPI (standard normal quantiles)
        spi_values = self._probability_to_spi(probabilities)

        return spi_values, accumulated, probabilities

    def _accumulate_precipitation(self, precipitation: np.ndarray) -> np.ndarray:
        """Accumulate precipitation over timescale."""
        n = len(precipitation)
        accumulated = np.convolve(
            precipitation, np.ones(self.timescale), mode="valid"
        )

        # Pad with NaN for initial values
        accumulated = np.concatenate(
            [np.full(self.timescale - 1, np.nan), accumulated]
        )

        return accumulated

    def _fit_gamma(self, data: np.ndarray) -> Tuple[float, float]:
        """Fit gamma distribution to precipitation data."""
        # Filter out zero/NaN values
        valid_data = data[(~np.isnan(data)) & (data > 0)]

        if len(valid_data) < 10:
            raise ValidationError(
                "Insufficient valid data for gamma fitting"
            )

        # Use method of moments
        mean = np.mean(valid_data)
        var = np.var(valid_data)

        shape = mean**2 / var
        scale = var / mean

        return shape, scale

    def _fit_lognormal(self, data: np.ndarray) -> Tuple[float, float]:
        """Fit lognormal distribution to precipitation data."""
        valid_data = data[(~np.isnan(data)) & (data > 0)]

        if len(valid_data) < 10:
            raise ValidationError(
                "Insufficient valid data for lognormal fitting"
            )

        log_data = np.log(valid_data)
        mu = np.mean(log_data)
        sigma = np.std(log_data)

        return mu, sigma

    def _gamma_cdf(
        self, data: np.ndarray, shape: float, scale: float
    ) -> np.ndarray:
        """Compute gamma CDF."""
        probabilities = np.zeros_like(data, dtype=float)

        for i, val in enumerate(data):
            if np.isnan(val) or val == 0:
                probabilities[i] = np.nan
            else:
                probabilities[i] = gamma.cdf(val, shape, scale=scale)

        return probabilities

    def _lognormal_cdf(
        self, data: np.ndarray, mu: float, sigma: float
    ) -> np.ndarray:
        """Compute lognormal CDF."""
        from scipy.stats import lognorm

        probabilities = np.zeros_like(data, dtype=float)

        for i, val in enumerate(data):
            if np.isnan(val) or val == 0:
                probabilities[i] = np.nan
            else:
                probabilities[i] = lognorm.cdf(val, sigma, scale=np.exp(mu))

        return probabilities

    def _probability_to_spi(self, probabilities: np.ndarray) -> np.ndarray:
        """Convert probabilities to SPI values using inverse normal."""
        from scipy.stats import norm

        spi_values = np.zeros_like(probabilities)

        for i, prob in enumerate(probabilities):
            if np.isnan(prob):
                spi_values[i] = np.nan
            else:
                # Clamp to avoid extreme values
                prob = np.clip(prob, 0.0001, 0.9999)
                spi_values[i] = norm.ppf(prob)

        return spi_values
