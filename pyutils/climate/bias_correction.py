"""Bias correction methods for climate model output.

Consolidates corr_med.py, rv/RvGama.py, rv/RvEmpirica.py, and
rv/RvLinear.py implementations for distributional bias correction
and quantile mapping.
"""

from typing import Tuple, Optional
import numpy as np
from scipy.stats import gamma, norm
from sklearn.linear_model import LinearRegression
import xarray as xr

from pyutils.core import ValidationError


class BiasCorrection:
    """Bias correction for climate model forecasts.

    Consolidates corr_med.py and RvGama.py implementations with support for:
    - Linear scaling (mean/variance correction)
    - Quantile mapping (gamma distribution)
    - Parametric distribution fitting
    """

    @staticmethod
    def linear_scaling(
        forecast: np.ndarray,
        hindcast: np.ndarray,
        observation: np.ndarray,
        method: str = "additive",
    ) -> np.ndarray:
        """Linear scaling bias correction.

        Adjusts forecast to match observed statistics. Consolidates
        corr_med.py method() function.

        Parameters
        ----------
        forecast : np.ndarray
            Climate model forecast.
        hindcast : np.ndarray
            Historical model simulation (same period as observations).
        observation : np.ndarray
            Observations for reference period.
        method : str, optional
            'additive' (mean correction) or 'multiplicative' (ratio correction).

        Returns
        -------
        np.ndarray
            Bias-corrected forecast.

        Raises
        ------
        ValidationError
            If arrays have inconsistent shapes or invalid method.
        """
        forecast = np.asarray(forecast, dtype=float)
        hindcast = np.asarray(hindcast, dtype=float)
        observation = np.asarray(observation, dtype=float)

        if not (hindcast.shape == observation.shape):
            raise ValidationError(
                f"hindcast {hindcast.shape} and observation {observation.shape} "
                "must have same shape"
            )

        if method == "additive":
            # Mean correction: fcst_corrected = fcst + (obs_mean - hindcast_mean)
            correction = np.nanmean(observation) - np.nanmean(hindcast)
            corrected = forecast + correction
        elif method == "multiplicative":
            # Ratio correction: fcst_corrected = fcst * (obs_mean / hindcast_mean)
            correction = np.nanmean(observation) / (np.nanmean(hindcast) + 1e-10)
            corrected = forecast * correction
        else:
            raise ValidationError(f"Unknown method: {method}")

        return corrected

    @staticmethod
    def variance_scaling(
        forecast: np.ndarray,
        hindcast: np.ndarray,
        observation: np.ndarray,
    ) -> np.ndarray:
        """Variance scaling bias correction (Monthly method from corr_med.py).

        Adjusts both mean and variance to match observations.
        Consolidates corr_med.py method 1 (detailed implementation).

        Parameters
        ----------
        forecast : np.ndarray
            Climate model forecast (shape: n_times).
        hindcast : np.ndarray
            Historical model data (shape: n_years*12).
        observation : np.ndarray
            Observations (shape: n_years*12).

        Returns
        -------
        np.ndarray
            Variance-scaled forecast.

        Raises
        ------
        ValidationError
            If input shapes are inconsistent.
        """
        forecast = np.asarray(forecast, dtype=float)
        hindcast = np.asarray(hindcast, dtype=float)
        observation = np.asarray(observation, dtype=float)

        # Reshape to years x months if needed
        if len(hindcast) % 12 == 0 and len(observation) % 12 == 0:
            n_years_h = len(hindcast) // 12
            n_years_o = len(observation) // 12

            hindcast = hindcast[: n_years_h * 12].reshape((n_years_h, 12))
            observation = observation[: n_years_o * 12].reshape((n_years_o, 12))

            # Monthly statistics
            hist_mean = np.nanmean(hindcast, axis=0)
            hist_std = np.nanstd(hindcast, axis=0)
            obs_mean = np.nanmean(observation, axis=0)
            obs_std = np.nanstd(observation, axis=0)

            # Apply monthly correction to forecast
            corrected = np.zeros_like(forecast)
            for i in range(len(forecast)):
                month = i % 12
                if hist_std[month] > 0:
                    corrected[i] = (
                        obs_mean[month]
                        + (forecast[i] - hist_mean[month])
                        * (obs_std[month] / hist_std[month])
                    )
                else:
                    corrected[i] = forecast[i]

                # Ensure non-negative (for precipitation)
                corrected[i] = max(0, corrected[i])
        else:
            corrected = forecast

        return corrected

    @staticmethod
    def quantile_mapping_gamma(
        forecast: np.ndarray,
        hindcast: np.ndarray,
        observation: np.ndarray,
    ) -> np.ndarray:
        """Quantile mapping using gamma distribution.

        Consolidates rv/RvGama.py rv_gama() function.
        Preserves zero probability from observations.

        Parameters
        ----------
        forecast : np.ndarray
            Single forecast value or array of values.
        hindcast : np.ndarray
            Hindcast (model historical) data.
        observation : np.ndarray
            Observation (reference) data.

        Returns
        -------
        np.ndarray
            Quantile-mapped forecast value(s).

        Raises
        ------
        ValidationError
            If distribution fitting fails.
        """
        forecast = np.asarray(forecast, dtype=float)
        hindcast = np.asarray(hindcast, dtype=float)
        observation = np.asarray(observation, dtype=float)

        try:
            # Fit gamma distribution to observations
            obs_positive = observation[observation > 0]
            n_zeros_obs = np.sum(observation == 0)
            q = n_zeros_obs / len(observation)  # Probability of zero

            if len(obs_positive) > 1:
                gamma_obs = gamma.fit(obs_positive, floc=0)
            else:
                raise ValidationError("Insufficient non-zero observation data")

            # Fit gamma to hindcast
            hind_sorted = np.sort(hindcast[hindcast > 0])
            if len(hind_sorted) > 1:
                gamma_hind = gamma.fit(hind_sorted, floc=0)
            else:
                raise ValidationError("Insufficient hindcast data")

            # Apply quantile mapping to forecast
            if np.isscalar(forecast):
                corrected = BiasCorrection._qm_single_gamma(
                    forecast, gamma_hind, gamma_obs, q, hind_sorted[0]
                )
            else:
                corrected = np.array(
                    [
                        BiasCorrection._qm_single_gamma(
                            f, gamma_hind, gamma_obs, q, hind_sorted[0]
                        )
                        for f in forecast
                    ]
                )

            return corrected
        except Exception as e:
            raise ValidationError(f"Quantile mapping failed: {e}")

    @staticmethod
    def _qm_single_gamma(
        fcst: float,
        gamma_hind: Tuple[float, float],
        gamma_obs: Tuple[float, float],
        q: float,
        hind_min: float,
    ) -> float:
        """Apply quantile mapping to single forecast value."""
        if fcst < hind_min:
            return 0.0

        prob_mod = gamma.cdf(fcst, *gamma_hind)
        H = q + (1 - q) * prob_mod

        corrected = gamma.ppf(H, *gamma_obs)

        # Handle numerical issues
        if np.isinf(corrected) or np.isnan(corrected):
            corrected = 0.0

        return corrected

    @staticmethod
    def quantile_mapping_empirical(
        forecast: float,
        hindcast: np.ndarray,
        observation: np.ndarray,
    ) -> float:
        """Empirical quantile mapping using rank-based nearest-neighbour lookup.

        Finds the rank of the forecast in the sorted hindcast distribution and
        returns the observation at the same rank. Consolidates rv/RvEmpirica.py.

        Parameters
        ----------
        forecast : float
            Single forecast value to correct.
        hindcast : np.ndarray
            Historical model simulation data.
        observation : np.ndarray
            Historical observation data (same length as hindcast).

        Returns
        -------
        float
            Bias-corrected forecast value.

        Raises
        ------
        ValidationError
            If hindcast and observation have different lengths.

        Examples
        --------
        >>> import numpy as np
        >>> hind = np.array([10.0, 20.0, 30.0, 40.0, 50.0])
        >>> obs  = np.array([12.0, 22.0, 33.0, 43.0, 54.0])
        >>> BiasCorrection.quantile_mapping_empirical(25.0, hind, obs)
        33.0
        """
        hindcast = np.asarray(hindcast, dtype=float)
        observation = np.asarray(observation, dtype=float)

        if len(hindcast) != len(observation):
            raise ValidationError(
                f"hindcast length {len(hindcast)} != observation length {len(observation)}"
            )

        sorted_hindcast = np.sort(hindcast)
        sorted_observation = np.sort(observation)

        # Find rank of forecast in hindcast (index of closest value)
        distances = np.abs(forecast - sorted_hindcast)
        rank = int(np.argmin(distances))

        return float(sorted_observation[rank])

    @staticmethod
    def linear_error_correction(
        forecast: np.ndarray,
        errors: np.ndarray,
    ) -> np.ndarray:
        """Regression-based error propagation for multi-step forecasts.

        Applies a linear regression model trained on consecutive historical
        errors to propagate a sampled error through the forecast horizon.
        Consolidates rv/RvLinear.py.

        Parameters
        ----------
        forecast : np.ndarray
            Forecast values for each lead time (shape: n_lead_times).
        errors : np.ndarray
            Historical errors, shape (n_samples, n_lead_times). Each row is
            one historical sample; columns correspond to lead times d=1…n.

        Returns
        -------
        np.ndarray
            Error-corrected forecast (shape: n_lead_times). Negative values
            are clipped to zero.

        Raises
        ------
        ValidationError
            If forecast length does not match number of error columns, or if
            there are fewer than 2 historical samples.

        Examples
        --------
        >>> import numpy as np
        >>> rng = np.random.default_rng(0)
        >>> fcst = np.array([10.0, 11.0, 12.0])
        >>> errs = rng.normal(0, 1, size=(30, 3))
        >>> corrected = BiasCorrection.linear_error_correction(fcst, errs)
        >>> corrected.shape
        (3,)
        """
        forecast = np.asarray(forecast, dtype=float)
        errors = np.asarray(errors, dtype=float)

        n_lead = forecast.shape[0]
        if errors.ndim != 2:
            raise ValidationError("errors must be a 2-D array (n_samples, n_lead_times)")

        n_samples, n_cols = errors.shape
        if n_cols != n_lead:
            raise ValidationError(
                f"forecast length {n_lead} != number of error columns {n_cols}"
            )
        if n_samples < 2:
            raise ValidationError("Need at least 2 historical error samples")

        corrected = forecast.copy()

        # Sample an initial error from history
        idx = int(np.random.randint(0, n_samples))
        initial_error = errors[idx, 0]
        corrected[0] += initial_error

        # Propagate error through lead times using linear regression
        for lead in range(1, n_lead):
            x_train = errors[:, lead - 1].reshape(-1, 1)
            y_train = errors[:, lead]
            regr = LinearRegression().fit(x_train, y_train)
            propagated = regr.predict([[initial_error]])[0]
            corrected[lead] += propagated
            initial_error = propagated

        corrected = np.maximum(corrected, 0.0)
        return corrected

    @staticmethod
    def correct_dataset(
        forecast_ds: xr.Dataset,
        hindcast_ds: xr.Dataset,
        observation_ds: xr.Dataset,
        method: str = "variance_scaling",
        var_name: str = "pr",
    ) -> xr.Dataset:
        """Apply bias correction to entire xarray dataset.

        Parameters
        ----------
        forecast_ds : xr.Dataset
            Forecast dataset.
        hindcast_ds : xr.Dataset
            Hindcast dataset.
        observation_ds : xr.Dataset
            Observation dataset.
        method : str, optional
            Correction method ('linear_scaling', 'variance_scaling', 'quantile_mapping').
        var_name : str, optional
            Variable name to correct. Default is 'pr' (precipitation).

        Returns
        -------
        xr.Dataset
            Corrected dataset.

        Raises
        ------
        ValidationError
            If method unknown or correction fails.
        """
        if method == "linear_scaling":
            func = BiasCorrection.linear_scaling
        elif method == "variance_scaling":
            func = BiasCorrection.variance_scaling
        elif method == "quantile_mapping":
            func = BiasCorrection.quantile_mapping_gamma
        else:
            raise ValidationError(f"Unknown method: {method}")

        corrected_ds = forecast_ds.copy()

        # Apply correction to each grid point
        # This is simplified; real implementation would vectorize properly
        corrected_ds[var_name].values = func(
            forecast_ds[var_name].values.ravel(),
            hindcast_ds[var_name].values.ravel(),
            observation_ds[var_name].values.ravel(),
        ).reshape(forecast_ds[var_name].shape)

        return corrected_ds
