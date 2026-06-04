"""Data validation utilities."""

from typing import Union, Tuple, Optional
import numpy as np
from pyutils.core import ValidationError, constants


class DataValidator:
    """Centralized validation for common data patterns."""

    @staticmethod
    def validate_precipitation(
        data: np.ndarray,
        allow_zero: bool = True,
        allow_missing: bool = False,
    ) -> bool:
        """Validate precipitation data.

        Parameters
        ----------
        data : np.ndarray
            Precipitation values (mm).
        allow_zero : bool, optional
            Allow zero values. Default is True.
        allow_missing : bool, optional
            Allow NaN values. Default is False.

        Returns
        -------
        bool
            True if valid.

        Raises
        ------
        ValidationError
            If validation fails.
        """
        data = np.asarray(data)

        if not allow_missing and np.any(np.isnan(data)):
            raise ValidationError("Missing values (NaN) found")

        valid_data = data[~np.isnan(data)]

        if np.any(valid_data < 0):
            raise ValidationError("Negative precipitation values found")

        if not allow_zero and np.any(valid_data == 0):
            raise ValidationError("Zero precipitation values found")

        if np.any(valid_data > constants.MAX_VALID_PRECIPITATION):
            raise ValidationError(
                f"Precipitation exceeds max ({constants.MAX_VALID_PRECIPITATION} mm)"
            )

        return True

    @staticmethod
    def validate_temperature(
        data: np.ndarray,
        allow_missing: bool = False,
    ) -> bool:
        """Validate temperature data.

        Parameters
        ----------
        data : np.ndarray
            Temperature values (°C).
        allow_missing : bool, optional
            Allow NaN values.

        Returns
        -------
        bool
            True if valid.

        Raises
        ------
        ValidationError
            If validation fails.
        """
        data = np.asarray(data)

        if not allow_missing and np.any(np.isnan(data)):
            raise ValidationError("Missing values (NaN) found")

        valid_data = data[~np.isnan(data)]

        if np.any(valid_data < constants.MIN_VALID_TEMPERATURE):
            raise ValidationError(
                f"Temperature below minimum ({constants.MIN_VALID_TEMPERATURE}°C)"
            )

        if np.any(valid_data > constants.MAX_VALID_TEMPERATURE):
            raise ValidationError(
                f"Temperature exceeds maximum ({constants.MAX_VALID_TEMPERATURE}°C)"
            )

        return True

    @staticmethod
    def validate_relative_humidity(data: np.ndarray) -> bool:
        """Validate relative humidity (0-100%).

        Parameters
        ----------
        data : np.ndarray
            Relative humidity (%).

        Returns
        -------
        bool
            True if valid.

        Raises
        ------
        ValidationError
            If validation fails.
        """
        data = np.asarray(data)
        valid_data = data[~np.isnan(data)]

        if np.any(valid_data < 0) or np.any(valid_data > 100):
            raise ValidationError("Relative humidity must be 0-100%")

        return True

    @staticmethod
    def validate_wind_speed(data: np.ndarray) -> bool:
        """Validate wind speed data.

        Parameters
        ----------
        data : np.ndarray
            Wind speed (m/s).

        Returns
        -------
        bool
            True if valid.

        Raises
        ------
        ValidationError
            If validation fails.
        """
        data = np.asarray(data)
        valid_data = data[~np.isnan(data)]

        if np.any(valid_data < 0):
            raise ValidationError("Wind speed cannot be negative")

        if np.any(valid_data > 50):  # Reasonable upper bound
            raise ValidationError("Wind speed unreasonably high (> 50 m/s)")

        return True

    @staticmethod
    def validate_spatial_data(
        x: np.ndarray,
        y: np.ndarray,
        z: np.ndarray,
    ) -> bool:
        """Validate spatial data consistency.

        Parameters
        ----------
        x, y, z : np.ndarray
            Spatial coordinates and values.

        Returns
        -------
        bool
            True if valid.

        Raises
        ------
        ValidationError
            If validation fails.
        """
        x = np.asarray(x)
        y = np.asarray(y)
        z = np.asarray(z)

        if not (len(x) == len(y) == len(z)):
            raise ValidationError("x, y, z must have same length")

        if len(x) < 2:
            raise ValidationError("Need at least 2 points for interpolation")

        return True

    @staticmethod
    def validate_time_series_length(
        data: np.ndarray,
        min_length: int = 12,
    ) -> bool:
        """Validate time series has minimum length.

        Parameters
        ----------
        data : np.ndarray
            Time series data.
        min_length : int, optional
            Minimum required length.

        Returns
        -------
        bool
            True if valid.

        Raises
        ------
        ValidationError
            If series too short.
        """
        if len(data) < min_length:
            raise ValidationError(
                f"Time series too short: {len(data)} < {min_length}"
            )
        return True
