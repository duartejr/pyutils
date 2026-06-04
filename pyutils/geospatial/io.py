"""NetCDF I/O and climate data handling using xarray.

Consolidates read_nc.py functionality with modern xarray interface.
Enables efficient handling of multi-dimensional climate data.
"""

from typing import Optional, Dict, Any, List
import numpy as np
import xarray as xr
from pathlib import Path
from pyutils.core import ClimateDataHandler, ValidationError


class XarrayNetCDFHandler(ClimateDataHandler):
    """Handle NetCDF files using xarray for efficient N-dimensional operations.

    Replaces the previous netCDF4-based approach with xarray for:
    - Better handling of multi-dimensional data
    - Built-in coordinate management
    - Efficient lazy loading and chunked processing
    - Seamless integration with dask for large files

    Examples
    --------
    >>> handler = XarrayNetCDFHandler()
    >>> ds = handler.read('climate_data.nc')
    >>> climatology = handler.compute_climatology()
    >>> anomaly = handler.compute_anomaly(climatology)
    """

    def __init__(self, name: str = "xarray_handler", chunks: Optional[Dict[str, int]] = None):
        """Initialize xarray NetCDF handler.

        Parameters
        ----------
        name : str, optional
            Handler identifier.
        chunks : dict, optional
            Dask chunk sizes for lazy loading. Example: {'time': 100, 'lat': 50, 'lon': 50}
        """
        super().__init__(
            name=name,
            description="xarray-based NetCDF handler for climate data",
        )
        self.chunks = chunks or {}

    def read(self, path: str) -> xr.Dataset:
        """Read NetCDF file using xarray.

        Parameters
        ----------
        path : str
            Path to NetCDF file.

        Returns
        -------
        xr.Dataset
            Loaded dataset with all variables and coordinates.

        Raises
        ------
        ValidationError
            If file cannot be read or is invalid.
        """
        path = Path(path)
        if not path.exists():
            raise ValidationError(f"File not found: {path}")

        try:
            if self.chunks:
                self._dataset = xr.open_dataset(path, chunks=self.chunks)
            else:
                self._dataset = xr.open_dataset(path)
            return self._dataset
        except Exception as e:
            raise ValidationError(f"Failed to read NetCDF file: {e}")

    def write(
        self,
        path: str,
        dataset: Optional[xr.Dataset] = None,
        encoding: Optional[Dict[str, Any]] = None,
    ) -> None:
        """Write dataset to NetCDF file.

        Parameters
        ----------
        path : str
            Output file path.
        dataset : xr.Dataset, optional
            Dataset to write. If None, uses current dataset.
        encoding : dict, optional
            Variable encoding (e.g., compression, dtype). Example:
            {'temperature': {'zlib': True, 'complevel': 4}}

        Raises
        ------
        ValidationError
            If dataset is None or write fails.
        """
        if dataset is None:
            if self._dataset is None:
                raise ValidationError("No dataset loaded or provided")
            dataset = self._dataset

        path = Path(path)
        path.parent.mkdir(parents=True, exist_ok=True)

        try:
            dataset.to_netcdf(path, encoding=encoding)
        except Exception as e:
            raise ValidationError(f"Failed to write NetCDF file: {e}")

    def compute_climatology(
        self,
        time_dim: str = "time",
        group_by: str = "month",
    ) -> xr.Dataset:
        """Compute climatological statistics (mean, std, etc.).

        Parameters
        ----------
        time_dim : str, optional
            Name of time dimension.
        group_by : str, optional
            Grouping method: 'month', 'season', 'year'.

        Returns
        -------
        xr.Dataset
            Dataset with climatological statistics.

        Raises
        ------
        ValidationError
            If dataset not loaded.
        """
        if self._dataset is None:
            raise ValidationError("No dataset loaded")

        try:
            if group_by == "month":
                # Group by month across all years
                grouper = self._dataset[time_dim].dt.month
                climatology = self._dataset.groupby(grouper).mean()
                climatology = climatology.rename({time_dim: "month"})
            elif group_by == "season":
                # Group by season (DJF, MAM, JJA, SON)
                grouper = self._dataset[time_dim].dt.season
                climatology = self._dataset.groupby(grouper).mean()
            elif group_by == "year":
                # Annual climatology
                climatology = self._dataset.mean(dim=time_dim)
            else:
                raise ValidationError(f"Unknown grouping method: {group_by}")

            # Add statistics
            climatology_std = self._dataset.groupby(grouper).std()
            climatology_min = self._dataset.groupby(grouper).min()
            climatology_max = self._dataset.groupby(grouper).max()

            climatology["std"] = climatology_std
            climatology["min"] = climatology_min
            climatology["max"] = climatology_max

            return climatology
        except Exception as e:
            raise ValidationError(f"Failed to compute climatology: {e}")

    def compute_anomaly(
        self,
        climatology: Optional[xr.Dataset] = None,
        time_dim: str = "time",
        method: str = "absolute",
    ) -> xr.Dataset:
        """Compute anomalies relative to climatology.

        Parameters
        ----------
        climatology : xr.Dataset, optional
            Climatological reference. If None, computes internally.
        time_dim : str, optional
            Name of time dimension.
        method : str, optional
            'absolute' (value - climatology) or 'relative' (100 * (value - climatology) / climatology)

        Returns
        -------
        xr.Dataset
            Dataset with anomalies.

        Raises
        ------
        ValidationError
            If dataset not loaded or computation fails.
        """
        if self._dataset is None:
            raise ValidationError("No dataset loaded")

        if climatology is None:
            climatology = self.compute_climatology()

        try:
            # Add month field if using monthly climatology
            if "month" in climatology.dims:
                months = self._dataset[time_dim].dt.month
                clim_aligned = climatology.sel(month=months)
            else:
                clim_aligned = climatology

            # Align dimensions for broadcasting
            clim_aligned = clim_aligned.drop_vars(
                [v for v in clim_aligned.dims if v not in self._dataset.dims],
                errors="ignore"
            )

            if method == "absolute":
                anomaly = self._dataset - clim_aligned
            elif method == "relative":
                anomaly = 100.0 * (self._dataset - clim_aligned) / clim_aligned
            else:
                raise ValidationError(f"Unknown anomaly method: {method}")

            return anomaly
        except Exception as e:
            raise ValidationError(f"Failed to compute anomaly: {e}")

    def select_by_bounds(
        self,
        lat_min: float,
        lat_max: float,
        lon_min: float,
        lon_max: float,
        lat_dim: str = "lat",
        lon_dim: str = "lon",
    ) -> xr.Dataset:
        """Select spatial region by bounding box.

        Parameters
        ----------
        lat_min, lat_max : float
            Latitude range.
        lon_min, lon_max : float
            Longitude range.
        lat_dim : str, optional
            Name of latitude dimension.
        lon_dim : str, optional
            Name of longitude dimension.

        Returns
        -------
        xr.Dataset
            Subset dataset within bounds.

        Raises
        ------
        ValidationError
            If dataset not loaded or bounds invalid.
        """
        if self._dataset is None:
            raise ValidationError("No dataset loaded")

        try:
            subset = self._dataset.sel(
                {
                    lat_dim: slice(lat_min, lat_max),
                    lon_dim: slice(lon_min, lon_max),
                }
            )
            return subset
        except Exception as e:
            raise ValidationError(f"Failed to select bounds: {e}")

    def time_slice(
        self,
        start_date: str,
        end_date: str,
        time_dim: str = "time",
    ) -> xr.Dataset:
        """Select time period.

        Parameters
        ----------
        start_date : str
            Start date (ISO format: YYYY-MM-DD).
        end_date : str
            End date (ISO format: YYYY-MM-DD).
        time_dim : str, optional
            Name of time dimension.

        Returns
        -------
        xr.Dataset
            Subset for time period.

        Raises
        ------
        ValidationError
            If dataset not loaded or dates invalid.
        """
        if self._dataset is None:
            raise ValidationError("No dataset loaded")

        try:
            subset = self._dataset.sel(
                {time_dim: slice(start_date, end_date)}
            )
            return subset
        except Exception as e:
            raise ValidationError(f"Failed to select time period: {e}")

    def validate(self) -> bool:
        """Validate that dataset is loaded and valid."""
        return self._dataset is not None and isinstance(self._dataset, xr.Dataset)
