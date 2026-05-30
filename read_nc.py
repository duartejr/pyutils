"""
NetCDF4 climate data I/O utilities.

Provides helpers for reading gridded climate variables from NetCDF files,
with optional spatial subsetting to a shapefile bounding box.
"""

from datetime import datetime as dt
from typing import Optional, Tuple

import numpy as np
import geopandas as gpd
from dateutil.relativedelta import relativedelta
from netCDF4 import Dataset


def _in_poly(
    shp: str, lat: np.ndarray, lon: np.ndarray, buffer: float = 0
) -> Tuple[np.ndarray, np.ndarray]:
    """Return grid-index arrays for points inside a shapefile's bounding box.

    Parameters
    ----------
    shp : str
        Path to the shapefile.
    lat, lon : np.ndarray
        1-D arrays of grid latitudes / longitudes.
    buffer : float, optional
        Extra degrees to expand the bounding box. Doubled on each retry if no
        grid points fall inside.

    Returns
    -------
    pos_lon, pos_lat : np.ndarray
        Integer index arrays into *lon* and *lat*.
    """
    gpd_shp = gpd.read_file(shp)["geometry"].set_crs("EPSG:4326", allow_override=True)
    gpd_shp = gpd_shp.buffer(buffer)
    bounds = gpd_shp.bounds
    pos_lon = np.where((lon >= bounds.minx.iloc[0]) & (lon <= bounds.maxx.iloc[0]))[0]
    pos_lat = np.where((lat >= bounds.miny.iloc[0]) & (lat <= bounds.maxy.iloc[0]))[0]
    buffer = buffer * 2 if buffer > 0 else 0.5
    if len(pos_lon) == 0 or len(pos_lat) == 0:
        pos_lon, pos_lat = _in_poly(shp, lat, lon, buffer)
    return pos_lon, pos_lat


def _read_var(
    file: str,
    var_name: str,
    plat: np.ndarray,
    plon: np.ndarray,
    ptime: np.ndarray,
) -> np.ndarray:
    """Read a variable slice from a NetCDF file.

    Handles 3-D (time, lat, lon) and 4-D (level, time, lat, lon) variables.
    """
    ds = Dataset(file)
    n_dim = len(ds.variables[var_name].shape)
    if n_dim == 3:
        return ds.variables[var_name][ptime, plat, plon]
    elif n_dim == 4:
        return ds.variables[var_name][:, ptime, plat, plon]
    return ds.variables[var_name][plat, plon]


def _parse_dates(time_units: str, times: np.ndarray) -> np.ndarray:
    """Convert a NetCDF time axis to an array of datetime objects.

    Parameters
    ----------
    time_units : str
        The ``units`` attribute of the NetCDF time variable (e.g. ``'days
        since 1900-01-01'``).
    times : np.ndarray
        Numeric time values from the NetCDF file.

    Returns
    -------
    np.ndarray of datetime
    """
    parts = time_units.split(" ")
    # Try several date-string positions in the units string
    for candidate in (parts[-1], parts[-2], parts[-1].split("T")[0]):
        try:
            origin = dt.strptime(candidate, "%Y-%m-%d")
            break
        except ValueError:
            continue
    else:
        raise ValueError(f"Cannot parse origin date from units: {time_units!r}")

    freq = parts[0].lower()
    if freq == "months":
        return np.array([origin + relativedelta(months=int(x)) for x in times])
    elif freq[:3] == "day":
        return np.array([origin + relativedelta(days=int(x)) for x in times])
    elif freq == "hours":
        return np.array([origin + relativedelta(hours=int(x)) for x in times])
    raise ValueError(f"Unsupported time frequency: {freq!r}")


def _dates_in(
    dates: np.ndarray, start_date: dt, end_date: dt
) -> Tuple[np.ndarray, np.ndarray]:
    """Return indices and date values within [start_date, end_date]."""
    pos = np.where((dates >= start_date) & (dates <= end_date))[0]
    return pos, dates[pos]


def _get_coords(file: str, variables: list) -> Tuple[np.ndarray, np.ndarray, str, np.ndarray]:
    """Extract lat, lon, and time axis from a NetCDF file."""
    lat = lon = time_units = times = None
    for var in variables:
        vl = var.lower()
        if vl in ("lat", "latitude", "y", "nlat"):
            lat = Dataset(file).variables[var][:]
        elif vl in ("lon", "longitude", "x", "nlon"):
            lon = Dataset(file).variables[var][:]
            lon[lon > 180] -= 360
        elif vl in ("time", "gc_t", "gc", "s"):
            time_units = Dataset(file).variables[var].units
            times = Dataset(file).variables[var][:]
    return lat, lon, time_units, times


def read_nc(
    file: str,
    shape: Optional[str] = None,
    var_name: str = "pre",
    freq: str = "month",
    start_date: Optional[dt] = None,
    end_date: Optional[dt] = None,
    buffer: float = 0,
):
    """Read a variable from a NetCDF file, optionally clipped to a shapefile.

    Parameters
    ----------
    file : str
        Path to the NetCDF file.
    shape : str, optional
        Path to shapefile. When provided, only grid points within the
        shapefile's bounding box are returned.
    var_name : str, optional
        Name of the variable to read. Default ``'pre'``.
    freq : str, optional
        Temporal frequency. When ``'month'`` and the variable units contain
        ``'day'``, values are multiplied by 30 to convert to monthly totals.
    start_date, end_date : datetime, optional
        Date range filter. Both must be provided together.
    buffer : float, optional
        Extra degrees to add to the shapefile bounding box. Default 0.

    Returns
    -------
    var : np.ndarray
        Variable values with shape (time, lat, lon) or (lat, lon).
    dates : np.ndarray of datetime
        Timestamps for each time step (only when a time dimension exists).
    lat_in : np.ndarray
        Latitudes of the returned grid subset.
    lon_in : np.ndarray
        Longitudes of the returned grid subset.

    Notes
    -----
    When the file has no time dimension, the date tuple element is omitted
    and only ``(var, lat_in, lon_in)`` is returned.
    """
    variables = list(Dataset(file).variables)

    if var_name not in variables:
        raise ValueError(f"Variable '{var_name}' not found. Available: {variables}")

    lat, lon, time_units, times = _get_coords(file, variables)

    if time_units is not None:
        dates = _parse_dates(time_units, times)
    else:
        dates = None

    if start_date is not None and end_date is not None and dates is not None:
        pos_time_in, dates = _dates_in(dates, start_date, end_date)
    elif dates is not None and len(dates) > 1:
        pos_time_in = np.arange(len(dates))
    else:
        pos_time_in = np.array([0])

    if shape:
        pos_lon_in, pos_lat_in = _in_poly(shape, lat, lon, buffer)
        lat_in = lat[pos_lat_in]
        lon_in = lon[pos_lon_in]
    else:
        lat_in = lat
        lon_in = lon
        pos_lon_in = np.arange(len(lon))
        pos_lat_in = np.arange(len(lat))

    var = _read_var(file, var_name, pos_lat_in, pos_lon_in, pos_time_in)

    try:
        if freq == "month" and "day" in Dataset(file).variables[var_name].units:
            var = var * 30
    except AttributeError:
        pass

    if dates is not None:
        return var, dates, lat_in, lon_in
    return np.array(var), lat_in, lon_in


def save_nc(
    var: np.ndarray,
    lat: np.ndarray,
    lon: np.ndarray,
    dates: np.ndarray,
    filename: str,
    time_units: str,
    var_name: str,
    var_units: str,
    var_long_name: str,
) -> None:
    """Write a (time, lat, lon) array to a NetCDF4 file.

    Parameters
    ----------
    var : np.ndarray
        Data array with shape (time, lat, lon).
    lat, lon : np.ndarray
        1-D coordinate arrays.
    dates : np.ndarray
        Numeric time values matching *time_units*.
    filename : str
        Output file path.
    time_units : str
        CF-compliant time units string (e.g. ``'months since 1900-01-01'``).
    var_name : str
        Short variable name in the output file.
    var_units : str
        Units of the variable (e.g. ``'mm'``).
    var_long_name : str
        Descriptive long name for the variable.
    """
    output = Dataset(filename, "w", format="NETCDF4")
    output.createDimension("time", len(dates))
    time_var = output.createVariable("time", "f4", "time")
    time_var.units = time_units
    time_var.calendar = "standard"
    time_var[:] = dates

    output.createDimension("lat", len(lat))
    lat_nc = output.createVariable("lat", "f4", "lat")
    lat_nc.units = "degrees_north"
    lat_nc.long_name = "latitude"
    lat_nc.standard_name = "latitude"
    lat_nc[:] = lat

    output.createDimension("lon", len(lon))
    lon_nc = output.createVariable("lon", "f4", "lon")
    lon_nc.units = "degrees_east"
    lon_nc.long_name = "longitude"
    lon_nc.standard_name = "longitude"
    lon_nc[:] = lon

    var_nc = output.createVariable(var_name, "f4", ["time", "lat", "lon"])
    var_nc.units = var_units
    var_nc.long_name = var_long_name
    var_nc.missing_value = np.float32(np.nan)
    var_nc[:] = var

    output.close()
