# NetCDF I/O (xarray)

`XarrayNetCDFHandler` reads and writes climate NetCDF files using [xarray](https://docs.xarray.dev), enabling lazy loading, multi-dimensional operations, and seamless dask integration for large datasets.

## Reading a file

```python
from pyutils.geospatial import XarrayNetCDFHandler

handler = XarrayNetCDFHandler()
ds = handler.read("era5_precip_1991_2020.nc")

print(ds)
# <xarray.Dataset>
# Dimensions: (time: 360, lat: 90, lon: 180)
# Variables: pr (time, lat, lon)
```

## Lazy loading with dask (large files)

```python
handler = XarrayNetCDFHandler(chunks={"time": 60, "lat": 45, "lon": 90})
ds = handler.read("large_climate_file.nc")
# Data is NOT loaded into memory yet — computations are deferred
```

## Compute climatology

```python
# Monthly climatology (mean across all years for each calendar month)
climatology = handler.compute_climatology(group_by="month")

# Seasonal climatology
clim_seasonal = handler.compute_climatology(group_by="season")

# Annual mean
clim_annual = handler.compute_climatology(group_by="year")
```

## Compute anomalies

```python
# Absolute anomaly: value − monthly_mean
anomaly = handler.compute_anomaly(method="absolute")

# Relative anomaly (%): 100 × (value − monthly_mean) / monthly_mean
anomaly_pct = handler.compute_anomaly(method="relative")
```

## Spatial and temporal subsetting

```python
# Clip to Northeast Brazil
northeast = handler.select_by_bounds(
    lat_min=-18.0, lat_max=5.0,
    lon_min=-48.0, lon_max=-34.0,
)

# Select a time period
period = handler.time_slice("1991-01-01", "2000-12-31")
```

## Writing to disk

```python
# Save with compression
handler.write(
    "output.nc",
    encoding={"pr": {"zlib": True, "complevel": 4}},
)
```

::: pyutils.geospatial.io.XarrayNetCDFHandler
    options:
      show_signature_annotations: true
      members: [read, write, compute_climatology, compute_anomaly, select_by_bounds, time_slice]
