# pyutils.geospatial

Spatial data handling — NetCDF climate I/O, spatial interpolation, shapefile operations, and cartographic visualization.

---

## Classes

| Class | Description |
|-------|-------------|
| [`XarrayNetCDFHandler`](io.md#xarraynetcdfhandler) | Read/write NetCDF files with xarray; climatology, anomaly, spatial subsetting |
| [`InverseDistanceWeighting`](interpolation.md#inversedistanceweighting) | IDW spatial interpolation from station data to a regular grid |
| [`ThiessenPolygon`](interpolation.md#thiessenpolygon) | Nearest-station (Thiessen/Voronoi) interpolation |
| [`ShapefileHandler`](shapefile.md#shapefilehandler) | Shapefile I/O, area computation, buffering, clipping, spatial joins |
| [`MapRenderer`](visualization.md#maprenderer) | Thematic map creation with matplotlib; colorbars, titles, export |

---

## Quick Start

### NetCDF Climate Data

```python
from pyutils.geospatial import XarrayNetCDFHandler

handler = XarrayNetCDFHandler()
ds = handler.read("era5_precip_1991_2020.nc")

# Monthly climatology
climatology = handler.compute_climatology(group_by="month")

# Absolute anomaly
anomaly = handler.compute_anomaly()

# Clip to region
northeast = handler.select_by_bounds(
    lat_min=-18.0, lat_max=5.0,
    lon_min=-48.0, lon_max=-34.0,
)
```

### Spatial Interpolation

```python
import numpy as np
from pyutils.geospatial import InverseDistanceWeighting, ThiessenPolygon

x_stations = np.array([-45.0, -44.0, -43.0, -44.0])
y_stations = np.array([ -5.0,  -5.0,  -5.0,  -6.0])
precip     = np.array([  85.0, 120.0,  95.0, 110.0])

x_grid = np.linspace(-46.0, -42.0, 40)
y_grid = np.linspace( -7.0,  -4.0, 30)

# IDW
idw_grid = InverseDistanceWeighting(power=2).interpolate(
    x_stations, y_stations, precip, x_grid, y_grid
)

# Thiessen (nearest station)
th_grid = ThiessenPolygon().interpolate(
    x_stations, y_stations, precip, x_grid, y_grid
)
```

### Shapefile Operations

```python
from pyutils.geospatial import ShapefileHandler
import numpy as np

shp = ShapefileHandler(crs="EPSG:4674")
gdf = shp.read_shapefile("basin.shp")

# Compute area in km²
area_m2 = shp.compute_area()
print(f"Basin area: {area_m2 / 1e6:.1f} km²")

# Buffer by 0.1 degrees
gdf_buffered = shp.buffer_geometries(distance=0.1)

# Points inside polygon
lats = np.linspace(-10.0, -5.0, 50)
lons = np.linspace(-48.0, -44.0, 50)
inside = shp.points_in_polygon(lat=lats, lon=lons)
print(f"{len(inside)} grid points inside basin")
```

### Map Visualization

```python
from pyutils.geospatial import MapRenderer, ShapefileHandler

shp = ShapefileHandler()
gdf = shp.read_shapefile("brazil_states.shp")

renderer = MapRenderer(figsize=(12, 10))
fig, ax = renderer.create_figure()
renderer.plot_geodataframe(gdf, column="population", cmap="YlOrRd")
renderer.set_title("Population by State")
renderer.save("brazil_population.png", dpi=150)
renderer.close()
```

---

## API Reference

- [NetCDF I/O (xarray)](io.md) — `XarrayNetCDFHandler`
- [Interpolation](interpolation.md) — `InverseDistanceWeighting`, `ThiessenPolygon`
- [Shapefile Operations](shapefile.md) — `ShapefileHandler`
- [Visualization](visualization.md) — `MapRenderer`
