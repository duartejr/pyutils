# Shapefile Operations

`ShapefileHandler` wraps geopandas with a clean interface for the most common GIS operations needed in hydroclimate workflows.

## Read and inspect

```python
from pyutils.geospatial import ShapefileHandler

shp = ShapefileHandler()
gdf = shp.read_shapefile("basins/northeast_brazil.shp")

print(gdf.head())
print(f"CRS: {gdf.crs}")
print(f"Features: {len(gdf)}")
```

## Compute polygon areas

```python
# Area in m² using equal-area projection (EPSG:102033 — Albers Americas)
areas_m2 = shp.compute_areas()
areas_km2 = areas_m2 / 1e6

print(f"Basin areas: {areas_km2.round(1)} km²")

# Single polygon
area = shp.compute_area(gdf.geometry.iloc[0])
print(f"First polygon: {area / 1e6:.2f} km²")
```

## Create buffers

```python
# 10 km buffer around each polygon (units depend on CRS)
gdf_buffered = shp.buffer_geometries(distance=10000)
shp.write_shapefile("basins_buffered.shp", gdf_buffered)
```

## Reproject

```python
# Convert from SIRGAS 2000 geographic to UTM zone 24S
gdf_utm = shp.reproject("EPSG:31984")
```

## Spatial join — points to polygons

```python
import geopandas as gpd

stations = gpd.read_file("stations.shp")
basins   = gpd.read_file("basins.shp")

# Which basin does each station fall in?
joined = shp.spatial_join(stations, basins, how="left", predicate="within")
print(joined[["station_id", "basin_name"]])
```

## Dissolve

```python
# Merge all sub-basins into single geometry per region
regions = shp.dissolve(by="region_id", gdf=gdf)
```

## Clip raster extent to watershed

```python
from shapely.ops import unary_union

watershed = unary_union(gdf.geometry)
gdf_clipped = shp.clip(stations, watershed)
```

::: pyutils.geospatial.shapefile.ShapefileHandler
    options:
      show_signature_annotations: true
