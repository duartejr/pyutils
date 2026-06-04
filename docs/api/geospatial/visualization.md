# Visualization

`MapRenderer` wraps matplotlib and geopandas for creating thematic maps with custom colorbars — the most common cartographic output in hydroclimate research.

## Choropleth map from GeoDataFrame

```python
import geopandas as gpd
from pyutils.geospatial import MapRenderer

gdf = gpd.read_file("basins.shp")
gdf["annual_precip_mm"] = [1200, 850, 600, 1450, 920]

renderer = MapRenderer(figsize=(10, 8))
renderer.create_figure()
renderer.plot_geodataframe(gdf, column="annual_precip_mm", cmap="Blues")
renderer.set_title("Annual Precipitation by Basin")
renderer.set_labels(xlabel="Longitude", ylabel="Latitude")
renderer.save("precip_map.png", dpi=300)
renderer.show()
```

## Add a standard colorbar

```python
renderer.create_figure()
renderer.plot_geodataframe(gdf, column="annual_precip_mm", cmap="Blues", legend=False)
renderer.add_colorbar(
    cmap="Blues",
    bounds=[400, 600, 800, 1000, 1200, 1400, 1600],
    label="Precipitation (mm/year)",
    orientation="vertical",
    extend="both",
)
renderer.save("precip_colorbar.png")
```

## Discrete custom colorbar (legacy style)

```python
import matplotlib.pyplot as plt

renderer.create_figure()
renderer.plot_geodataframe(gdf, column="annual_precip_mm", cmap="RdBu", legend=False)
renderer.customize_colorbar(
    cmap=plt.get_cmap("RdBu"),
    bounds=[-200, -100, 0, 100, 200],
    label="Anomaly (mm)",
    ticks_label=["-200", "-100", "Normal", "+100", "+200"],
    orientation="vertical",
)
renderer.save("anomaly_map.png")
renderer.close()
```

## Clean map without axes

```python
renderer.create_figure()
renderer.plot_geodataframe(gdf, cmap="terrain")
renderer.remove_axis()   # no tick marks or spines
renderer.save("clean_map.png")
renderer.close()
```

::: pyutils.geospatial.visualization.MapRenderer
    options:
      show_signature_annotations: true
