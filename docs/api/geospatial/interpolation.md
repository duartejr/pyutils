# Spatial Interpolation

Two complementary interpolation methods for mapping station observations to regular grids.

## InverseDistanceWeighting

IDW assigns weights proportional to the inverse of distance raised to a power. Higher `power` values make the interpolation more local (sharp transitions near stations); lower values give smoother fields.

### Basic usage

```python
import numpy as np
from pyutils.geospatial import InverseDistanceWeighting

# Station positions (longitude, latitude or UTM X, Y)
x_stations = np.array([-45.0, -44.0, -43.0, -44.0])
y_stations = np.array([ -5.0,  -5.0,  -5.0,  -6.0])
z_values   = np.array([  85.0, 120.0,  95.0, 110.0])  # precipitation (mm)

# Target grid
x_grid = np.linspace(-46.0, -42.0, 40)
y_grid = np.linspace( -7.0,  -4.0, 30)

idw = InverseDistanceWeighting(power=2)
grid = idw.interpolate(x_stations, y_stations, z_values, x_grid, y_grid)

print(grid.shape)  # (30, 40)
```

### Plotting the result

```python
import matplotlib.pyplot as plt

fig, ax = plt.subplots(figsize=(8, 5))
im = ax.pcolormesh(x_grid, y_grid, grid, cmap="Blues", shading="auto")
ax.scatter(x_stations, y_stations, c=z_values, cmap="Blues",
           edgecolors="k", s=60, zorder=5)
plt.colorbar(im, ax=ax, label="Precipitation (mm)")
ax.set_title("IDW Interpolated Precipitation")
plt.tight_layout()
plt.show()
```

### Choosing the power parameter

```python
for power in [1, 2, 3]:
    idw = InverseDistanceWeighting(power=power)
    grid = idw.interpolate(x_stations, y_stations, z_values, x_grid, y_grid)
    print(f"power={power}  std={grid.std():.2f}")
# power=1  std=8.1   (smoothest)
# power=2  std=10.4
# power=3  std=12.9  (most local)
```

::: pyutils.geospatial.interpolation.InverseDistanceWeighting
    options:
      show_signature_annotations: true
      members: [interpolate]

---

## ThiessenPolygon

Assigns each grid point the value of the nearest station (Voronoi / nearest-neighbour). Useful for precipitation when spatial smoothing is not desired.

### Basic usage

```python
from pyutils.geospatial import ThiessenPolygon

th = ThiessenPolygon()
grid = th.interpolate(x_stations, y_stations, z_values, x_grid, y_grid)

# Only station values appear in the result
import numpy as np
print(np.unique(grid))  # [85.0, 95.0, 110.0, 120.0]
```

### Polygon areas (areal rainfall weight)

```python
areas = th.compute_polygon_areas(x_stations, y_stations)
weights = areas / areas.sum()

areal_mean = (z_values * weights).sum()
print(f"Area-weighted mean precipitation: {areal_mean:.1f} mm")
```

::: pyutils.geospatial.interpolation.ThiessenPolygon
    options:
      show_signature_annotations: true
      members: [interpolate, compute_polygon_areas]
