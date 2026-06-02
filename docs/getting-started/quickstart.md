# Quick Start

This page walks through the most common workflows in pyutils.

---

## Evapotranspiration

### Hargreaves-Samani (monthly)

```python
import numpy as np
from evapotranspiration import hargreaves

# Monthly mean/max/min temperatures (°C) for one year
mean_temps = np.array([25.0, 25.5, 26.0, 25.5, 25.0, 24.5,
                       24.0, 24.0, 24.5, 25.0, 25.5, 25.5])
max_temps = mean_temps + 5.0
min_temps = mean_temps - 5.0

# Latitude -15° (Southern Hemisphere), monthly aggregation
etp_monthly = hargreaves(
    mean_temps, max_temps, min_temps,
    y=-15,
    months=np.arange(1, 13),
    data_freq="monthly",
)
print(etp_monthly)  # mm/month for each of the 12 months
```

### Thornthwaite (monthly)

```python
from evapotranspiration import thornthwaite

sunshine_hours = np.array([11, 11, 12, 12, 11, 10, 10, 11, 12, 12, 12, 11])

etp = thornthwaite(mean_temps, sunshine_hours)
print(etp)  # mm/month
```

---

## Standardized Precipitation Index (SPI)

```python
import numpy as np
from spi import spi

# Generate a synthetic 10-year monthly precipitation series
rng = np.random.default_rng(seed=0)
precipitation = rng.exponential(scale=80, size=120)  # mm/month

# 3-month SPI (captures short-term drought)
drought_index = spi(precipitation, accumulation_scale=3, num_seasons=12)
print(drought_index.shape)    # (118,)
print(drought_index[:12])     # First year's 3-month SPI values

# SPI < -1.5 → severe drought
drought_months = np.where(drought_index < -1.5)[0]
print(f"Severe drought in months: {drought_months}")
```

!!! tip "Interpreting SPI"
    | SPI value | Category |
    |-----------|----------|
    | > 2.0     | Extremely wet |
    | 1.5–2.0   | Very wet |
    | 1.0–1.5   | Moderately wet |
    | -1.0–1.0  | Near normal |
    | -1.5 to -1.0 | Moderately dry |
    | -2.0 to -1.5 | Severely dry |
    | < -2.0    | Extremely dry |

---

## IDW Spatial Interpolation

```python
import numpy as np
from idw import idw

# Four weather stations in a 100×100 grid
station_rows = np.array([10.0, 10.0, 90.0, 90.0])
station_cols = np.array([10.0, 90.0, 10.0, 90.0])
station_values = np.array([120.0, 80.0, 150.0, 60.0])  # mm precipitation

grid = np.zeros((100, 100))
idw(station_rows, station_cols, station_values, grid, power=2)

print(f"Grid range: {grid.min():.1f} – {grid.max():.1f} mm")
```

!!! note "Power parameter"
    - `power=1` — linear distance decay (smooth gradients)
    - `power=2` — squared decay (most common, balanced)
    - `power=3+` — very local influence (sharp transitions near stations)

---

## FAO-56 Evapotranspiration (pyeto)

The `pyeto` sub-package implements the complete FAO-56 Penman-Monteith procedure.

```python
from pyeto import fao, convert

latitude_rad = convert.deg2rad(-15)  # Convert degrees to radians

# Solar parameters for day 15 of the year
day_of_year = 15
sol_dec = fao.sol_dec(day_of_year)
sha = fao.sunset_hour_angle(latitude_rad, sol_dec)
ird = fao.inv_rel_dist_earth_sun(day_of_year)
et_rad = fao.et_rad(latitude_rad, sol_dec, sha, ird)

# Saturation vapour pressure
svp = fao.mean_svp(tmin=18.0, tmax=28.0)

# Actual vapour pressure (from dew point)
avp = fao.avp_from_tdew(tdew=16.0)

# Penman-Monteith ETP
et0 = fao.fao56_penman_monteith(
    net_rad=12.5,  # MJ m-2 day-1
    t=23.0,        # °C
    ws=2.0,        # m/s at 2 m height
    svp=svp,
    avp=avp,
    delta_svp=fao.delta_svp(23.0),
    psy=fao.psy_const(fao.atm_pressure(altitude=500)),
)
print(f"FAO-56 ET₀: {et0:.2f} mm/day")
```

---

## Reading NetCDF Files

```python
from read_nc import read_nc

# Read monthly precipitation from a NetCDF4 file
data, lats, lons, times = read_nc(
    filepath="precip_2000_2020.nc",
    variable="precip",
)
print(data.shape)   # (time, lat, lon)
```
