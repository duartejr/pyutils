# pyutils

A Python toolkit for hydroclimate data processing: evapotranspiration, geospatial interpolation, NetCDF I/O, and cartographic visualization.

[![CI](https://github.com/duartejr/pyutils/actions/workflows/ci.yml/badge.svg)](https://github.com/duartejr/pyutils/actions/workflows/ci.yml)

---

## Features

| Module | What it does |
|--------|-------------|
| `evapotranspiration` | Hargreaves, Thornthwaite ETP calculations |
| `pyeto` | FAO-56 Penman-Monteith (full implementation) |
| `read_nc` | Read climate variables from NetCDF4 files |
| `thiessen` | Thiessen polygon rainfall interpolation |
| `idw` | Inverse Distance Weighting interpolation |
| `plot_maps` | Cartopy-based regional map plotting |
| `spi` | Standardized Precipitation Index |
| `rv` | Rainfall distribution correction (Gamma, empirical, linear) |
| `bh_thorthwaite` | Thornthwaite-Mather water balance |
| `corr_med` | Bias correction utilities |
| `vazoes_minimas` | Low-flow statistics (Q7,10, permanence curves) |

---

## Installation

```bash
git clone https://github.com/duartejr/pyutils.git
cd pyutils
pip install -e .
```

**System requirements** (for geospatial modules):

```bash
# Ubuntu / Debian
sudo apt-get install libgeos-dev libproj-dev proj-data proj-bin
```

### Dependencies

Core: `numpy`, `pandas`, `scipy`, `geopandas`, `shapely`, `netCDF4`, `matplotlib`, `scikit-learn`

Optional (geospatial): `cartopy`, `geovoronoi`, `mapclassify`

---

## Quick Start

### Evapotranspiration

```python
import numpy as np
from evapotranspiration import hargreaves, thornthwaite

# Hargreaves ETP (mm/day)
t_med = np.array([25.0, 26.0, 24.5])
t_max = np.array([31.0, 32.0, 30.5])
t_min = np.array([19.0, 20.0, 18.5])

etp = hargreaves(t_med, t_max, t_min, y=-5.0, months=1)
print(etp)  # array of daily ETP values

# Thornthwaite ETP (mm/month)
etp_monthly = thornthwaite(t_med_monthly, latitude=-5.0)
```

### FAO-56 Penman-Monteith (via pyeto)

```python
from pyeto import fao

# Solar radiation from sunshine hours
et_rad = fao.et_rad(latitude_rad, solar_dec, sunset_hour_angle, inv_rel_dist_earth_sun)
sol_rad = fao.sol_rad_from_sun_hours(sunshine_hours, sol_rad_clear_sky, coastal=False)

# Reference ETP
eto = fao.fao56_penman_monteith(
    net_rad=sol_rad,
    t=25.0,
    ws=2.0,
    svp=3.17,
    avp=2.00,
    delta_svp=0.189,
    psy=0.0665,
)
```

### Read NetCDF Climate Data

```python
from read_nc import read_nc, read_nc_shp

# Read entire variable
lats, lons, dates, data = read_nc("era5_precip.nc", "tp")

# Clip to a shapefile extent
lats, lons, dates, data = read_nc_shp("era5_precip.nc", "tp", "basin.shp")
```

### Thiessen Polygon Interpolation

```python
from datetime import date
from thiessen import thiessen

pr_files = ["station_01.csv", "station_02.csv", "station_03.csv"]

# Daily average (Jan 2020)
pr_daily = thiessen(
    shp="basin.shp",
    pr_files=pr_files,
    start_date=date(2020, 1, 1),
    end_date=date(2020, 1, 31),
    kind="daily",
)
print(pr_daily.head())
#         data  pr_media
# 0 2020-01-01      12.3
# 1 2020-01-02       0.0
```

### Standardized Precipitation Index

```python
from spi import spi

# SPI-3 for a monthly precipitation series
spi_values = spi(precip_monthly, scale=3)
```

### Water Balance (Thornthwaite-Mather)

```python
from bh_thorthwaite import bhthorthwaite

deficit, excess = bh_thorthwaite2.bhthorthwaite(
    pr=[80, 60, 40, 20],
    etp=[70, 75, 80, 85],
    cad=100,
)
```

---

## Module Reference

### `evapotranspiration`

| Function | Description |
|----------|-------------|
| `hargreaves(t_med, t_max, t_min, y, months)` | Hargreaves-Samani ETP |
| `thornthwaite(t_med, latitude)` | Thornthwaite monthly ETP |
| `penman_monteith(...)` | Simplified Penman-Monteith |

### `read_nc`

| Function | Description |
|----------|-------------|
| `read_nc(file, var)` | Read NetCDF variable; returns lats, lons, dates, data |
| `read_nc_shp(file, var, shp)` | Same, clipped to shapefile bounding box |

### `pyeto.fao`

Full FAO-56 implementation. See [FAO Irrigation Paper No. 56](http://www.fao.org/3/X0490E/x0490e00.htm).

Key functions: `fao56_penman_monteith`, `hargreaves`, `et_rad`, `sol_rad_from_sun_hours`, `net_rad`, `svp_from_t`, `avp_from_rhmin_rhmax`.

### `rv` subpackage

| Module | Description |
|--------|-------------|
| `rv.RvGama` | Gamma distribution quantile mapping |
| `rv.RvEmpirica` | Empirical distribution correction |
| `rv.RvLinear` | Linear regression correction |
| `rv.RvMedia` | Mean-based correction |
| `rv.RvOns` | ONS (Brazilian grid operator) correction |

---

## Project Structure

```
pyutils/
├── pyutils.py            # Package metadata (__version__, __author__)
├── evapotranspiration.py # Hargreaves, Thornthwaite, Penman-Monteith
├── read_nc.py            # NetCDF4 I/O utilities
├── thiessen.py           # Thiessen polygon interpolation
├── idw.py                # Inverse distance weighting
├── plot_maps.py          # Cartopy map plotting
├── spi.py                # Standardized Precipitation Index
├── bh_thorthwaite.py     # Water balance (simple)
├── bh_thorthwaite2.py    # Water balance (Thornthwaite-Mather)
├── ccsthorthwaite.py     # Thornthwaite for climate change scenarios
├── corr_med.py           # Bias correction
├── vazoes_minimas.py     # Low-flow statistics
├── buffer.py             # Shapefile buffer operations
├── isinpoly.py           # Point-in-polygon test
├── isinshp.py            # Point-in-shapefile test
├── shp_area.py           # Shapefile area calculation
├── getvertshp.py         # Extract shapefile vertices
├── custom_colorbar.py    # Custom matplotlib colorbars
├── rvGamma.py            # Gamma correction for ensemble forecasts
├── eto.py                # Simplified ETP calculations
├── pyeto/                # FAO-56 Penman-Monteith (full)
│   ├── fao.py
│   ├── convert.py
│   └── thornthwaite.py
└── rv/                   # Rainfall distribution corrections
    ├── RvGama.py
    ├── RvEmpirica.py
    ├── RvLinear.py
    ├── RvMedia.py
    └── RvOns.py
```

---

## Development

```bash
# Clone and install in editable mode with dev dependencies
git clone https://github.com/duartejr/pyutils.git
cd pyutils
pip install -e ".[dev]"

# Run verification tests
python test_phase1_verification.py
python test_phase2_verification.py
python test_phase3_verification.py
```

---

## License

GNU General Public License v3.0 — see [LICENSE](LICENSE) for details.

## Author

Duarte Junior — duarte.jr105@gmail.com
