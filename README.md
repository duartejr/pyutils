# pyutils

A Python toolkit for hydroclimate data processing: evapotranspiration, water balance, spatial interpolation, drought indices, NetCDF I/O with xarray, and cartographic visualization.

[![CI](https://github.com/duartejr/pyutils/actions/workflows/ci.yml/badge.svg)](https://github.com/duartejr/pyutils/actions/workflows/ci.yml)
[![Python](https://img.shields.io/badge/python-3.11%2B-blue)](https://www.python.org/)
[![License](https://img.shields.io/badge/license-GPL--3.0--or--later-green)](LICENSE)
[![Docs](https://img.shields.io/badge/docs-online-blue)](https://duartejr.github.io/pyutils)

📖 **[Full documentation → duartejr.github.io/pyutils](https://duartejr.github.io/pyutils)**

---

## Overview

pyutils is organised into domain-specific sub-packages. Each one consolidates several related modules behind a clean, class-based API.

| Sub-package | Classes | What it covers |
|-------------|---------|----------------|
| `pyutils.hydrology` | `Hargreaves`, `Thornthwaite`, `PenmanMonteith`, `ThornthwaiteMather`, `StandardizedPrecipitationIndex`, `FlowAnalyzer`, `TimeOfConcentration` | ET models, water balance, SPI drought index, stream flow statistics, time of concentration |
| `pyutils.geospatial` | `XarrayNetCDFHandler`, `InverseDistanceWeighting`, `ThiessenPolygon`, `ShapefileHandler`, `MapRenderer` | NetCDF I/O via xarray, IDW/Thiessen interpolation, shapefile GIS ops, thematic maps |
| `pyutils.climate` | `BiasCorrection` | Linear scaling, variance scaling, quantile mapping for model post-processing |
| `pyutils.utils` | `DataValidator`, `UnitConverter` | Input validation, unit conversions |
| `pyutils.core` | `BaseModel`, `ValidationError`, constants | Shared base classes and physical constants |

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

**Dependencies:** `numpy`, `pandas`, `scipy`, `geopandas`, `shapely`, `netCDF4`, `xarray`, `dask`, `matplotlib`, `scikit-learn`

Optional: `cartopy` (advanced map projections), `geovoronoi`, `mapclassify`

---

## Quick Start

### Evapotranspiration

```python
import numpy as np
from pyutils.hydrology import Hargreaves, Thornthwaite

t_med = np.array([22, 23, 24, 25, 26, 27, 27, 26, 25, 24, 23, 22], dtype=float)
t_max = t_med + 5
t_min = t_med - 5

# Hargreaves-Samani (mm/month) — only needs T_min, T_max, T_mean
et = Hargreaves().compute(t_med=t_med, t_max=t_max, t_min=t_min, y=-5.0, months=1)
print(et.round(1))
# [11.7 11.9 12.1 12.4 12.6 12.8 12.8 12.6 12.4 12.1 11.9 11.7]

# Thornthwaite (mm/month)
pet = Thornthwaite().compute(t_med=t_med, y=-5.0)
```

### Water Balance (Thornthwaite-Mather)

```python
from pyutils.hydrology import ThornthwaiteMather

precipitation = np.array([200, 150, 100,  50,  30,  20,
                            15,  20,  40,  80, 150, 200], dtype=float)
potential_et  = np.array([ 40,  45,  60,  80, 100, 110,
                           120, 110,  90,  70,  50,  40], dtype=float)

model = ThornthwaiteMather(cap=100.0)
result = model.compute(precipitation=precipitation, potential_et=potential_et)

print("Annual runoff  :", result["runoff"].sum(), "mm")
print("Annual deficit :", result["water_deficit"].sum(), "mm")
```

### Drought Index (SPI)

```python
import numpy as np
from pyutils.hydrology import StandardizedPrecipitationIndex

monthly_precip = np.random.gamma(shape=2, scale=60, size=360)   # 30 years

spi3 = StandardizedPrecipitationIndex(timescale=3)
spi_values, _, _ = spi3.compute(monthly_precip)

# SPI < -1.5 → severe drought
drought_months = np.sum(spi_values < -1.5)
print(f"Severe drought months: {drought_months}")
```

### Stream Flow Analysis

```python
from pyutils.hydrology import FlowAnalyzer

daily_flow = ...  # np.ndarray of daily discharge (m³/s)

stats = FlowAnalyzer.low_flow_statistics(daily_flow)
print(f"Q90: {stats['q90']:.3f} m³/s")
print(f"Q95: {stats['q95']:.3f} m³/s")

q7_10 = FlowAnalyzer.q7_10_extreme_value(daily_flow)
print(f"Q7,10 = {q7_10:.3f} m³/s")
```

### Time of Concentration

```python
from pyutils.hydrology import TimeOfConcentration

results = TimeOfConcentration.compute_all(
    area_km2=25.0, length_km=8.0, height_m=120.0, slope_pct=1.5,
)
for method, tc in sorted(results.items(), key=lambda kv: kv[1]):
    print(f"{method:<18} {tc:6.2f} h")
```

### NetCDF Climate Data (xarray)

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

### Bias Correction

```python
import numpy as np
from pyutils.climate import BiasCorrection

hindcast    = np.random.gamma(2, 10, 120)   # 10-year model history
observation = np.random.gamma(2, 14, 120)   # 10-year observations
forecast    = np.random.gamma(2, 10, 12)    # 1-year model forecast

# Quantile mapping (recommended for precipitation)
corrected = BiasCorrection.quantile_mapping_gamma(forecast, hindcast, observation)

# Linear scaling (additive — good for temperature)
corrected_t = BiasCorrection.linear_scaling(forecast, hindcast, observation, method="additive")
```

### Unit Conversions

```python
from pyutils.utils import UnitConverter

UnitConverter.celsius_to_kelvin(25.0)            # 298.15
UnitConverter.degrees_to_radians(-15.0)          # -0.2618

# Wind at 10 m → 2 m (FAO-56 standard)
wind_2m = UnitConverter.wind_speed_height_correction(5.0, measurement_height=10.0)
```

---

## Project Structure

```
pyutils/
│
├── pyutils/                   # Class-based sub-packages
│   ├── core/                  # Abstract base classes, exceptions, constants
│   ├── hydrology/             # ET, water balance, SPI, flow analysis
│   ├── geospatial/            # xarray NetCDF, IDW, Thiessen, shapefiles, maps
│   ├── climate/               # Bias correction (linear, variance, quantile)
│   └── utils/                 # Validators, unit converters
│
├── pyeto/                     # FAO-56 Penman-Monteith engine (internal)
└── tests/                     # Test suite
```

---

## Development

```bash
git clone https://github.com/duartejr/pyutils.git
cd pyutils
pip install -e ".[dev]"

# Run tests
pytest tests/ -v

# Build docs locally
pip install -e ".[docs]"
mkdocs serve
# → http://127.0.0.1:8000
```

---

## License

GNU General Public License v3.0 — see [LICENSE](LICENSE) for details.

---

## Contact

**Duarte Junior**

- Email: [duarte.jr105@gmail.com](mailto:duarte.jr105@gmail.com)
- LinkedIn: [linkedin.com/in/duarte-junior-9530595a](https://www.linkedin.com/in/duarte-junior-9530595a)
- GitHub: [github.com/duartejr](https://github.com/duartejr)
- Documentation: [duartejr.github.io/pyutils](https://duartejr.github.io/pyutils)
