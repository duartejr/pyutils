# pyutils

**Hydroclimate analysis toolkit for Python** — evapotranspiration, spatial interpolation, precipitation indices, NetCDF I/O, and cartographic utilities.

[![CI](https://github.com/duartejr/pyutils/actions/workflows/ci.yml/badge.svg)](https://github.com/duartejr/pyutils/actions/workflows/ci.yml)
[![Python](https://img.shields.io/badge/python-3.11%2B-blue)](https://www.python.org/)
[![License](https://img.shields.io/badge/license-GPL--3.0--or--later-green)](https://github.com/duartejr/pyutils/blob/master/LICENSE)

---

## Overview

pyutils is organised into domain-specific sub-packages. Each one consolidates related functionality behind a clean, class-based API.

| Sub-package | Classes | What it covers |
|-------------|---------|----------------|
| [`pyutils.hydrology`](api/hydrology/index.md) | `Hargreaves`, `Thornthwaite`, `PenmanMonteith`, `ThornthwaiteMather`, `StandardizedPrecipitationIndex`, `FlowAnalyzer` | ET models, water balance, SPI drought index, stream flow statistics |
| [`pyutils.geospatial`](api/geospatial/index.md) | `XarrayNetCDFHandler`, `InverseDistanceWeighting`, `ThiessenPolygon`, `ShapefileHandler`, `MapRenderer` | NetCDF I/O via xarray, IDW/Thiessen interpolation, shapefile GIS ops, thematic maps |
| [`pyutils.climate`](api/climate/index.md) | `BiasCorrection` | Linear scaling, variance scaling, quantile mapping for model post-processing |
| [`pyutils.utils`](api/utils/index.md) | `DataValidator`, `UnitConverter` | Input validation, unit conversions |

---

## Quick Install

```bash
git clone https://github.com/duartejr/pyutils.git
cd pyutils
pip install -e .
```

Requires **Python 3.11+**. See [Installation](getting-started/installation.md) for full system requirements.

---

## Highlights

### Evapotranspiration

```python
import numpy as np
from pyutils.hydrology import Hargreaves, PenmanMonteith

t_med = np.array([22, 23, 24, 25, 26, 27, 27, 26, 25, 24, 23, 22], dtype=float)
t_max = t_med + 5
t_min = t_med - 5

# Hargreaves-Samani (mm/month)
et = Hargreaves().compute(t_med=t_med, t_max=t_max, t_min=t_min, y=-5.0, months=1)

# FAO-56 Penman-Monteith (mm/day)
et0 = PenmanMonteith().compute(
    t_mean=25.0, t_max=32.0, t_min=18.0,
    rh_mean=60.0, u2=2.0,
    latitude=-15.0, day_of_year=180,
    elevation=500.0,
)
```

### Drought Monitoring with SPI

```python
import numpy as np
from pyutils.hydrology import StandardizedPrecipitationIndex

monthly_precip = np.random.gamma(shape=2, scale=60, size=360)  # 30 years
spi_values, _, _ = StandardizedPrecipitationIndex(timescale=3).compute(monthly_precip)
```

### Spatial Interpolation

```python
import numpy as np
from pyutils.geospatial import InverseDistanceWeighting, ThiessenPolygon

x = np.array([-45.0, -44.0, -43.0, -44.0])
y = np.array([ -5.0,  -5.0,  -5.0,  -6.0])
z = np.array([  85.0, 120.0,  95.0, 110.0])

xi = np.linspace(-46.0, -42.0, 40)
yi = np.linspace( -7.0,  -4.0, 30)

grid = InverseDistanceWeighting(power=2).interpolate(x, y, z, xi, yi)
```

### NetCDF Climate Data

```python
from pyutils.geospatial import XarrayNetCDFHandler

handler = XarrayNetCDFHandler()
ds = handler.read("era5_precip_1991_2020.nc")
climatology = handler.compute_climatology(group_by="month")
anomaly = handler.compute_anomaly()
```

---

## Changelog

See the [Changelog](changelog.md) for a history of releases.
