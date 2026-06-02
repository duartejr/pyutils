# pyutils

**Hydroclimate analysis toolkit for Python** — evapotranspiration, spatial interpolation, precipitation indices, NetCDF I/O, and cartographic utilities.

[![CI](https://github.com/duartejr/pyutils/actions/workflows/ci.yml/badge.svg)](https://github.com/duartejr/pyutils/actions/workflows/ci.yml)
[![Python](https://img.shields.io/badge/python-3.11%2B-blue)](https://www.python.org/)
[![License](https://img.shields.io/badge/license-GPL--3.0--or--later-green)](https://github.com/duartejr/pyutils/blob/master/LICENSE)

---

## Overview

pyutils provides a collection of scientific Python utilities for hydroclimate data processing. All modules are designed for use in research workflows involving atmospheric science, hydrology, and geospatial analysis.

| Module | Purpose |
|--------|---------|
| [`evapotranspiration`](api/evapotranspiration.md) | Hargreaves, Thornthwaite, and Penman-Monteith ETP |
| [`spi`](api/spi.md) | Standardized Precipitation Index |
| [`idw`](api/idw.md) | Inverse Distance Weighting spatial interpolation |
| [`pyeto`](api/pyeto/fao.md) | FAO-56 reference evapotranspiration (full implementation) |
| [`read_nc`](api/read_nc.md) | Read climate variables from NetCDF4 files |
| [`thiessen`](api/thiessen.md) | Thiessen polygon rainfall interpolation |
| [`plot_maps`](api/plot_maps.md) | Cartopy-based regional map plotting |
| [`corr_med`](api/corr_med.md) | Bias correction utilities |

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

### Vectorized IDW Interpolation

Fills a spatial grid from station observations using Inverse Distance Weighting.
The implementation uses NumPy broadcasting — **30–45× faster** than a loop-based approach.

```python
import numpy as np
from idw import idw

station_rows = np.array([0.0, 0.0, 9.0, 9.0])
station_cols = np.array([0.0, 9.0, 0.0, 9.0])
station_values = np.array([10.0, 20.0, 30.0, 40.0])
grid = np.zeros((10, 10))

idw(station_rows, station_cols, station_values, grid, power=2)
```

### Drought Monitoring with SPI

Computes the Standardized Precipitation Index — positive values are wet, negative are drought.

```python
import numpy as np
from spi import spi

precipitation = np.random.exponential(scale=80, size=120)  # 10-year monthly series
drought_index = spi(precipitation, accumulation_scale=3, num_seasons=12)
```

### Multi-Method Evapotranspiration

```python
import numpy as np
from evapotranspiration import hargreaves

mean_temps = np.array([25, 25.5, 26, 25.5, 25, 24.5, 24, 24, 24.5, 25, 25.5, 25.5])
max_temps = mean_temps + 5
min_temps = mean_temps - 5

etp = hargreaves(mean_temps, max_temps, min_temps, y=-15, months=np.arange(1, 13))
```

---

## Changelog

See the [Changelog](changelog.md) for a history of releases.
