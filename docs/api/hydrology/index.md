# pyutils.hydrology

Hydroclimate analysis — evapotranspiration models, water balance, drought indices, and stream flow statistics.

---

## Classes

| Class | Description |
|-------|-------------|
| [`Hargreaves`](evapotranspiration.md#hargreaves) | Temperature-based ET (Hargreaves-Samani) |
| [`Thornthwaite`](evapotranspiration.md#thornthwaite) | Temperature-based potential ET (Thornthwaite 1948) |
| [`PenmanMonteith`](evapotranspiration.md#penmanmonteith) | Full FAO-56 reference ET with radiation and humidity inputs |
| [`ThornthwaiteMather`](water_balance.md#thornthwaitemather) | Monthly water balance (standard, variant 2, and CCS methods) |
| [`StandardizedPrecipitationIndex`](indices.md#standardizedprecipitationindex) | SPI drought index (gamma or lognormal distribution) |
| [`FlowAnalyzer`](flow_analysis.md#flowanalyzer) | Flow duration curves, Q7,10, low/high flow statistics |
| [`TimeOfConcentration`](time_concentration.md#timeofconcentration) | Watershed time of concentration via 13 empirical formulas |

---

## Quick Start

### Evapotranspiration

```python
import numpy as np
from pyutils.hydrology import Hargreaves, Thornthwaite, PenmanMonteith

t_med = np.array([22, 23, 24, 25, 26, 27, 27, 26, 25, 24, 23, 22], dtype=float)
t_max = t_med + 5
t_min = t_med - 5

# Hargreaves-Samani (mm/month)
et_harg = Hargreaves().compute(t_med=t_med, t_max=t_max, t_min=t_min, y=-5.0, months=1)

# Thornthwaite (mm/month)
et_thorn = Thornthwaite().compute(t_med=t_med, y=-5.0)

# FAO-56 Penman-Monteith (mm/day — single day)
et0_pm = PenmanMonteith().compute(
    t_mean=25.0, t_max=32.0, t_min=18.0,
    rh_mean=60.0, u2=2.0,
    latitude=-15.0, day_of_year=180,
    elevation=500.0,
)
```

### Water Balance

```python
from pyutils.hydrology import ThornthwaiteMather

precipitation = np.array([200, 150, 100,  50,  30,  20,
                            15,  20,  40,  80, 150, 200], dtype=float)
potential_et  = np.array([ 40,  45,  60,  80, 100, 110,
                           120, 110,  90,  70,  50,  40], dtype=float)

result = ThornthwaiteMather(cap=100.0).compute(
    precipitation=precipitation, potential_et=potential_et
)
print("Annual runoff:", result["runoff"].sum(), "mm")
print("Annual deficit:", result["water_deficit"].sum(), "mm")
```

### Drought Index (SPI)

```python
from pyutils.hydrology import StandardizedPrecipitationIndex

monthly_precip = np.random.gamma(shape=2, scale=60, size=360)   # 30 years

spi3 = StandardizedPrecipitationIndex(timescale=3)
spi_values, _, _ = spi3.compute(monthly_precip)

drought_months = (spi_values < -1.5).sum()
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

---

## API Reference

- [Evapotranspiration](evapotranspiration.md) — `Hargreaves`, `Thornthwaite`, `PenmanMonteith`
- [Water Balance](water_balance.md) — `ThornthwaiteMather`
- [Drought Indices](indices.md) — `StandardizedPrecipitationIndex`
- [Flow Analysis](flow_analysis.md) — `FlowAnalyzer`
- [Time of Concentration](time_concentration.md) — `TimeOfConcentration`
