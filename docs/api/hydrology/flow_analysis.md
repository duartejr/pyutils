# Flow Analysis

`FlowAnalyzer` computes standard hydrological indices for stream flow series: flow duration curves, minimum flow quantiles, Q7,10, and comprehensive low/high flow statistics.

## Flow Duration Curve

```python
import numpy as np
import matplotlib.pyplot as plt
from pyutils.hydrology import FlowAnalyzer

np.random.seed(7)
# Simulate 10 years of daily flow (m³/s)
daily_flow = np.abs(np.random.gamma(shape=2, scale=5, size=365 * 10))

exceedance, flow_sorted = FlowAnalyzer.flow_duration_curve(daily_flow)

plt.figure(figsize=(8, 4))
plt.plot(exceedance, flow_sorted, color="steelblue")
plt.yscale("log")
plt.xlabel("Exceedance probability (%)")
plt.ylabel("Flow (m³/s)")
plt.title("Flow Duration Curve")
plt.grid(True, alpha=0.3)
plt.tight_layout()
plt.show()
```

## Low-flow Indices

```python
from pyutils.hydrology import FlowAnalyzer

stats = FlowAnalyzer.low_flow_statistics(daily_flow)

print(f"Q50  (median flow) : {stats['q50']:.2f} m³/s")
print(f"Q90  (low flow)    : {stats['q90']:.2f} m³/s")
print(f"Q95               : {stats['q95']:.2f} m³/s")
print(f"Q99               : {stats['q99']:.2f} m³/s")
print(f"Mean annual flow   : {stats['mean_annual']:.2f} m³/s")
print(f"CV                 : {stats['cv']:.2f}")
```

## Q7,10 — 7-day, 10-year minimum flow

The Q7,10 is the minimum 7-day average flow with a 10-year return period. It is used to define minimum environmental flows and water availability thresholds.

```python
q7_10 = FlowAnalyzer.q7_10_extreme_value(daily_flow, return_period=10.0)
print(f"Q7,10 = {q7_10:.3f} m³/s")

# Change return period
q7_20 = FlowAnalyzer.q7_10_extreme_value(daily_flow, return_period=20.0)
print(f"Q7,20 = {q7_20:.3f} m³/s")  # More conservative (lower)
```

## High-flow Statistics

```python
high = FlowAnalyzer.high_flow_statistics(daily_flow)

print(f"Q1  (1% exceedance)  : {high['q1']:.2f} m³/s")
print(f"Q5  (5% exceedance)  : {high['q5']:.2f} m³/s")
print(f"Q10 (10% exceedance) : {high['q10']:.2f} m³/s")
print(f"Maximum flow         : {high['max_flow']:.2f} m³/s")
```

::: pyutils.hydrology.flow_analysis.FlowAnalyzer
    options:
      show_signature_annotations: true
