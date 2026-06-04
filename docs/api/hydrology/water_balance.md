# Water Balance

`ThornthwaiteMather` computes the classical monthly water balance — soil moisture, actual ET, runoff, water deficit and surplus. It consolidates three previous separate implementations (`bh_thorthwaite.py`, `bh_thorthwaite2.py`, `ccsthorthwaite.py`) into a single class with a `method` parameter.

## Quick start

```python
import numpy as np
from pyutils.hydrology import ThornthwaiteMather

# Monthly precipitation and potential ET (mm)
precipitation = np.array([200, 150, 100,  50,  30,  20,
                            15,  20,  40,  80, 150, 200], dtype=float)
potential_et  = np.array([ 40,  45,  60,  80, 100, 110,
                           120, 110,  90,  70,  50,  40], dtype=float)

model = ThornthwaiteMather(cap=100.0, method="standard")
result = model.compute(precipitation=precipitation, potential_et=potential_et)

print("Annual runoff  :", result["runoff"].sum().round(1), "mm")
print("Annual AET     :", result["actual_et"].sum().round(1), "mm")
print("Annual deficit :", result["water_deficit"].sum().round(1), "mm")
print("Annual surplus :", result["water_surplus"].sum().round(1), "mm")
```

## Choosing a method

| `method` | Description |
|----------|-------------|
| `"standard"` | Classic Thornthwaite-Mather (default) |
| `"variant2"` | Slightly relaxed runoff threshold (+10% capacity) |
| `"ccs"` | Reduced ET under water stress (climate change scenarios) |

## Visualising the water balance

```python
import matplotlib.pyplot as plt

months = ["Jan","Feb","Mar","Apr","May","Jun",
          "Jul","Aug","Sep","Oct","Nov","Dec"]

fig, ax = plt.subplots(figsize=(10, 4))
ax.bar(months, result["water_surplus"], label="Surplus", color="steelblue")
ax.bar(months, result["water_deficit"], bottom=-result["water_deficit"],
       label="Deficit", color="tomato")
ax.plot(months, result["soil_moisture"], "k--o", label="Soil moisture")
ax.set_ylabel("mm")
ax.legend()
plt.tight_layout()
plt.show()
```

::: pyutils.hydrology.water_balance.ThornthwaiteMather
    options:
      show_signature_annotations: true
      members: [compute]
