# Evapotranspiration

Three evapotranspiration models covering the full spectrum from data-sparse (Hargreaves) to data-rich (Penman-Monteith) conditions. All classes live in `pyutils.hydrology`.

## Quick comparison

| Class | Inputs required | Accuracy | Use case |
|-------|----------------|----------|----------|
| `Hargreaves` | T_mean, T_max, T_min, latitude | Moderate | Remote areas, limited data |
| `Thornthwaite` | T_mean, latitude | Low–moderate | Regional climate classification |
| `PenmanMonteith` | Full met. data | High | Precision irrigation, research |

---

## Hargreaves

Temperature-based model requiring only minimum, maximum, and mean air temperature.

### Example

```python
import numpy as np
from pyutils.hydrology import Hargreaves

model = Hargreaves()

t_med = np.array([22, 23, 24, 25, 26, 27, 27, 26, 25, 24, 23, 22], dtype=float)
t_max = t_med + 5
t_min = t_med - 5

# ET for a station at latitude -5° (near Equator, Northern Brazil)
et = model.compute(t_med=t_med, t_max=t_max, t_min=t_min, y=-5.0, months=1)

print(et.round(1))
# [11.7 11.9 12.1 12.4 12.6 12.8 12.8 12.6 12.4 12.1 11.9 11.7]  mm/month
```

::: pyutils.hydrology.evapotranspiration.Hargreaves
    options:
      show_signature_annotations: true
      members: [compute]

---

## Thornthwaite

Potential evapotranspiration from mean monthly temperature and latitude only.

### Example

```python
import numpy as np
from pyutils.hydrology import Thornthwaite

model = Thornthwaite()

t_med = np.array([28, 28, 27, 26, 25, 24, 24, 25, 26, 27, 28, 28], dtype=float)

pet = model.compute(t_med=t_med, y=-15.0)
print(pet.round(1))
```

::: pyutils.hydrology.evapotranspiration.Thornthwaite
    options:
      show_signature_annotations: true
      members: [compute]

---

## PenmanMonteith

FAO-56 reference evapotranspiration using the full set of meteorological variables.

### Example

```python
from pyutils.hydrology import PenmanMonteith

model = PenmanMonteith()

et0 = model.compute(
    t_mean=25.0,
    t_max=31.0,
    t_min=19.0,
    rh_mean=70.0,
    u2=2.5,           # wind speed at 2 m (m/s)
    elevation=500.0,  # metres above sea level
)
print(f"ET0 = {et0:.2f} mm/day")
```

::: pyutils.hydrology.evapotranspiration.PenmanMonteith
    options:
      show_signature_annotations: true
      members: [compute]
