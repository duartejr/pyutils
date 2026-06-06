# pyutils.utils

Validation and conversion utilities — input validation and physical unit conversions for hydroclimate data.

---

## Classes

| Class | Description |
|-------|-------------|
| [`DataValidator`](validators.md#datavalidator) | Validate precipitation, temperature, humidity, wind speed, spatial data, and time series length |
| [`UnitConverter`](conversions.md#unitconverter) | Temperature, angle, wind speed, precipitation depth, and solar radiation conversions |

---

## Quick Start

### Input Validation

```python
import numpy as np
from pyutils.utils import DataValidator
from pyutils.core import ValidationError

precip = np.array([10.0, 25.0, 5.0, 80.0, 0.0])

try:
    DataValidator.validate_precipitation(precip, allow_zero=False)
except ValidationError as e:
    print(f"Invalid data: {e}")

# Temperature validation
temp = np.array([22.0, 23.5, 21.0, 24.0])
DataValidator.validate_temperature(temp)  # raises if outside -90 to 60 °C

# Relative humidity validation
rh = np.array([65.0, 70.0, 80.0])
DataValidator.validate_relative_humidity(rh)  # raises if outside 0–100 %
```

### Unit Conversions

```python
from pyutils.utils import UnitConverter

# Temperature
UnitConverter.celsius_to_kelvin(25.0)          # 298.15 K
UnitConverter.celsius_to_fahrenheit(100.0)     # 212.0 °F

# Angles (latitude / solar calculations)
UnitConverter.degrees_to_radians(-15.0)        # -0.2618 rad

# Wind speed — measured at 10 m → FAO-56 standard 2 m
wind_2m = UnitConverter.wind_speed_height_correction(
    wind_speed=5.0, measurement_height=10.0
)

# Precipitation depth
UnitConverter.mm_to_meters(1000.0)             # 1.0 m

# Solar radiation
UnitConverter.mj_per_m2_to_watts_per_m2(15.0) # 173.6 W/m²
```

---

## API Reference

- [Validators](validators.md) — `DataValidator`
- [Unit Conversions](conversions.md) — `UnitConverter`
