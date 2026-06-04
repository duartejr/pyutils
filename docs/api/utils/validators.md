# Data Validation

`DataValidator` provides centralized input validation for all pyutils modules. Every method raises `ValidationError` on failure, making it easy to catch and report problems early.

## Precipitation

```python
import numpy as np
from pyutils.utils import DataValidator
from pyutils.core import ValidationError

# Valid series
DataValidator.validate_precipitation(np.array([0.0, 10.5, 200.0]))  # OK

# Catch invalid data
try:
    DataValidator.validate_precipitation(np.array([-5.0, 10.0]))
except ValidationError as e:
    print(f"Error: {e}")  # Error: Negative precipitation values found

# Allow NaN (missing data)
DataValidator.validate_precipitation(
    np.array([np.nan, 10.0, 50.0]),
    allow_missing=True,
)
```

## Temperature

```python
# Valid range: −100 °C to +60 °C
DataValidator.validate_temperature(np.array([-10.0, 0.0, 35.0]))  # OK

try:
    DataValidator.validate_temperature(np.array([70.0]))
except ValidationError as e:
    print(f"Error: {e}")  # Error: Temperature exceeds maximum (60.0°C)
```

## Other variables

```python
# Relative humidity (0–100%)
DataValidator.validate_relative_humidity(np.array([20.0, 60.0, 95.0]))

# Wind speed (0–50 m/s)
DataValidator.validate_wind_speed(np.array([0.0, 3.5, 12.0]))

# Spatial coordinates: x, y, z must all have the same length and ≥ 2 points
DataValidator.validate_spatial_data(
    x=np.array([0.0, 10.0, 20.0]),
    y=np.array([0.0, 10.0, 20.0]),
    z=np.array([100.0, 150.0, 120.0]),
)

# Time series minimum length check
DataValidator.validate_time_series_length(np.ones(60), min_length=30)
```

::: pyutils.utils.validators.DataValidator
    options:
      show_signature_annotations: true
