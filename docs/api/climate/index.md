# pyutils.climate

Climate model post-processing — bias correction methods for precipitation and temperature forecasts.

---

## Classes

| Class | Description |
|-------|-------------|
| [`BiasCorrection`](bias_correction.md#biascorrection) | Linear scaling, variance scaling, gamma quantile mapping, empirical quantile mapping, regression-based error correction |

---

## Quick Start

### Linear Scaling

```python
import numpy as np
from pyutils.climate import BiasCorrection

hindcast    = np.random.gamma(2, 10, 120)   # 10-year model history
observation = np.random.gamma(2, 14, 120)   # 10-year observations
forecast    = np.random.gamma(2, 10, 12)    # 1-year model forecast

# Additive correction (good for temperature)
corrected_t = BiasCorrection.linear_scaling(
    forecast, hindcast, observation, method="additive"
)

# Multiplicative correction (good for precipitation)
corrected_p = BiasCorrection.linear_scaling(
    forecast, hindcast, observation, method="multiplicative"
)
```

### Quantile Mapping (Gamma Distribution)

```python
# Recommended for precipitation
corrected = BiasCorrection.quantile_mapping_gamma(forecast, hindcast, observation)
```

### Empirical Quantile Mapping

```python
# Nearest-rank lookup — no distributional assumption
corrected_value = BiasCorrection.quantile_mapping_empirical(
    forecast=25.0,
    hindcast=hindcast,
    observation=observation,
)
```

### Bias-Correct Entire xarray Dataset

```python
from pyutils.climate import BiasCorrection

corrected_ds = BiasCorrection.correct_dataset(
    forecast_ds, hindcast_ds, observation_ds,
    method="variance_scaling",
    var_name="pr",
)
```

---

## API Reference

- [Bias Correction](bias_correction.md) — `BiasCorrection`
