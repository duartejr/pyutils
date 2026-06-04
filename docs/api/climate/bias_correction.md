# Bias Correction

`BiasCorrection` corrects systematic errors in climate model output (hindcast or forecast) by mapping its distribution toward observed statistics. All methods are static — no instantiation needed.

## When to use each method

| Method | Best for | Preserves |
|--------|----------|-----------|
| `linear_scaling` | Temperature, wind | Mean |
| `variance_scaling` | Temperature, precipitation | Mean + variance |
| `quantile_mapping_gamma` | Precipitation | Full distribution + zero probability |

---

## Linear Scaling

Shifts the forecast mean to match observations. Two variants:

```python
import numpy as np
from pyutils.climate import BiasCorrection

# 10 years of monthly historical data
hindcast    = np.random.normal(loc=15, scale=3, size=120)
observation = np.random.normal(loc=18, scale=4, size=120)
forecast    = np.random.normal(loc=14, scale=3, size=12)

# Additive — good for temperature
corrected = BiasCorrection.linear_scaling(
    forecast, hindcast, observation, method="additive"
)
print(f"Raw mean   : {forecast.mean():.1f}")
print(f"Corrected  : {corrected.mean():.1f}")
print(f"Obs mean   : {observation.mean():.1f}")

# Multiplicative — preserves ratios (avoid with data near zero)
corrected_mult = BiasCorrection.linear_scaling(
    forecast, hindcast, observation, method="multiplicative"
)
```

## Variance Scaling

Corrects both mean and month-to-month variance. The correction is applied month-by-month, so January hindcast statistics are used only for January forecasts.

```python
corrected = BiasCorrection.variance_scaling(forecast, hindcast, observation)
print(f"Raw std      : {forecast.std():.2f}")
print(f"Corrected std: {corrected.std():.2f}")
print(f"Observed std : {observation.std():.2f}")
```

## Quantile Mapping (Gamma)

Remaps the full precipitation distribution, preserving the probability of zero rain from observations. Suitable for sub-monthly to seasonal precipitation forecasts.

```python
# Precipitation data (non-negative)
hindcast_pr    = np.abs(np.random.gamma(shape=2, scale=10, size=120))
observation_pr = np.abs(np.random.gamma(shape=2, scale=14, size=120))
forecast_pr    = np.abs(np.random.gamma(shape=2, scale=10, size=12))

corrected_pr = BiasCorrection.quantile_mapping_gamma(
    forecast_pr, hindcast_pr, observation_pr
)

print(f"Raw mean   : {forecast_pr.mean():.1f} mm")
print(f"Corrected  : {corrected_pr.mean():.1f} mm")
print(f"Obs mean   : {observation_pr.mean():.1f} mm")
```

## Applying to xarray datasets

```python
import xarray as xr
from pyutils.climate import BiasCorrection

forecast_ds    = xr.open_dataset("gcm_forecast.nc")
hindcast_ds    = xr.open_dataset("gcm_hindcast.nc")
observation_ds = xr.open_dataset("observed.nc")

corrected_ds = BiasCorrection.correct_dataset(
    forecast_ds, hindcast_ds, observation_ds,
    method="variance_scaling",
    var_name="pr",
)
corrected_ds.to_netcdf("corrected_forecast.nc")
```

::: pyutils.climate.bias_correction.BiasCorrection
    options:
      show_signature_annotations: true
      members: [linear_scaling, variance_scaling, quantile_mapping_gamma, correct_dataset]
