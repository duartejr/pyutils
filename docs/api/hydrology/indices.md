# Hydrological Indices

## StandardizedPrecipitationIndex

Quantifies precipitation deficits at multiple timescales. Positive values indicate wet conditions; negative values indicate drought. The McKee (1993) classification is:

| SPI range | Category |
|-----------|----------|
| ≥ 2.0 | Extremely wet |
| 1.5 – 1.99 | Very wet |
| 1.0 – 1.49 | Moderately wet |
| −0.99 – 0.99 | Near normal |
| −1.0 – −1.49 | Moderately dry |
| −1.5 – −1.99 | Severely dry |
| ≤ −2.0 | Extremely dry |

## Example

```python
import numpy as np
import matplotlib.pyplot as plt
from pyutils.hydrology import StandardizedPrecipitationIndex

np.random.seed(42)
# Simulate 30 years of monthly precipitation (mm)
monthly_precip = np.random.gamma(shape=2, scale=60, size=360)

spi3 = StandardizedPrecipitationIndex(timescale=3)
spi_values, accumulated, probabilities = spi3.compute(monthly_precip)

# Plot SPI-3
fig, ax = plt.subplots(figsize=(12, 3))
x = np.arange(len(spi_values))
valid = ~np.isnan(spi_values)
ax.bar(x[valid & (spi_values >= 0)], spi_values[valid & (spi_values >= 0)],
       color="steelblue", label="Wet")
ax.bar(x[valid & (spi_values < 0)], spi_values[valid & (spi_values < 0)],
       color="tomato", label="Dry")
ax.axhline(-1.5, color="darkred", linestyle="--", linewidth=0.8, label="Severe drought")
ax.set_xlabel("Month")
ax.set_ylabel("SPI-3")
ax.legend()
plt.tight_layout()
plt.show()
```

## Multi-timescale comparison

```python
results = {}
for scale in [1, 3, 6, 12]:
    spi_obj = StandardizedPrecipitationIndex(timescale=scale)
    spi_vals, _, _ = spi_obj.compute(monthly_precip)
    results[f"SPI-{scale}"] = spi_vals

# Shorter timescales respond faster; longer ones capture multi-year droughts
```

::: pyutils.hydrology.indices.StandardizedPrecipitationIndex
    options:
      show_signature_annotations: true
      members: [compute]
