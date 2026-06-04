# Unit Conversions

`UnitConverter` provides static conversion methods for all common units used in hydroclimate analysis: temperature, angles, wind speed, precipitation depth, and radiation.

## Temperature

```python
from pyutils.utils import UnitConverter

# Celsius ↔ Kelvin
UnitConverter.celsius_to_kelvin(25.0)    # 298.15
UnitConverter.kelvin_to_celsius(298.15)  # 25.0

# Celsius ↔ Fahrenheit
UnitConverter.celsius_to_fahrenheit(100.0)  # 212.0
UnitConverter.fahrenheit_to_celsius(32.0)   # 0.0
```

## Angles (latitude / solar)

```python
import math

UnitConverter.degrees_to_radians(-15.0)   # -0.2618  (latitude -15° South)
UnitConverter.radians_to_degrees(math.pi) # 180.0
```

## Wind speed height correction

Wind measured at 10 m (standard AWS height) needs to be corrected to 2 m (FAO-56 standard) before using it in Penman-Monteith:

```python
wind_2m = UnitConverter.wind_speed_height_correction(
    wind_speed=5.0,          # measured at 10 m
    measurement_height=10.0,
    target_height=2.0,
)
print(f"Wind at 2 m: {wind_2m:.2f} m/s")  # 3.87 m/s
```

## Precipitation depth

```python
UnitConverter.mm_to_meters(1000.0)  # 1.0
UnitConverter.meters_to_mm(0.025)   # 25.0
```

## Solar radiation

```python
# MJ m⁻² day⁻¹ ↔ W m⁻²
UnitConverter.mj_per_m2_to_watts_per_m2(15.0)    # 173.6 W/m²
UnitConverter.watts_per_m2_to_mj_per_m2(173.6)   # ≈ 15.0 MJ/m²/day
```

::: pyutils.utils.conversions.UnitConverter
    options:
      show_signature_annotations: true
