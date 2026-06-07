# Time of Concentration

`TimeOfConcentration` estimates how long runoff takes to travel from the most hydraulically distant point of a watershed to its outlet — a key input for the rational method, flood routing, and drainage design. Thirteen widely used empirical formulas are provided as static methods, plus `compute_all` to compare them side by side.

Each formula was calibrated on different watershed types, sizes, and regions, so estimates can vary considerably for the same watershed — comparing several methods is good practice.

## Single Method

```python
from pyutils.hydrology import TimeOfConcentration

# Watershed characteristics
area_km2 = 25.0      # drainage area
length_km = 8.0      # main channel / talweg length
height_m = 120.0     # elevation difference (most distant point -> outlet)
slope_pct = 1.5      # average channel slope

tc = TimeOfConcentration.kirpich(length_km=length_km, slope_pct=slope_pct)
print(f"Kirpich Tc = {tc:.2f} hours")
```

## Comparing All Methods

```python
from pyutils.hydrology import TimeOfConcentration

results = TimeOfConcentration.compute_all(
    area_km2=25.0,
    length_km=8.0,
    height_m=120.0,
    slope_pct=1.5,
)

for method, tc in sorted(results.items(), key=lambda kv: kv[1]):
    print(f"{method:<18} {tc:6.2f} h")
```

## Choosing a Method

| Method | Required inputs | Typical use |
|--------|-----------------|-------------|
| `kirpich` | length, slope | Small, steep rural basins |
| `ventura` | area, slope | Empirical (Italian) general-purpose |
| `giandotti` | area, length, slope | Medium/large rural basins |
| `temez` | length, slope | Natural basins up to ~3000 km² |
| `pasini` | area, length, slope | Empirical (Italian) general-purpose |
| `picking` | length, slope | Small to medium basins |
| `pickering` | length, height | Metric Kirpich-type formula |
| `bransby_williams` | area, length, slope | Small rural basins, overland flow dominant |
| `chpw` | length, height | California Culverts Practice; small dam drainage |
| `ven_te_chow` | length, slope | Rural basins |
| `corps_engineers` | length, slope | Rural basins up to ~11,000 km² |
| `dooge` | area, slope | General-purpose |
| `espey` | area, length, slope | Simplified rural adaptation of Espey-Winslow |

All methods return Tc in **hours**. Inputs follow this convention unless noted in the method's docstring:

- `area_km2` — drainage area (km²)
- `length_km` — main channel / talweg length (km)
- `height_m` — elevation difference between the most distant point and the outlet (m)
- `slope_pct` — average watershed or channel slope (%)

::: pyutils.hydrology.time_concentration.TimeOfConcentration
    options:
      show_signature_annotations: true
