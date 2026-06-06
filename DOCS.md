# Documentation Structure

This document maps the documentation to ensure all classes and modules are discoverable.

## File Organization

```
docs/
├── index.md                          # Overview, features, installation link
├── changelog.md                      # Release history
├── getting-started/
│   ├── installation.md               # System requirements, pip install
│   └── quickstart.md                 # Code examples for each module
└── api/
    ├── hydrology/
    │   ├── evapotranspiration.md     # Hargreaves, Thornthwaite, PenmanMonteith
    │   ├── water_balance.md          # ThornthwaiteMather
    │   ├── indices.md                # StandardizedPrecipitationIndex (SPI)
    │   └── flow_analysis.md          # FlowAnalyzer
    ├── geospatial/
    │   ├── io.md                     # XarrayNetCDFHandler
    │   ├── interpolation.md          # InverseDistanceWeighting, ThiessenPolygon
    │   ├── shapefile.md              # ShapefileHandler
    │   └── visualization.md          # MapRenderer
    ├── climate/
    │   └── bias_correction.md        # BiasCorrection
    └── utils/
        ├── validators.md             # DataValidator
        └── conversions.md            # UnitConverter
```

## Class Coverage Checklist

| Class | File | Status |
|-------|------|--------|
| `Hargreaves` | `api/hydrology/evapotranspiration.md` | ✓ Documented |
| `Thornthwaite` | `api/hydrology/evapotranspiration.md` | ✓ Documented |
| `PenmanMonteith` | `api/hydrology/evapotranspiration.md` | ✓ Documented |
| `ThornthwaiteMather` | `api/hydrology/water_balance.md` | ✓ Documented |
| `StandardizedPrecipitationIndex` | `api/hydrology/indices.md` | ✓ Documented |
| `FlowAnalyzer` | `api/hydrology/flow_analysis.md` | ✓ Documented |
| `XarrayNetCDFHandler` | `api/geospatial/io.md` | ✓ Documented |
| `InverseDistanceWeighting` | `api/geospatial/interpolation.md` | ✓ Documented |
| `ThiessenPolygon` | `api/geospatial/interpolation.md` | ✓ Documented |
| `ShapefileHandler` | `api/geospatial/shapefile.md` | ✓ Documented |
| `MapRenderer` | `api/geospatial/visualization.md` | ✓ Documented |
| `BiasCorrection` | `api/climate/bias_correction.md` | ✓ Documented |
| `DataValidator` | `api/utils/validators.md` | ✓ Documented |
| `UnitConverter` | `api/utils/conversions.md` | ✓ Documented |

## Navigation Structure (mkdocs.yml)

All files are properly registered in the sidebar navigation:
- **Home** → `index.md`
- **Getting Started**
  - Installation → `getting-started/installation.md`
  - Quick Start → `getting-started/quickstart.md`
- **API Reference** (organized by sub-package)
  - Hydrology (4 pages)
  - Geospatial (4 pages)
  - Climate (1 page)
  - Utils (2 pages)
- **Changelog** → `changelog.md`

## Build Status

- ✅ `mkdocs build --strict` passes with zero warnings
- ✅ All 15 markdown files build successfully
- ✅ All internal links are valid
- ✅ No broken references in nav
- ✅ All classes mentioned in index.md have dedicated pages

## Maintenance Guidelines

When adding new classes:
1. Create the markdown doc file in the appropriate `api/*/` subdirectory
2. Add a reference to `mkdocs.yml` nav under the appropriate section
3. Update `docs/index.md` overview table if it's a top-level class
4. Run `mkdocs build --strict` to verify no warnings
5. Check that the class name appears in its dedicated doc page
