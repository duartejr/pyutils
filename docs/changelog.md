# Changelog

All notable changes to pyutils are documented here.

---

## [0.4.0] — 2026-06-02

### Added
- **MkDocs documentation site** with Material theme and mkdocstrings auto-generation
- API reference pages for all modules: `evapotranspiration`, `spi`, `idw`, `pyeto.fao`, `pyeto.convert`, `pyeto._check`, `read_nc`, `plot_maps`, `thiessen`, `corr_med`
- Getting-started guides: Installation and Quick Start
- Changelog page
- `docs/` directory with full MkDocs project structure
- `mkdocs.yml` configuration
- Updated `pyproject.toml` with `[project.optional-dependencies] docs` group

---

## [0.3.0] — 2026-05-30

### Added
- Google-style docstrings and type hints for `evapotranspiration.py` and `read_nc.py`
- Comprehensive `README.md` (243 lines) with installation, Quick Start, and module reference

### Fixed
- CI root causes: hardcoded `pyutils-0.1.0.dist-info` path in `test_phase1_verification.py`
- CI `fail-fast: true` was cancelling sibling Python-version jobs — changed to `fail-fast: false`
- `cartopy` build failures no longer block CI (made best-effort)
- Test paths now absolute (using `Path(__file__).parent`) for CI compatibility

---

## [0.2.0] — 2026-05-30

### Added
- GitHub Actions CI workflow (`ci.yml`)
- Auto-release workflow (`release.yml`) — creates a GitHub Release on every PR merge to master
- `pyproject.toml`, `setup.py`, and `__init__.py` for installable Python package

### Removed
- Duplicate files: `plot_maps_old.py`, `anom_pr.py`, `kriging.py`, `test_plot_maps.py`
- Embedded `PyETo-master` directory (31 files) — replaced by `pyeto` package

### Fixed
- `thiessen.py` — `shapely.ops.cascaded_union` → `unary_union` (deprecation)
- `read_nc.py` — removed `print()` debug statement; fixed `df.iloc` deprecation warning
- `bh_thorthwaite.py` — removed broken import of non-existent module
- `corr_med.py` — refactored broken `from … import *` into reusable function

---

## [0.1.0] — 2026-05-29

### Added
- Initial Python package structure: `pyproject.toml`, `setup.py`, `__init__.py`
- Installable as `pip install -e .`
- Registered all 21 Python modules in `[tool.setuptools]`
- Python 3.11 and 3.12 classifiers
- `pyeto` sub-package (FAO-56 reference ETP)
- `rv` sub-package (rainfall distribution correction)
