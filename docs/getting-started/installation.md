# Installation

## Requirements

- **Python 3.11** or later
- For geospatial modules (`plot_maps`, `thiessen`): GEOS and PROJ system libraries

## Clone and install

```bash
git clone https://github.com/duartejr/pyutils.git
cd pyutils
pip install -e .
```

The `-e` flag installs the package in editable mode so local changes take effect immediately.

## System dependencies (geospatial)

Some modules depend on C extensions that require system libraries. Install them before running `pip install`:

=== "Ubuntu / Debian"
    ```bash
    sudo apt-get install libgeos-dev libproj-dev proj-data proj-bin
    ```

=== "macOS (Homebrew)"
    ```bash
    brew install geos proj
    ```

=== "Conda"
    ```bash
    conda install -c conda-forge cartopy geopandas
    ```

## Optional extras

```bash
# Install dev tools (pytest, black, mypy, flake8)
pip install -e ".[dev]"

# Install cartopy visualization support
pip install -e ".[viz]"
```

## Verifying the install

```python
import pyutils
import idw
import spi
import evapotranspiration

print("pyutils ready")
```

If you see no import errors, the core installation succeeded.

!!! note "NetCDF4 module"
    `read_nc` requires the `netCDF4` package which needs HDF5 system libraries.
    If you only need SPI, IDW, or evapotranspiration, you can skip this.
