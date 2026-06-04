"""Pyutils: Geospatial analysis and hydrology utilities for Python

A comprehensive package for hydrological and geospatial analysis including:
- Evapotranspiration calculations (FAO-56 Penman-Monteith, Thornthwaite, Hargreaves)
- Geospatial operations (Thiessen polygons, IDW interpolation, shapefile utilities)
- Hydrological analysis (water balance, random variable modeling, rainfall statistics)
- Cartographic visualization and climate data processing

Usage
-----
>>> import pyutils
>>> from pyutils.hydrology import Hargreaves
>>> from pyutils.geospatial import InverseDistanceWeighting
"""

__version__ = "0.1.0"
__author__ = "Duarte Junior"
__email__ = "duarte.jr105@gmail.com"
__license__ = "GPL-3.0-or-later"

from . import core
from . import hydrology
from . import geospatial
from . import climate
from . import utils

__all__ = [
    "core",
    "hydrology",
    "geospatial",
    "climate",
    "utils",
]
