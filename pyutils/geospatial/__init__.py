"""Geospatial analysis module.

Consolidates spatial interpolation, NetCDF I/O using xarray,
shapefile operations, and cartographic visualization.
"""

from .io import XarrayNetCDFHandler
from .interpolation import InverseDistanceWeighting, ThiessenPolygon
from .shapefile import ShapefileHandler
from .visualization import MapRenderer

__all__ = [
    "XarrayNetCDFHandler",
    "InverseDistanceWeighting",
    "ThiessenPolygon",
    "ShapefileHandler",
    "MapRenderer",
]
