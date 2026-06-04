"""Hydrological analysis module.

Consolidates evapotranspiration models, water balance calculations,
and precipitation-based drought indices.
"""

from .evapotranspiration import Hargreaves, Thornthwaite, PenmanMonteith
from .water_balance import ThornthwaiteMather
from .indices import StandardizedPrecipitationIndex

__all__ = [
    "Hargreaves",
    "Thornthwaite",
    "PenmanMonteith",
    "ThornthwaiteMather",
    "StandardizedPrecipitationIndex",
]
