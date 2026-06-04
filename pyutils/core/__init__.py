"""Core infrastructure for pyutils package.

Provides base classes, exceptions, and constants for all pyutils modules.
"""

from .base import (
    BaseModel,
    EvapotranspirationModel,
    WaterBalanceModel,
    SpatialInterpolator,
    ClimateDataHandler,
    ClimateAnalysisModel,
)
from .exceptions import (
    PyutilsError,
    ValidationError,
    DataError,
    InterpolationError,
    GeospatialError,
    ClimateError,
)
from .constants import (
    SOLAR_CONSTANT,
    STEFAN_BOLTZMANN,
    LATENT_HEAT_VAPORIZATION,
    PSYCHROMETRIC_CONSTANT,
    EARTH_RADIUS,
    DEGREES_TO_RADIANS,
    RADIANS_TO_DEGREES,
)

__all__ = [
    "BaseModel",
    "EvapotranspirationModel",
    "WaterBalanceModel",
    "SpatialInterpolator",
    "ClimateDataHandler",
    "ClimateAnalysisModel",
    "PyutilsError",
    "ValidationError",
    "DataError",
    "InterpolationError",
    "GeospatialError",
    "ClimateError",
    "SOLAR_CONSTANT",
    "STEFAN_BOLTZMANN",
    "LATENT_HEAT_VAPORIZATION",
    "PSYCHROMETRIC_CONSTANT",
    "EARTH_RADIUS",
    "DEGREES_TO_RADIANS",
    "RADIANS_TO_DEGREES",
]
