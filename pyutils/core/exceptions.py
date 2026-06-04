"""Custom exceptions for pyutils package."""


class PyutilsError(Exception):
    """Base exception for all pyutils errors."""

    pass


class ValidationError(PyutilsError):
    """Raised when input validation fails."""

    pass


class DataError(PyutilsError):
    """Raised when data processing fails."""

    pass


class InterpolationError(PyutilsError):
    """Raised when interpolation fails."""

    pass


class GeospatialError(PyutilsError):
    """Raised when geospatial operation fails."""

    pass


class ClimateError(PyutilsError):
    """Raised when climate analysis fails."""

    pass
