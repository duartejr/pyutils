"""Abstract base classes for pyutils package."""

from abc import ABC, abstractmethod
from typing import Any, Dict, Optional
import numpy as np
import xarray as xr


class BaseModel(ABC):
    """Base class for all pyutils models.

    Provides common functionality for initialization, validation, and parameter management.
    """

    def __init__(self, name: str, description: str = ""):
        """Initialize base model.

        Parameters
        ----------
        name : str
            Model name or identifier.
        description : str, optional
            Brief description of the model.
        """
        self._name = name
        self._description = description
        self._parameters: Dict[str, Any] = {}

    @property
    def name(self) -> str:
        """Get model name."""
        return self._name

    @property
    def description(self) -> str:
        """Get model description."""
        return self._description

    @property
    def parameters(self) -> Dict[str, Any]:
        """Get model parameters."""
        return self._parameters.copy()

    def set_parameter(self, key: str, value: Any) -> None:
        """Set a model parameter.

        Parameters
        ----------
        key : str
            Parameter name.
        value : Any
            Parameter value.
        """
        self._parameters[key] = value

    @abstractmethod
    def validate(self) -> bool:
        """Validate model state and parameters.

        Returns
        -------
        bool
            True if model is valid, False otherwise.
        """
        pass

    def __repr__(self) -> str:
        """Return string representation."""
        return f"{self.__class__.__name__}(name='{self._name}')"


class EvapotranspirationModel(BaseModel):
    """Abstract base for evapotranspiration models."""

    def __init__(self, name: str, description: str = ""):
        """Initialize evapotranspiration model.

        Parameters
        ----------
        name : str
            Model name.
        description : str, optional
            Model description.
        """
        super().__init__(name, description)

    @abstractmethod
    def compute(self, **kwargs) -> np.ndarray:
        """Compute evapotranspiration.

        Returns
        -------
        np.ndarray
            Evapotranspiration values (mm/day or mm/month).
        """
        pass

    def validate(self) -> bool:
        """Validate model parameters."""
        return True


class WaterBalanceModel(BaseModel):
    """Abstract base for water balance models."""

    def __init__(self, name: str, description: str = ""):
        """Initialize water balance model.

        Parameters
        ----------
        name : str
            Model name.
        description : str, optional
            Model description.
        """
        super().__init__(name, description)
        self._state: Dict[str, Any] = {}

    @property
    def state(self) -> Dict[str, Any]:
        """Get model state variables."""
        return self._state.copy()

    @abstractmethod
    def compute(self, **kwargs) -> Dict[str, np.ndarray]:
        """Compute water balance components.

        Returns
        -------
        dict
            Dictionary with water balance components (runoff, soil_moisture, etc.).
        """
        pass

    def validate(self) -> bool:
        """Validate model parameters and state."""
        return True


class SpatialInterpolator(BaseModel):
    """Abstract base for spatial interpolation methods."""

    def __init__(self, name: str, description: str = ""):
        """Initialize spatial interpolator.

        Parameters
        ----------
        name : str
            Interpolator name.
        description : str, optional
            Interpolator description.
        """
        super().__init__(name, description)

    @abstractmethod
    def interpolate(
        self,
        x: np.ndarray,
        y: np.ndarray,
        z: np.ndarray,
        xi: np.ndarray,
        yi: np.ndarray,
    ) -> np.ndarray:
        """Interpolate values at target locations.

        Parameters
        ----------
        x : np.ndarray
            X-coordinates of source points.
        y : np.ndarray
            Y-coordinates of source points.
        z : np.ndarray
            Values at source points.
        xi : np.ndarray
            X-coordinates of target grid.
        yi : np.ndarray
            Y-coordinates of target grid.

        Returns
        -------
        np.ndarray
            Interpolated values at target locations.
        """
        pass

    def validate(self) -> bool:
        """Validate interpolator parameters."""
        return True


class ClimateDataHandler(BaseModel):
    """Abstract base for climate data I/O and processing."""

    def __init__(self, name: str, description: str = ""):
        """Initialize climate data handler.

        Parameters
        ----------
        name : str
            Handler name.
        description : str, optional
            Handler description.
        """
        super().__init__(name, description)
        self._dataset: Optional[xr.Dataset] = None

    @property
    def dataset(self) -> Optional[xr.Dataset]:
        """Get current dataset."""
        return self._dataset

    @abstractmethod
    def read(self, path: str) -> xr.Dataset:
        """Read climate data from file.

        Parameters
        ----------
        path : str
            Path to data file.

        Returns
        -------
        xr.Dataset
            Loaded dataset.
        """
        pass

    @abstractmethod
    def write(self, path: str, dataset: Optional[xr.Dataset] = None) -> None:
        """Write climate data to file.

        Parameters
        ----------
        path : str
            Output file path.
        dataset : xr.Dataset, optional
            Dataset to write. If None, uses current dataset.
        """
        pass

    @abstractmethod
    def compute_climatology(self) -> xr.Dataset:
        """Compute climatological statistics.

        Returns
        -------
        xr.Dataset
            Dataset with climatological statistics.
        """
        pass

    @abstractmethod
    def compute_anomaly(self, climatology: Optional[xr.Dataset] = None) -> xr.Dataset:
        """Compute anomalies relative to climatology.

        Parameters
        ----------
        climatology : xr.Dataset, optional
            Climatological reference. If None, computes internally.

        Returns
        -------
        xr.Dataset
            Dataset with anomalies.
        """
        pass

    def validate(self) -> bool:
        """Validate that dataset is loaded."""
        return self._dataset is not None


class ClimateAnalysisModel(BaseModel):
    """Abstract base for climate analysis models."""

    def __init__(self, name: str, description: str = ""):
        """Initialize climate analysis model.

        Parameters
        ----------
        name : str
            Model name.
        description : str, optional
            Model description.
        """
        super().__init__(name, description)

    @abstractmethod
    def analyze(self, dataset: xr.Dataset) -> Dict[str, Any]:
        """Perform climate analysis on dataset.

        Parameters
        ----------
        dataset : xr.Dataset
            Input climate data.

        Returns
        -------
        dict
            Analysis results.
        """
        pass

    def validate(self) -> bool:
        """Validate model parameters."""
        return True
