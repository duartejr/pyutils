"""Water balance models consolidating multiple Thornthwaite implementations."""

from typing import Dict, Tuple
import numpy as np
from pyutils.core import WaterBalanceModel, ValidationError


class ThornthwaiteMather(WaterBalanceModel):
    """Thornthwaite-Mather water balance model.

    Computes monthly water balance components including soil moisture,
    runoff, and actual evapotranspiration. Consolidates three previous
    implementations (bh_thorthwaite.py, bh_thorthwaite2.py, ccsthorthwaite.py).

    Reference
    ---------
    Thornthwaite, C.W. and J.R. Mather. 1955. The water balance.
    Publications in Climatology, 8(1): 1-104.
    """

    def __init__(
        self,
        cap: float = 100.0,
        method: str = "standard",
    ):
        """Initialize Thornthwaite-Mather water balance model.

        Parameters
        ----------
        cap : float, optional
            Water holding capacity of soil (mm). Default is 100mm.
        method : str, optional
            Calculation method: 'standard', 'variant2', 'ccs'.
            Different implementations for different regional requirements.
        """
        super().__init__(
            name="Thornthwaite-Mather",
            description="Comprehensive water balance model",
        )
        if cap <= 0:
            raise ValidationError("Water holding capacity must be positive")
        self.cap = cap
        if method not in {"standard", "variant2", "ccs"}:
            raise ValidationError(f"Unknown method: {method}")
        self.method = method
        self._cap_soil = cap

    def compute(
        self,
        precipitation: np.ndarray,
        potential_et: np.ndarray,
        initial_soil_moisture: float = None,
    ) -> Dict[str, np.ndarray]:
        """Compute water balance components.

        Parameters
        ----------
        precipitation : np.ndarray
            Monthly precipitation (mm).
        potential_et : np.ndarray
            Monthly potential evapotranspiration (mm).
        initial_soil_moisture : float, optional
            Initial soil moisture (mm). Default is water holding capacity.

        Returns
        -------
        dict
            Water balance components:
            - 'runoff': Surface runoff (mm)
            - 'actual_et': Actual evapotranspiration (mm)
            - 'soil_moisture': Soil moisture (mm)
            - 'water_deficit': Water deficit (mm)
            - 'water_surplus': Water surplus (mm)

        Raises
        ------
        ValidationError
            If input data is invalid.
        """
        precipitation = np.asarray(precipitation, dtype=float)
        potential_et = np.asarray(potential_et, dtype=float)

        if len(precipitation) != len(potential_et):
            raise ValidationError(
                "precipitation and potential_et must have same length"
            )

        if not np.all(precipitation >= 0):
            raise ValidationError("Precipitation must be non-negative")

        if not np.all(potential_et >= 0):
            raise ValidationError("Potential ET must be non-negative")

        n = len(precipitation)
        if initial_soil_moisture is None:
            initial_soil_moisture = self._cap_soil

        # Initialize arrays
        soil_moisture = np.zeros(n)
        actual_et = np.zeros(n)
        runoff = np.zeros(n)
        water_deficit = np.zeros(n)
        water_surplus = np.zeros(n)

        sm = initial_soil_moisture

        for i in range(n):
            p = precipitation[i]
            pet = potential_et[i]

            if self.method == "standard":
                soil_moisture[i], actual_et[i], runoff[i] = (
                    self._compute_standard(p, pet, sm)
                )
            elif self.method == "variant2":
                soil_moisture[i], actual_et[i], runoff[i] = (
                    self._compute_variant2(p, pet, sm)
                )
            else:  # ccs
                soil_moisture[i], actual_et[i], runoff[i] = (
                    self._compute_ccs(p, pet, sm)
                )

            water_deficit[i] = max(0, pet - actual_et[i])
            water_surplus[i] = max(0, p - actual_et[i] - runoff[i])

            sm = soil_moisture[i]

        self._state = {
            "soil_moisture": soil_moisture[-1],
            "accumulated_runoff": np.sum(runoff),
            "accumulated_aet": np.sum(actual_et),
        }

        return {
            "runoff": runoff,
            "actual_et": actual_et,
            "soil_moisture": soil_moisture,
            "water_deficit": water_deficit,
            "water_surplus": water_surplus,
        }

    def _compute_standard(
        self, p: float, pet: float, sm_previous: float
    ) -> Tuple[float, float, float]:
        """Standard Thornthwaite-Mather calculation."""
        # Available water = previous soil moisture + precipitation
        available = sm_previous + p

        # Check if we can meet potential ET
        if available >= pet:
            actual_et = pet
            sm = available - pet
            runoff = 0
        else:
            # Water is limited
            actual_et = available
            sm = 0
            runoff = 0

        # Limit soil moisture to capacity
        if sm > self._cap_soil:
            runoff = sm - self._cap_soil
            sm = self._cap_soil

        return sm, actual_et, runoff

    def _compute_variant2(
        self, p: float, pet: float, sm_previous: float
    ) -> Tuple[float, float, float]:
        """Variant 2 calculation with different runoff threshold."""
        available = sm_previous + p

        if available >= pet:
            actual_et = pet
            sm = available - pet
        else:
            actual_et = available
            sm = 0

        # Higher runoff threshold
        runoff = 0
        if sm > self._cap_soil * 1.1:
            runoff = sm - self._cap_soil * 1.1
            sm = self._cap_soil * 1.1

        return sm, actual_et, runoff

    def _compute_ccs(
        self, p: float, pet: float, sm_previous: float
    ) -> Tuple[float, float, float]:
        """CCS variant calculation with storage capacity dynamics."""
        available = sm_previous + p

        if available >= pet:
            actual_et = pet
            sm = available - pet
        else:
            # Reduced ET under water stress
            stress_factor = available / pet if pet > 0 else 0
            actual_et = available * stress_factor
            sm = available * (1 - stress_factor)

        # Variable capacity based on vegetation
        effective_cap = self._cap_soil * 0.95
        runoff = 0
        if sm > effective_cap:
            runoff = sm - effective_cap
            sm = effective_cap

        return sm, actual_et, runoff

    def validate(self) -> bool:
        """Validate model parameters."""
        return self.cap > 0 and self.method in {"standard", "variant2", "ccs"}
