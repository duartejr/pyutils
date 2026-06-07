"""Time of concentration estimation using empirical formulas."""

from typing import Dict
import numpy as np
from pyutils.core import ValidationError


class TimeOfConcentration:
    """Estimate watershed time of concentration with empirical formulas.

    Time of concentration (Tc) is the time required for runoff to travel
    from the most hydraulically distant point of a watershed to its outlet.
    It is a key input for the rational method, flood routing, and drainage
    design.

    Thirteen widely used empirical formulas are provided. Each one was
    calibrated on different watershed types and sizes, so results can vary
    considerably between methods for the same watershed; ``compute_all``
    helps compare them side by side.

    Notes
    -----
    Unless noted otherwise, formulas expect:

    - ``area_km2``: drainage area (km^2)
    - ``length_km``: main channel / talweg length (km)
    - ``height_m``: elevation difference between the most distant point and
      the outlet (m)
    - ``slope_pct``: average watershed or channel slope (%)

    All methods return Tc in hours.
    """

    @staticmethod
    def kirpich(length_km: float, slope_pct: float) -> float:
        """Estimate Tc with the Kirpich (1940) formula.

        Suitable for small, steep rural watersheds (originally calibrated on
        basins of 0.004-0.5 km^2 in Tennessee, USA).

        Parameters
        ----------
        length_km : float
            Main channel length (km).
        slope_pct : float
            Average channel slope (%).

        Returns
        -------
        float
            Time of concentration (hours).

        Raises
        ------
        ValidationError
            If inputs are not positive.
        """
        TimeOfConcentration._validate_positive(length_km=length_km, slope_pct=slope_pct)
        slope = slope_pct / 100.0
        return 0.0663 * length_km**0.77 * slope**-0.385

    @staticmethod
    def ventura(area_km2: float, slope_pct: float) -> float:
        """Estimate Tc with the Ventura formula.

        Empirical Italian formula relating Tc to drainage area and slope.

        Parameters
        ----------
        area_km2 : float
            Drainage area (km^2).
        slope_pct : float
            Average watershed slope (%).

        Returns
        -------
        float
            Time of concentration (hours).

        Raises
        ------
        ValidationError
            If inputs are not positive.
        """
        TimeOfConcentration._validate_positive(area_km2=area_km2, slope_pct=slope_pct)
        slope = slope_pct / 100.0
        return 0.127 * area_km2**0.5 * slope**-0.5

    @staticmethod
    def giandotti(area_km2: float, length_km: float, slope_pct: float) -> float:
        """Estimate Tc with the Giandotti (1934) formula.

        Widely used for medium to large rural basins, particularly in
        Mediterranean climates.

        Parameters
        ----------
        area_km2 : float
            Drainage area (km^2).
        length_km : float
            Main channel length (km).
        slope_pct : float
            Average watershed slope (%).

        Returns
        -------
        float
            Time of concentration (hours).

        Raises
        ------
        ValidationError
            If inputs are not positive.
        """
        TimeOfConcentration._validate_positive(
            area_km2=area_km2, length_km=length_km, slope_pct=slope_pct
        )
        slope = slope_pct / 100.0
        return (
            0.0559
            * (4.0 * area_km2**0.5 + 1.5 * length_km)
            * length_km**-0.5
            * slope**-0.5
        )

    @staticmethod
    def temez(length_km: float, slope_pct: float) -> float:
        """Estimate Tc with the Temez (1978) formula.

        Recommended by the Spanish drainage manual (IEP) for natural basins
        up to roughly 3000 km^2.

        Parameters
        ----------
        length_km : float
            Main channel length (km).
        slope_pct : float
            Average channel slope (%).

        Returns
        -------
        float
            Time of concentration (hours).

        Raises
        ------
        ValidationError
            If inputs are not positive.
        """
        TimeOfConcentration._validate_positive(length_km=length_km, slope_pct=slope_pct)
        slope = slope_pct / 100.0
        return 0.30 * (length_km / slope**0.25) ** 0.76

    @staticmethod
    def pasini(area_km2: float, length_km: float, slope_pct: float) -> float:
        """Estimate Tc with the Pasini (1914) formula.

        Empirical Italian formula combining area, channel length, and slope.

        Parameters
        ----------
        area_km2 : float
            Drainage area (km^2).
        length_km : float
            Main channel length (km).
        slope_pct : float
            Average watershed slope (%).

        Returns
        -------
        float
            Time of concentration (hours).

        Raises
        ------
        ValidationError
            If inputs are not positive.
        """
        TimeOfConcentration._validate_positive(
            area_km2=area_km2, length_km=length_km, slope_pct=slope_pct
        )
        slope = slope_pct / 100.0
        return 0.107 * (area_km2 * length_km) ** (1.0 / 3.0) * slope**-0.5

    @staticmethod
    def picking(length_km: float, slope_pct: float) -> float:
        """Estimate Tc with the Picking formula.

        Parameters
        ----------
        length_km : float
            Main channel length (km).
        slope_pct : float
            Average channel slope (%).

        Returns
        -------
        float
            Time of concentration (hours).

        Raises
        ------
        ValidationError
            If inputs are not positive.
        """
        TimeOfConcentration._validate_positive(length_km=length_km, slope_pct=slope_pct)
        slope = slope_pct / 100.0
        return 0.0883 * length_km**0.667 * slope**-0.333

    @staticmethod
    def pickering(length_km: float, height_m: float) -> float:
        """Estimate Tc with the Pickering formula.

        Metric adaptation of the Kirpich equation based on channel length
        and elevation drop rather than slope.

        Parameters
        ----------
        length_km : float
            Main channel length (km).
        height_m : float
            Elevation difference between the most distant point and the
            outlet (m).

        Returns
        -------
        float
            Time of concentration (hours).

        Raises
        ------
        ValidationError
            If inputs are not positive.
        """
        TimeOfConcentration._validate_positive(length_km=length_km, height_m=height_m)
        return (0.871 * length_km**3 / height_m) ** 0.385

    @staticmethod
    def bransby_williams(area_km2: float, length_km: float, slope_pct: float) -> float:
        """Estimate Tc with the Bransby Williams (1922) formula.

        Suited to small rural catchments without well-defined drainage
        channels, where overland flow dominates. The original formula uses
        area in hectares and slope in m/km; both are converted internally.

        Parameters
        ----------
        area_km2 : float
            Drainage area (km^2).
        length_km : float
            Main channel length (km).
        slope_pct : float
            Average channel slope (%).

        Returns
        -------
        float
            Time of concentration (hours).

        Raises
        ------
        ValidationError
            If inputs are not positive.
        """
        TimeOfConcentration._validate_positive(
            area_km2=area_km2, length_km=length_km, slope_pct=slope_pct
        )
        area_ha = area_km2 * 100.0
        slope_m_per_km = slope_pct * 10.0
        tc_minutes = 92.5 * length_km / (area_ha**0.1 * slope_m_per_km**0.2)
        return tc_minutes / 60.0

    @staticmethod
    def chpw(length_km: float, height_m: float) -> float:
        """Estimate Tc with the CHPW (California Culverts Practice) formula.

        Published in *California Highways and Public Works* (1942); a
        metric, Kirpich-type formula based on channel length and elevation
        drop. Adopted in Brazil (DER/SP) for small dam drainage design.

        Parameters
        ----------
        length_km : float
            Main channel length (km).
        height_m : float
            Elevation difference between the most distant point and the
            outlet (m).

        Returns
        -------
        float
            Time of concentration (hours).

        Raises
        ------
        ValidationError
            If inputs are not positive.
        """
        TimeOfConcentration._validate_positive(length_km=length_km, height_m=height_m)
        return 0.9608 * (length_km**3 / height_m) ** 0.385

    @staticmethod
    def ven_te_chow(length_km: float, slope_pct: float) -> float:
        """Estimate Tc with the Ven Te Chow (1962) formula.

        Recommended for rural basins.

        Parameters
        ----------
        length_km : float
            Main channel length (km).
        slope_pct : float
            Average channel slope (%).

        Returns
        -------
        float
            Time of concentration (hours).

        Raises
        ------
        ValidationError
            If inputs are not positive.
        """
        TimeOfConcentration._validate_positive(length_km=length_km, slope_pct=slope_pct)
        slope = slope_pct / 100.0
        return 0.160 * length_km**0.64 * slope**-0.32

    @staticmethod
    def corps_engineers(length_km: float, slope_pct: float) -> float:
        """Estimate Tc with the US Army Corps of Engineers formula.

        Indicated for rural basins of up to roughly 11,000 km^2
        (Silveira, 2005).

        Parameters
        ----------
        length_km : float
            Main channel length (km).
        slope_pct : float
            Average channel slope (%).

        Returns
        -------
        float
            Time of concentration (hours).

        Raises
        ------
        ValidationError
            If inputs are not positive.
        """
        TimeOfConcentration._validate_positive(length_km=length_km, slope_pct=slope_pct)
        slope = slope_pct / 100.0
        return 0.191 * length_km**0.76 * slope**-0.19

    @staticmethod
    def dooge(area_km2: float, slope_pct: float) -> float:
        """Estimate Tc with the Dooge (1973) formula.

        Parameters
        ----------
        area_km2 : float
            Drainage area (km^2).
        slope_pct : float
            Average watershed slope (%).

        Returns
        -------
        float
            Time of concentration (hours).

        Raises
        ------
        ValidationError
            If inputs are not positive.
        """
        TimeOfConcentration._validate_positive(area_km2=area_km2, slope_pct=slope_pct)
        slope = slope_pct / 100.0
        return 0.365 * area_km2**0.41 * slope**-0.17

    @staticmethod
    def espey(area_km2: float, length_km: float, slope_pct: float) -> float:
        """Estimate Tc with the Espey-Winslow (1968) formula.

        Simplified rural-basin adaptation that drops the original
        impervious-area term (assumes a natural, non-urbanized basin).

        Parameters
        ----------
        area_km2 : float
            Drainage area (km^2).
        length_km : float
            Main channel length (km).
        slope_pct : float
            Average channel slope (%).

        Returns
        -------
        float
            Time of concentration (hours).

        Raises
        ------
        ValidationError
            If inputs are not positive.
        """
        TimeOfConcentration._validate_positive(
            area_km2=area_km2, length_km=length_km, slope_pct=slope_pct
        )
        slope = slope_pct / 100.0
        return 0.108 * (area_km2 * length_km) ** 0.33 * slope**-0.5

    @classmethod
    def compute_all(
        cls,
        area_km2: float,
        length_km: float,
        height_m: float,
        slope_pct: float,
    ) -> Dict[str, float]:
        """Compute Tc with every available method.

        Parameters
        ----------
        area_km2 : float
            Drainage area (km^2).
        length_km : float
            Main channel length (km).
        height_m : float
            Elevation difference between the most distant point and the
            outlet (m).
        slope_pct : float
            Average watershed/channel slope (%).

        Returns
        -------
        dict
            Mapping of method name to time of concentration (hours).
        """
        return {
            "Kirpich": cls.kirpich(length_km, slope_pct),
            "Ventura": cls.ventura(area_km2, slope_pct),
            "Giandotti": cls.giandotti(area_km2, length_km, slope_pct),
            "Temez": cls.temez(length_km, slope_pct),
            "Pasini": cls.pasini(area_km2, length_km, slope_pct),
            "Picking": cls.picking(length_km, slope_pct),
            "Pickering": cls.pickering(length_km, height_m),
            "Bransby Williams": cls.bransby_williams(area_km2, length_km, slope_pct),
            "CHPW": cls.chpw(length_km, height_m),
            "Ven te Chow": cls.ven_te_chow(length_km, slope_pct),
            "Corps Engineers": cls.corps_engineers(length_km, slope_pct),
            "Dooge": cls.dooge(area_km2, slope_pct),
            "Espey": cls.espey(area_km2, length_km, slope_pct),
        }

    @staticmethod
    def _validate_positive(**values: float) -> None:
        """Validate that all given values are finite and strictly positive."""
        for name, value in values.items():
            if not np.isfinite(value) or value <= 0:
                raise ValidationError(
                    f"'{name}' must be a positive, finite number, got {value}"
                )
