"""Unit conversion utilities for meteorological and hydrological quantities."""

import math
import numpy as np


class UnitConverter:
    """Convert between common units used in hydroclimate analysis.

    All methods are static — no instantiation needed.

    Examples
    --------
    >>> UnitConverter.celsius_to_kelvin(25.0)
    298.15
    >>> UnitConverter.degrees_to_radians(180.0)
    3.141592653589793
    """

    # ------------------------------------------------------------------ #
    # Temperature                                                          #
    # ------------------------------------------------------------------ #

    @staticmethod
    def celsius_to_kelvin(temperature: float) -> float:
        """Convert Celsius to Kelvin.

        Parameters
        ----------
        temperature : float
            Temperature in degrees Celsius.

        Returns
        -------
        float
            Temperature in Kelvin.

        Examples
        --------
        >>> UnitConverter.celsius_to_kelvin(0.0)
        273.15
        >>> UnitConverter.celsius_to_kelvin(25.0)
        298.15
        """
        return temperature + 273.15

    @staticmethod
    def kelvin_to_celsius(temperature: float) -> float:
        """Convert Kelvin to Celsius.

        Parameters
        ----------
        temperature : float
            Temperature in Kelvin.

        Returns
        -------
        float
            Temperature in degrees Celsius.

        Examples
        --------
        >>> UnitConverter.kelvin_to_celsius(273.15)
        0.0
        """
        return temperature - 273.15

    @staticmethod
    def fahrenheit_to_celsius(temperature: float) -> float:
        """Convert Fahrenheit to Celsius.

        Parameters
        ----------
        temperature : float
            Temperature in degrees Fahrenheit.

        Returns
        -------
        float
            Temperature in degrees Celsius.

        Examples
        --------
        >>> UnitConverter.fahrenheit_to_celsius(32.0)
        0.0
        >>> UnitConverter.fahrenheit_to_celsius(212.0)
        100.0
        """
        return (temperature - 32.0) * 5.0 / 9.0

    @staticmethod
    def celsius_to_fahrenheit(temperature: float) -> float:
        """Convert Celsius to Fahrenheit.

        Parameters
        ----------
        temperature : float
            Temperature in degrees Celsius.

        Returns
        -------
        float
            Temperature in degrees Fahrenheit.

        Examples
        --------
        >>> UnitConverter.celsius_to_fahrenheit(0.0)
        32.0
        >>> UnitConverter.celsius_to_fahrenheit(100.0)
        212.0
        """
        return temperature * 9.0 / 5.0 + 32.0

    # ------------------------------------------------------------------ #
    # Angles                                                               #
    # ------------------------------------------------------------------ #

    @staticmethod
    def degrees_to_radians(degrees: float) -> float:
        """Convert decimal degrees to radians.

        Parameters
        ----------
        degrees : float
            Angle in decimal degrees.

        Returns
        -------
        float
            Angle in radians.

        Examples
        --------
        >>> round(UnitConverter.degrees_to_radians(180.0), 6)
        3.141593
        """
        return degrees * math.pi / 180.0

    @staticmethod
    def radians_to_degrees(radians: float) -> float:
        """Convert radians to decimal degrees.

        Parameters
        ----------
        radians : float
            Angle in radians.

        Returns
        -------
        float
            Angle in decimal degrees.

        Examples
        --------
        >>> round(UnitConverter.radians_to_degrees(math.pi), 6)
        180.0
        """
        return radians * 180.0 / math.pi

    # ------------------------------------------------------------------ #
    # Wind speed                                                           #
    # ------------------------------------------------------------------ #

    @staticmethod
    def wind_speed_height_correction(
        wind_speed: float,
        measurement_height: float,
        target_height: float = 2.0,
    ) -> float:
        """Correct wind speed from measurement height to target height.

        Uses the logarithmic wind profile formula (FAO-56 Eq. 47).

        Parameters
        ----------
        wind_speed : float
            Measured wind speed (m/s).
        measurement_height : float
            Height at which wind was measured (m).
        target_height : float, optional
            Target height (m). Default is 2 m (FAO standard).

        Returns
        -------
        float
            Wind speed at target height (m/s).

        Examples
        --------
        >>> round(UnitConverter.wind_speed_height_correction(5.0, 10.0, 2.0), 4)
        3.8656
        """
        return wind_speed * (4.87 / math.log(67.8 * measurement_height - 5.42))

    # ------------------------------------------------------------------ #
    # Precipitation / depth                                                #
    # ------------------------------------------------------------------ #

    @staticmethod
    def mm_to_meters(depth_mm: float) -> float:
        """Convert millimetres to metres.

        Parameters
        ----------
        depth_mm : float
            Depth in millimetres.

        Returns
        -------
        float
            Depth in metres.

        Examples
        --------
        >>> UnitConverter.mm_to_meters(1000.0)
        1.0
        """
        return depth_mm / 1000.0

    @staticmethod
    def meters_to_mm(depth_m: float) -> float:
        """Convert metres to millimetres.

        Parameters
        ----------
        depth_m : float
            Depth in metres.

        Returns
        -------
        float
            Depth in millimetres.

        Examples
        --------
        >>> UnitConverter.meters_to_mm(1.0)
        1000.0
        """
        return depth_m * 1000.0

    # ------------------------------------------------------------------ #
    # Radiation                                                            #
    # ------------------------------------------------------------------ #

    @staticmethod
    def mj_per_m2_to_watts_per_m2(radiation_mj: float) -> float:
        """Convert MJ m⁻² day⁻¹ to W m⁻².

        Parameters
        ----------
        radiation_mj : float
            Radiation in MJ m⁻² day⁻¹.

        Returns
        -------
        float
            Radiation in W m⁻².

        Examples
        --------
        >>> round(UnitConverter.mj_per_m2_to_watts_per_m2(1.0), 4)
        11.5741
        """
        return radiation_mj * 1e6 / (24 * 3600)

    @staticmethod
    def watts_per_m2_to_mj_per_m2(radiation_w: float) -> float:
        """Convert W m⁻² to MJ m⁻² day⁻¹.

        Parameters
        ----------
        radiation_w : float
            Radiation in W m⁻².

        Returns
        -------
        float
            Radiation in MJ m⁻² day⁻¹.

        Examples
        --------
        >>> round(UnitConverter.watts_per_m2_to_mj_per_m2(11.5741), 2)
        1.0
        """
        return radiation_w * (24 * 3600) / 1e6
