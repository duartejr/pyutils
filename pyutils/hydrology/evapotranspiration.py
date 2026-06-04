"""Evapotranspiration calculation models."""

import math
from typing import Optional
import numpy as np
import pyeto.fao as _fao
from pyutils.core import EvapotranspirationModel, SOLAR_CONSTANT, ValidationError


class Hargreaves(EvapotranspirationModel):
    """Hargreaves evapotranspiration model.

    Computes reference evapotranspiration using mean, maximum and minimum
    temperature data. Suitable when only temperature data is available.

    Reference
    ---------
    Hargreaves, G.H., and Z.A. Samani. 1985. Reference crop evapotranspiration
    from temperature. Applied Engineering in Agriculture, 1(2): 96-99.
    """

    def __init__(self):
        """Initialize Hargreaves model."""
        super().__init__(
            name="Hargreaves",
            description="Temperature-based evapotranspiration model",
        )

    def compute(
        self,
        t_med: np.ndarray,
        t_max: np.ndarray,
        t_min: np.ndarray,
        y: float,
        months: int = 1,
    ) -> np.ndarray:
        """Compute reference evapotranspiration.

        Parameters
        ----------
        t_med : np.ndarray
            Mean monthly temperature (°C).
        t_max : np.ndarray
            Maximum monthly temperature (°C).
        t_min : np.ndarray
            Minimum monthly temperature (°C).
        y : float
            Latitude in decimal degrees (positive for N, negative for S).
        months : int
            Number of months to simulate (1-12).

        Returns
        -------
        np.ndarray
            Reference evapotranspiration (mm/month).

        Raises
        ------
        ValidationError
            If temperature data is invalid.
        """
        t_med = np.asarray(t_med)
        t_max = np.asarray(t_max)
        t_min = np.asarray(t_min)

        if not self._validate_temperatures(t_med, t_max, t_min):
            raise ValidationError("Invalid temperature data")

        # Temperature range
        dt = t_max - t_min

        # Solar radiation (simplified)
        ra = self._compute_solar_radiation(y, months)

        # Reference ET
        eto = 0.0023 * ra * (t_med + 17.8) * np.sqrt(dt)

        return eto

    def _validate_temperatures(
        self, t_med: np.ndarray, t_max: np.ndarray, t_min: np.ndarray
    ) -> bool:
        """Validate temperature data consistency."""
        return np.all(t_min <= t_med) and np.all(t_med <= t_max)

    def _compute_solar_radiation(self, y: float, month: int) -> np.ndarray:
        """Compute solar radiation for given latitude and month."""
        # Simplified calculation
        month_array = np.arange(1, 13)
        declination = 23.45 * np.sin(2 * np.pi * (month_array - 81) / 365.25)
        latitude_rad = np.radians(y)

        # Sunset hour angle
        omega = np.arccos(
            -np.tan(latitude_rad) * np.tan(np.radians(declination))
        )

        # Extraterrestrial radiation
        dr = 1 + 0.033 * np.cos(2 * np.pi * month_array / 365.25)
        ra = (24 * 60 / np.pi) * SOLAR_CONSTANT * dr * (
            omega * np.sin(latitude_rad) * np.sin(np.radians(declination))
            + np.cos(latitude_rad)
            * np.cos(np.radians(declination))
            * np.sin(omega)
        )

        return ra[month - 1] if hasattr(ra, "__len__") else ra

    def validate(self) -> bool:
        """Validate model is ready to compute."""
        return True


class Thornthwaite(EvapotranspirationModel):
    """Thornthwaite evapotranspiration model.

    Computes potential evapotranspiration using only temperature data.
    Suitable for large-scale applications with limited data availability.

    Reference
    ---------
    Thornthwaite, C.W. 1948. An approach toward a rational classification of
    climate. Geographical Review, 38(1): 55-94.
    """

    def __init__(self):
        """Initialize Thornthwaite model."""
        super().__init__(
            name="Thornthwaite",
            description="Temperature-based potential evapotranspiration model",
        )

    def compute(
        self,
        t_med: np.ndarray,
        y: float,
    ) -> np.ndarray:
        """Compute potential evapotranspiration.

        Parameters
        ----------
        t_med : np.ndarray
            Mean monthly temperature (°C).
        y : float
            Latitude in decimal degrees.

        Returns
        -------
        np.ndarray
            Potential evapotranspiration (mm/month).

        Raises
        ------
        ValidationError
            If temperature data is invalid.
        """
        t_med = np.asarray(t_med)

        # Filter out negative temperatures
        t_med_positive = np.maximum(t_med, 0)

        # Annual heat index
        i = np.sum(np.power(t_med_positive / 5, 1.514))

        # Exponent a
        a = (
            0.0000675 * i**3
            - 0.0077 * i**2
            + 0.01792 * i
            + 0.49239
        )

        # Monthly PET
        etp = np.zeros_like(t_med_positive)
        for month in range(12):
            if t_med_positive[month] > 0:
                month_index = np.arange(1, 13)
                n = 30  # Days in month (approximate)
                phi = 2 * np.pi * (month + 1) / 12
                latitude_rad = np.radians(y)

                # Day length correction (simplified)
                etp[month] = (
                    16
                    * (n / 30)
                    * np.power(t_med_positive[month] / i, a)
                    * (np.sin(latitude_rad) ** 0.5)
                )

        return etp

    def validate(self) -> bool:
        """Validate model is ready to compute."""
        return True


class PenmanMonteith(EvapotranspirationModel):
    """Penman-Monteith reference evapotranspiration model.

    Computes reference evapotranspiration using comprehensive meteorological data.
    More accurate than temperature-based methods when complete data is available.

    Reference
    ---------
    Allen, R.G., L.S. Pereira, D. Raes, M. Smith, and A.B. Alves. 1998.
    Crop evapotranspiration - Guidelines for computing crop water requirements.
    FAO Irrigation and Drainage Paper 56.
    """

    def __init__(self):
        """Initialize Penman-Monteith model."""
        super().__init__(
            name="Penman-Monteith",
            description="FAO-56 Reference evapotranspiration model",
        )

    def compute(
        self,
        t_mean: float,
        t_max: float,
        t_min: float,
        rh_mean: float,
        u2: float,
        latitude: float,
        day_of_year: int,
        rs: Optional[float] = None,
        elevation: float = 0.0,
    ) -> float:
        """Compute reference evapotranspiration using FAO-56 Penman-Monteith.

        Parameters
        ----------
        t_mean : float
            Mean daily temperature (°C).
        t_max : float
            Maximum daily temperature (°C).
        t_min : float
            Minimum daily temperature (°C).
        rh_mean : float
            Mean relative humidity (%).
        u2 : float
            Wind speed at 2 m height (m/s).
        latitude : float
            Latitude in decimal degrees (positive N, negative S).
        day_of_year : int
            Day of year (1–365).
        rs : float, optional
            Incoming solar radiation (MJ m⁻² day⁻¹). Estimated from
            temperature when not provided.
        elevation : float, optional
            Elevation above sea level (m). Default is 0.

        Returns
        -------
        float
            Reference evapotranspiration ET0 (mm/day).

        Raises
        ------
        ValidationError
            If input parameters are invalid.

        Examples
        --------
        >>> pm = PenmanMonteith()
        >>> et0 = pm.compute(
        ...     t_mean=25.0, t_max=32.0, t_min=18.0,
        ...     rh_mean=60.0, u2=2.0,
        ...     latitude=-15.0, day_of_year=180,
        ...     elevation=500,
        ... )
        >>> round(et0, 1)
        5.6
        """
        if not self._validate_parameters(t_mean, t_max, t_min, rh_mean, u2, elevation):
            raise ValidationError("Invalid Penman-Monteith parameters")

        if not 1 <= day_of_year <= 366:
            raise ValidationError(f"day_of_year must be 1–366, got {day_of_year}")

        # Atmospheric pressure and psychrometric constant
        atm_press = _fao.atm_pressure(elevation)
        psy = _fao.psy_const(atm_press)

        # Saturation and actual vapour pressure
        svp_tmax = _fao.svp_from_t(t_max)
        svp_tmin = _fao.svp_from_t(t_min)
        svp = (svp_tmax + svp_tmin) / 2.0
        avp = _fao.avp_from_rhmean(svp_tmin, svp_tmax, rh_mean)

        # Slope of saturation vapour pressure curve
        delta = _fao.delta_svp(t_mean)

        # Extraterrestrial and clear-sky radiation
        lat_rad = math.radians(latitude)
        sol_dec_val = _fao.sol_dec(day_of_year)
        sha = _fao.sunset_hour_angle(lat_rad, sol_dec_val)
        ird = _fao.inv_rel_dist_earth_sun(day_of_year)
        et_rad_val = _fao.et_rad(lat_rad, sol_dec_val, sha, ird)
        cs_rad_val = _fao.cs_rad(elevation, et_rad_val)

        # Solar radiation — use provided value or estimate from temperature
        if rs is None:
            rs = _fao.sol_rad_from_t(et_rad_val, cs_rad_val, t_min, t_max, coastal=False)

        # Net radiation components (net_out_lw_rad expects Kelvin)
        ni_sw = _fao.net_in_sol_rad(rs)
        no_lw = _fao.net_out_lw_rad(t_min + 273.15, t_max + 273.15, rs, cs_rad_val, avp)
        net_radiation = _fao.net_rad(ni_sw, no_lw)

        # FAO-56 Penman-Monteith (temperature in Kelvin)
        return _fao.fao56_penman_monteith(
            net_radiation, t_mean + 273.15, u2, svp, avp, delta, psy
        )

    def _validate_parameters(
        self,
        t_mean: float,
        t_max: float,
        t_min: float,
        rh_mean: float,
        u2: float,
        elevation: float,
    ) -> bool:
        """Validate Penman-Monteith input parameters."""
        return (
            t_min <= t_mean <= t_max
            and 0 <= rh_mean <= 100
            and u2 >= 0
            and elevation >= -500
        )

    def validate(self) -> bool:
        """Validate model is ready to compute."""
        return True
