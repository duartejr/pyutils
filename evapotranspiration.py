"""
Evapotranspiration calculation methods.

Implements Hargreaves-Samani, Thornthwaite, and Penman-Monteith
reference evapotranspiration (ETP) from basic meteorological inputs.

References
----------
Hargreaves, G.H. and Samani, Z.A. (1985). Reference crop evapotranspiration
    from temperature. Applied Engineering in Agriculture, 1(2), 96-99.
Thornthwaite, C.W. (1948). An approach toward a rational classification of
    climate. Geographical Review, 38(1), 55-94.
Allen, R.G. et al. (1998). FAO Irrigation and Drainage Paper No. 56. FAO, Rome.
"""

from typing import Union

import numpy as np
from numpy import arccos, arctan, cos, exp, log, pi, sin, sqrt, tan

_DAYS_PER_MONTH = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]


def hargreaves(
    t_med: np.ndarray,
    t_max: np.ndarray,
    t_min: np.ndarray,
    y: float,
    months: Union[np.ndarray, int],
    etp_daily: bool = False,
    data_freq: str = "daily",
) -> np.ndarray:
    """Hargreaves-Samani reference evapotranspiration.

    Parameters
    ----------
    t_med : np.ndarray
        Mean air temperature (°C). Monthly series when data_freq='daily',
        or daily series when etp_daily=True.
    t_max : np.ndarray
        Maximum air temperature (°C), same length as t_med.
    t_min : np.ndarray
        Minimum air temperature (°C), same length as t_med.
    y : float
        Latitude in decimal degrees (negative for Southern Hemisphere).
    months : np.ndarray or int
        Array of month numbers (1–12) corresponding to each entry when
        data_freq is not 'daily'. Not used when data_freq='daily'.
    etp_daily : bool, optional
        If True, t_med/t_max/t_min are daily series and output is daily ETP
        (mm/day). Default False.
    data_freq : str, optional
        'daily' (default) — t_med is a monthly series aggregated to daily
        output. Any other value treats t_med as a sparse monthly series
        indexed by *months*.

    Returns
    -------
    np.ndarray
        ETP in mm/day (etp_daily=True) or mm/month (etp_daily=False).
    """
    n = np.arange(1, 366)
    lat = np.radians(y)

    # Solar declination (rad)
    declination = 0.4093 * sin(2 * pi * n / 365 - 1.405)
    # Sunset hour angle (rad)
    w_s = arccos(-tan(lat) * tan(declination))
    # Inverse relative Earth-Sun distance
    d_r = 1 + 0.333 * cos(2 * pi * n / 365)
    # Extraterrestrial radiation (mm/day)
    S_0 = 15.392 * d_r * (
        w_s * sin(lat) * sin(declination)
        + cos(lat) * cos(declination) * sin(w_s)
    )

    n_end = [1] + list(np.cumsum(_DAYS_PER_MONTH))
    # Monthly mean of daily S_0
    S_0m = [
        np.mean(S_0[n_end[i + 1] - _DAYS_PER_MONTH[i]: n_end[i + 1] + 1])
        for i in range(12)
    ]

    if etp_daily:
        while len(S_0) < len(t_med):
            S_0 = np.concatenate((S_0, S_0))
        S_0 = S_0[: len(t_med)]
        HG1 = 0.0023 * (t_med + 17.8) * np.abs(t_max - t_min) ** 0.5 * S_0
        HG1[HG1 < 0] = 0.0
        return np.array(HG1)

    n_anos = len(t_med) // 12
    HG2 = []

    if data_freq == "daily":
        for i in range(n_anos):
            s = slice(12 * i, 12 * i + 12)
            aux = (
                0.0023
                * (t_med[s] + 17.8)
                * np.abs(t_max[s] - t_min[s]) ** 0.5
                * np.array(S_0m)
            )
            aux[aux < 0] = 0
            HG2.append(_DAYS_PER_MONTH * aux)
        return np.array(HG2).reshape(n_anos * 12)
    else:
        for i in range(len(t_med)):
            month = months[i] - 1
            aux = (
                0.0023
                * (t_med[i] + 17.8)
                * np.abs(t_max[i] - t_min[i]) ** 0.5
                * S_0m[month]
            )
            HG2.append(_DAYS_PER_MONTH[month] * max(aux, 0))
        return np.array(HG2)


def thornthwaite(Ti: np.ndarray, n: np.ndarray) -> np.ndarray:
    """Thornthwaite monthly reference evapotranspiration.

    Parameters
    ----------
    Ti : np.ndarray
        Monthly mean temperature series (°C). Length must be a multiple of 12.
    n : np.ndarray
        Actual sunshine hours per day for each month, same length as Ti.

    Returns
    -------
    np.ndarray
        Monthly ETP (mm/month), same length as Ti.
    """
    n_d = np.array(_DAYS_PER_MONTH * (len(Ti) // 12))
    n = n / n_d
    I = []
    for i in range(0, len(Ti), 12):
        I_m = [sum((0.2 * Ti[i: i + 12]) ** 1.514)] * 12
        I.extend(I_m)
    I = np.array(I)
    a = 6.75e-7 * I ** 3 - 7.71e-5 * I ** 2 + 1.7912e-2 * I + 0.49239
    etpp = 16 * (10 * Ti / I) ** a
    etp = etpp * n / 12 * n_d / 30
    return etp


def penman_monteith(
    Tmed=None,
    Tmax=None,
    Tmin=None,
    Rn=None,
    lat=None,
    U2=None,
    Patm=None,
    alt=None,
    UR=None,
    Rs=None,
    n=None,
    albedo: float = 0.23,
    Uz=None,
    G: float = 0.0,
    a: float = 0.25,
    b: float = 0.5,
    z: float = 10,
    date=None,
) -> np.ndarray:
    """FAO-56 Penman-Monteith reference evapotranspiration.

    Parameters
    ----------
    Tmed : array-like
        Mean daily temperature (°C).
    Tmax : array-like
        Maximum daily temperature (°C).
    Tmin : array-like
        Minimum daily temperature (°C).
    Rn : array-like, optional
        Net radiation (MJ m⁻² day⁻¹). Computed internally when not supplied.
    lat : float, optional
        Latitude (decimal degrees). Required when Rn is None.
    U2 : array-like, optional
        Wind speed at 2 m height (m/s). Computed from Uz when not supplied.
    Patm : float, optional
        Atmospheric pressure (kPa). Computed from alt when not supplied.
    alt : float, optional
        Elevation (m). Used to estimate Patm.
    UR : array-like, optional
        Mean relative humidity (%). Used to estimate actual vapour pressure.
    Rs : array-like, optional
        Incoming solar radiation (MJ m⁻² day⁻¹).
    n : array-like, optional
        Daily sunshine hours (h). Used to estimate Rs when Rs is None.
    albedo : float, optional
        Crop reflection coefficient. Default 0.23.
    Uz : array-like, optional
        Wind speed at height z (m/s). Converted to U2 internally.
    G : float, optional
        Soil heat flux density (MJ m⁻² day⁻¹). Default 0.
    a, b : float, optional
        Angström coefficients for Rs estimation. Defaults 0.25, 0.50.
    z : float, optional
        Height of wind measurement (m). Default 10.
    date : array-like, optional
        Sequence of datetime-like objects with a .dayofyear attribute.
        Required when Rn is None.

    Returns
    -------
    np.ndarray
        Reference ETP (mm/day).
    """
    if Patm is None:
        Patm = 101.3 * ((293 - 0.0065 * alt) / 293) ** 5.26

    es = 0.6108 * exp((17.25 * Tmed) / (Tmed + 237.3))

    if isinstance(UR, np.ndarray):
        ea = (es * UR) / 100
    else:
        ea = 0.61 * exp((17.27 * Tmin) / (Tmin + 237.3))

    if Rn is None:
        sigma = 4.903e-9
        J = np.array([x.dayofyear for x in date])
        dr = 1 + 0.033 * cos(J * (2 * pi) / 365)
        Ra = (118.08 / pi) * dr
        lat_rad = np.radians(lat)
        Rso = (0.75 + 2e-5 * alt) * Ra
        ds = 0.409 * sin(((2 * pi) / 365) * J - 1.39)
        X = 1 - (tan(lat_rad)) ** 2 * (tan(ds)) ** 2
        X = np.where(X <= 0, 1e-5, X)
        ws = pi / 2 - arctan((-tan(lat_rad) * tan(ds)) / X ** 0.5)

        if Rs is None:
            N = 24 / pi * ws
            Rs = (a + b * n / N) * Ra

        Rns = (1 - albedo) * Rs
        Ra = 37.6 * dr * (
            ws * sin(lat_rad) * sin(ds) + cos(lat_rad) * cos(ds) * sin(ws)
        )
        Rnl = (
            sigma
            * (((Tmax + 273.16) ** 4 + (Tmin + 273.16) ** 4) / 2)
            * (0.34 - 0.14 * sqrt(ea))
            * (1.35 * (Rs / Rso) - 0.35)
        )
        Rn = Rns - Rnl

    if U2 is None:
        U2 = Uz * 4.87 / log(67.8 * z - 5.42)

    delta = (4098 * (0.6108 * exp((17.27 * Tmed) / (Tmed + 237.3)))) / (Tmed + 237.3) ** 2
    gama = 0.665e-3 * Patm

    pet = (0.408 * delta * (Rn - G) + (gama * 900 * U2 * (es - ea)) / (Tmed + 273)) / (
        delta + gama * (1 + 0.34 * U2)
    )
    return pet
