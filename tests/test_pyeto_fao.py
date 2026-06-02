"""Unit tests for FAO-56 evapotranspiration calculation functions.

Tests cover atmospheric pressure, vapour pressure, solar radiation,
and temperature estimation functions from the FAO-56 methodology.
"""

import math

import pytest

from pyeto.fao import (
    atm_pressure,
    avp_from_rhmax,
    avp_from_rhmin_rhmax,
    avp_from_rhmean,
    avp_from_tdew,
    avp_from_tmin,
    cs_rad,
    daily_mean_t,
    daylight_hours,
    delta_svp,
    energy2evap,
    fao56_penman_monteith,
    hargreaves,
    inv_rel_dist_earth_sun,
    mean_svp,
    net_in_sol_rad,
)


class TestAtmosphericPressure:
    """Test atmospheric pressure estimation from altitude."""

    def test_atm_pressure_sea_level(self):
        """Atmospheric pressure at sea level should be ~101.3 kPa."""
        result = atm_pressure(0)
        assert result == pytest.approx(101.3, rel=0.01)

    def test_atm_pressure_1000m(self):
        """Atmospheric pressure decreases with altitude."""
        result_sea_level = atm_pressure(0)
        result_1000m = atm_pressure(1000)
        assert result_1000m < result_sea_level

    def test_atm_pressure_high_elevation(self):
        """High elevation (3000 m) should have much lower pressure."""
        result = atm_pressure(3000)
        assert result < atm_pressure(0)

    def test_atm_pressure_positive(self):
        """Atmospheric pressure should always be positive."""
        for elevation in [0, 500, 1000, 2000, 5000]:
            result = atm_pressure(elevation)
            assert result > 0


class TestVapourPressureFromTemperature:
    """Test vapour pressure calculations from temperature data."""

    def test_avp_from_tmin_zero_celsius(self):
        """Actual vapour pressure at 0°C should be ~0.611 kPa."""
        result = avp_from_tmin(0)
        assert result == pytest.approx(0.611, rel=0.01)

    def test_avp_from_tmin_25_celsius(self):
        """Actual vapour pressure at 25°C should be realistic."""
        result = avp_from_tmin(25)
        assert 3 < result < 4

    def test_avp_from_tmin_increases_with_temperature(self):
        """Vapour pressure increases with temperature."""
        result_10 = avp_from_tmin(10)
        result_20 = avp_from_tmin(20)
        assert result_20 > result_10

    def test_avp_from_tdew_equals_tmin(self):
        """Dew point vapour pressure at 15°C should be realistic."""
        result = avp_from_tdew(15)
        assert 1.5 < result < 2

    def test_avp_from_tdew_positive(self):
        """Dew point vapour pressure should be positive."""
        result = avp_from_tdew(10)
        assert result > 0


class TestVapourPressureFromHumidity:
    """Test actual vapour pressure from relative humidity data."""

    def test_avp_from_rhmax_100_percent(self):
        """At 100% RH, actual VP should equal saturation VP."""
        svp_tmin = 2.0  # Saturation vapour pressure at min temp
        result = avp_from_rhmax(svp_tmin, 100)
        assert result == pytest.approx(2.0)

    def test_avp_from_rhmax_50_percent(self):
        """At 50% RH, actual VP should be half saturation VP."""
        svp_tmin = 2.0
        result = avp_from_rhmax(svp_tmin, 50)
        assert result == pytest.approx(1.0)

    def test_avp_from_rhmax_zero_percent(self):
        """At 0% RH, actual VP should be zero."""
        svp_tmin = 2.0
        result = avp_from_rhmax(svp_tmin, 0)
        assert result == pytest.approx(0)

    def test_avp_from_rhmin_rhmax_symmetric(self):
        """Equal min/max RH should give same result as uniform RH."""
        svp_tmin = 2.0
        svp_tmax = 3.0
        result = avp_from_rhmin_rhmax(svp_tmin, svp_tmax, 60, 60)
        expected = (svp_tmin * 60 + svp_tmax * 60) / 100 / 2
        assert result == pytest.approx(expected)

    def test_avp_from_rhmean_50_percent(self):
        """At 50% mean RH, actual VP should be half saturation VP."""
        svp_tmin = 2.0
        svp_tmax = 3.0
        result = avp_from_rhmean(svp_tmin, svp_tmax, 50)
        expected = ((svp_tmin + svp_tmax) / 2) * 0.5
        assert result == pytest.approx(expected)


class TestSaturationVapourPressure:
    """Test saturation vapour pressure calculations."""

    def test_mean_svp_zero_celsius(self):
        """Mean SVP at 0°C should be realistic."""
        result = mean_svp(0, 0)
        assert 0.6 < result < 0.65

    def test_mean_svp_25_celsius(self):
        """Mean SVP at 25°C should be realistic."""
        result = mean_svp(25, 25)
        assert 3.1 < result < 3.3

    def test_mean_svp_increases_with_temperature(self):
        """Mean SVP increases with temperature."""
        result_10 = mean_svp(10, 10)
        result_20 = mean_svp(20, 20)
        assert result_20 > result_10

    def test_delta_svp_positive(self):
        """SVP slope (delta) should be positive."""
        for temp in [0, 10, 20, 30]:
            result = delta_svp(temp)
            assert result > 0

    def test_delta_svp_increases_with_temperature(self):
        """SVP slope increases with temperature."""
        delta_10 = delta_svp(10)
        delta_20 = delta_svp(20)
        assert delta_20 > delta_10


class TestTemperatureCalculations:
    """Test daily temperature-derived calculations."""

    def test_daily_mean_t_arithmetic_mean(self):
        """Daily mean should be arithmetic mean of max/min."""
        result = daily_mean_t(20, 30)
        assert result == pytest.approx(25)

    def test_daily_mean_t_equal_temps(self):
        """When Tmax equals Tmin, mean should equal both."""
        result = daily_mean_t(25, 25)
        assert result == pytest.approx(25)

    def test_daily_mean_t_large_range(self):
        """Large temperature range should be handled."""
        result = daily_mean_t(10, 40)
        assert result == pytest.approx(25)


class TestSolarRadiation:
    """Test solar radiation-related calculations."""

    def test_inv_rel_dist_earth_sun_jan_first(self):
        """Inverse relative Earth-Sun distance on Jan 1 (~1.033)."""
        result = inv_rel_dist_earth_sun(1)
        assert 1.03 < result < 1.04

    def test_inv_rel_dist_earth_sun_jan_4_perihelion(self):
        """Jan 4 (perihelion) should have smallest distance."""
        result_jan1 = inv_rel_dist_earth_sun(1)
        result_jan4 = inv_rel_dist_earth_sun(4)
        assert result_jan4 < result_jan1

    def test_inv_rel_dist_earth_sun_july_aphelion(self):
        """July (aphelion) has largest distance, smallest inverse distance."""
        result_jan = inv_rel_dist_earth_sun(1)
        result_jul = inv_rel_dist_earth_sun(183)
        assert result_jul < result_jan

    def test_inv_rel_dist_earth_sun_positive(self):
        """Inverse relative distance should always be positive."""
        for doy in [1, 100, 200, 300, 365]:
            result = inv_rel_dist_earth_sun(doy)
            assert result > 0

    def test_daylight_hours_equator_equinox(self):
        """At equator on equinox, daylight should be ~12 hours."""
        sunrise_hour_angle = math.pi / 2  # 90° = π/2 rad
        result = daylight_hours(sunrise_hour_angle)
        assert 11.5 < result < 12.5

    def test_daylight_hours_increases_with_hour_angle(self):
        """Daylight hours increase with sunset hour angle."""
        result_small = daylight_hours(math.pi / 4)
        result_large = daylight_hours(math.pi / 2)
        assert result_large > result_small

    def test_net_in_sol_rad_default_albedo(self):
        """Net insolation with default albedo (0.23)."""
        incoming = 20  # MJ m-2 day-1
        result = net_in_sol_rad(incoming)
        expected = incoming * (1 - 0.23)
        assert result == pytest.approx(expected)

    def test_net_in_sol_rad_custom_albedo(self):
        """Net insolation with custom albedo."""
        incoming = 20
        result = net_in_sol_rad(incoming, albedo=0.25)
        expected = incoming * 0.75
        assert result == pytest.approx(expected)

    def test_net_in_sol_rad_zero_incoming(self):
        """Zero incoming radiation should give zero net radiation."""
        result = net_in_sol_rad(0)
        assert result == pytest.approx(0)


class TestEnergyConversions:
    """Test energy to evaporation conversion."""

    def test_energy2evap_positive(self):
        """Energy conversion should produce positive evaporation."""
        energy = 20  # MJ m-2 day-1
        result = energy2evap(energy)
        assert result > 0

    def test_energy2evap_proportional(self):
        """Doubling energy should roughly double evaporation."""
        result_10 = energy2evap(10)
        result_20 = energy2evap(20)
        assert result_20 > result_10

    def test_energy2evap_zero(self):
        """Zero energy should give zero evaporation."""
        result = energy2evap(0)
        assert result == pytest.approx(0)


class TestHargreavesMethod:
    """Test FAO hargreaves function."""

    def test_hargreaves_positive(self):
        """Hargreaves ETP should be positive."""
        tmin = 20
        tmax = 30
        tmean = 25
        et_rad = 20
        result = hargreaves(tmin, tmax, tmean, et_rad)
        assert result > 0

    def test_hargreaves_increases_with_radiation(self):
        """Higher radiation should increase ETP."""
        tmin = 20
        tmax = 30
        tmean = 25
        result_low = hargreaves(tmin, tmax, tmean, 10)
        result_high = hargreaves(tmin, tmax, tmean, 20)
        assert result_high > result_low

    def test_hargreaves_zero_radiation(self):
        """Zero radiation should give zero ETP."""
        result = hargreaves(20, 30, 25, 0)
        assert result == pytest.approx(0)


class TestFAO56PenmanMonteith:
    """Test FAO-56 Penman-Monteith reference ET."""

    def test_fao56_penman_monteith_positive(self):
        """FAO-56 PM should produce positive ETP."""
        net_rad = 15  # MJ m-2 day-1
        temp = 25  # °C
        wind_speed = 2  # m/s
        svp = 3.17  # kPa
        avp = 2.0  # kPa
        delta_svp_val = 0.19  # kPa °C-1
        psy = 0.067  # kPa °C-1

        result = fao56_penman_monteith(net_rad, temp, wind_speed, svp, avp, delta_svp_val, psy)
        assert result > 0

    def test_fao56_penman_monteith_increases_with_radiation(self):
        """Higher net radiation should increase ETP."""
        temp = 25
        wind_speed = 2
        svp = 3.17
        avp = 2.0
        delta_svp_val = 0.19
        psy = 0.067

        result_low = fao56_penman_monteith(10, temp, wind_speed, svp, avp, delta_svp_val, psy)
        result_high = fao56_penman_monteith(20, temp, wind_speed, svp, avp, delta_svp_val, psy)
        assert result_high > result_low

    def test_fao56_penman_monteith_wind_effect(self):
        """Higher wind speed should increase ETP."""
        net_rad = 15
        temp = 25
        svp = 3.17
        avp = 2.0
        delta_svp_val = 0.19
        psy = 0.067

        result_low_wind = fao56_penman_monteith(
            net_rad, temp, 1, svp, avp, delta_svp_val, psy
        )
        result_high_wind = fao56_penman_monteith(
            net_rad, temp, 5, svp, avp, delta_svp_val, psy
        )
        assert result_high_wind > result_low_wind

    def test_fao56_penman_monteith_humidity_effect(self):
        """Lower actual vapour pressure (drier) should increase ETP."""
        net_rad = 15
        temp = 25
        wind_speed = 2
        svp = 3.17
        delta_svp_val = 0.19
        psy = 0.067

        result_wet = fao56_penman_monteith(
            net_rad, temp, wind_speed, svp, 2.5, delta_svp_val, psy
        )
        result_dry = fao56_penman_monteith(
            net_rad, temp, wind_speed, svp, 1.5, delta_svp_val, psy
        )
        assert result_dry > result_wet

    def test_fao56_penman_monteith_with_soil_heat_flux(self):
        """Soil heat flux parameter should affect ETP calculation."""
        net_rad = 15
        temp = 25
        wind_speed = 2
        svp = 3.17
        avp = 2.0
        delta_svp_val = 0.19
        psy = 0.067

        result_no_flux = fao56_penman_monteith(
            net_rad, temp, wind_speed, svp, avp, delta_svp_val, psy, shf=0
        )
        result_with_flux = fao56_penman_monteith(
            net_rad, temp, wind_speed, svp, avp, delta_svp_val, psy, shf=2
        )
        assert result_no_flux != result_with_flux


class TestFAOFunctionsEdgeCases:
    """Test edge cases for FAO functions."""

    def test_high_altitude_pressure(self):
        """High altitude (5000 m) should have realistic pressure."""
        result = atm_pressure(5000)
        assert 50 < result < 60

    def test_very_low_temperature(self):
        """Very low temperature should still compute."""
        result = avp_from_tmin(-40)
        assert result > 0

    def test_very_high_temperature(self):
        """Very high temperature should still compute."""
        result = avp_from_tmin(50)
        assert result > 0

    def test_wide_temperature_range(self):
        """Wide daily range should be handled."""
        result = daily_mean_t(5, 45)
        assert result == pytest.approx(25)
