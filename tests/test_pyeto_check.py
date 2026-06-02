"""Unit tests for pyeto validation and checking functions.

Tests cover input validation for meteorological variables used
in FAO-56 evapotranspiration calculations.
"""

import math

import pytest

from pyeto._check import (
    check_day_hours,
    check_doy,
    check_latitude_rad,
    check_sol_dec_rad,
    check_sunset_hour_angle_rad,
)


class TestCheckDayOfYear:
    """Test day-of-year (DOY) validation."""

    def test_valid_doy_january_first(self):
        """January 1st (DOY 1) should be valid."""
        check_doy(1)

    def test_valid_doy_december_31st(self):
        """December 31st (DOY 365) should be valid."""
        check_doy(365)

    def test_valid_doy_leap_day(self):
        """Leap day (DOY 366) should be valid."""
        check_doy(366)

    def test_valid_doy_mid_year(self):
        """Mid-year DOY should be valid."""
        check_doy(180)

    def test_invalid_doy_zero(self):
        """DOY 0 should raise ValueError."""
        with pytest.raises(ValueError):
            check_doy(0)

    def test_invalid_doy_negative(self):
        """Negative DOY should raise ValueError."""
        with pytest.raises(ValueError):
            check_doy(-5)

    def test_invalid_doy_too_large(self):
        """DOY 367 should raise ValueError."""
        with pytest.raises(ValueError):
            check_doy(367)

    def test_invalid_doy_exceeds_year(self):
        """DOY greater than 366 should raise ValueError."""
        with pytest.raises(ValueError):
            check_doy(400)


class TestCheckLatitudeRad:
    """Test latitude (in radians) validation."""

    def test_valid_latitude_equator(self):
        """Latitude 0 (equator) should be valid."""
        check_latitude_rad(0)

    def test_valid_latitude_north_pole(self):
        """Northern hemisphere pole (π/2) should be valid."""
        check_latitude_rad(math.pi / 2)

    def test_valid_latitude_south_pole(self):
        """Southern hemisphere pole (-π/2) should be valid."""
        check_latitude_rad(-math.pi / 2)

    def test_valid_latitude_northern_temperate(self):
        """Northern temperate latitude should be valid."""
        check_latitude_rad(0.5)

    def test_valid_latitude_southern_temperate(self):
        """Southern temperate latitude should be valid."""
        check_latitude_rad(-0.5)

    def test_invalid_latitude_exceeds_north(self):
        """Latitude greater than π/2 should raise ValueError."""
        with pytest.raises(ValueError):
            check_latitude_rad(math.pi / 2 + 0.01)

    def test_invalid_latitude_exceeds_south(self):
        """Latitude less than -π/2 should raise ValueError."""
        with pytest.raises(ValueError):
            check_latitude_rad(-math.pi / 2 - 0.01)

    def test_invalid_latitude_pi(self):
        """Latitude π (exceeds bounds) should raise ValueError."""
        with pytest.raises(ValueError):
            check_latitude_rad(math.pi)


class TestCheckSolarDeclination:
    """Test solar declination (in radians) validation."""

    def test_valid_declination_zero(self):
        """Zero declination (equinoxes) should be valid."""
        check_sol_dec_rad(0)

    def test_valid_declination_positive_max(self):
        """Maximum positive declination should be valid."""
        check_sol_dec_rad(0.41)

    def test_valid_declination_negative_max(self):
        """Maximum negative declination should be valid."""
        check_sol_dec_rad(-0.41)

    def test_valid_declination_summer_solstice(self):
        """Summer solstice declination (~23.44°) should be valid."""
        summer_decl = math.radians(23.44)
        check_sol_dec_rad(summer_decl)

    def test_valid_declination_winter_solstice(self):
        """Winter solstice declination (~-23.44°) should be valid."""
        winter_decl = math.radians(-23.44)
        check_sol_dec_rad(winter_decl)

    def test_invalid_declination_exceeds_positive(self):
        """Declination greater than ~23.45° should raise ValueError."""
        with pytest.raises(ValueError):
            check_sol_dec_rad(0.42)

    def test_invalid_declination_exceeds_negative(self):
        """Declination less than ~-23.45° should raise ValueError."""
        with pytest.raises(ValueError):
            check_sol_dec_rad(-0.42)

    def test_invalid_declination_pi_over_2(self):
        """Declination π/2 (exceeds bounds) should raise ValueError."""
        with pytest.raises(ValueError):
            check_sol_dec_rad(math.pi / 2)


class TestCheckSunsetHourAngle:
    """Test sunset hour angle (in radians) validation."""

    def test_valid_sunset_hour_angle_zero(self):
        """Sunset hour angle 0 (pole) should be valid."""
        check_sunset_hour_angle_rad(0)

    def test_valid_sunset_hour_angle_pi(self):
        """Sunset hour angle π (equator, equinox) should be valid."""
        check_sunset_hour_angle_rad(math.pi)

    def test_valid_sunset_hour_angle_half_pi(self):
        """Sunset hour angle π/2 should be valid."""
        check_sunset_hour_angle_rad(math.pi / 2)

    def test_valid_sunset_hour_angle_quarter_pi(self):
        """Sunset hour angle π/4 should be valid."""
        check_sunset_hour_angle_rad(math.pi / 4)

    def test_valid_sunset_hour_angle_three_pi_fourths(self):
        """Sunset hour angle 3π/4 should be valid."""
        check_sunset_hour_angle_rad(3 * math.pi / 4)

    def test_invalid_sunset_hour_angle_negative(self):
        """Negative sunset hour angle should raise ValueError."""
        with pytest.raises(ValueError):
            check_sunset_hour_angle_rad(-0.1)

    def test_invalid_sunset_hour_angle_exceeds_pi(self):
        """Sunset hour angle greater than π should raise ValueError."""
        with pytest.raises(ValueError):
            check_sunset_hour_angle_rad(math.pi + 0.01)

    def test_invalid_sunset_hour_angle_two_pi(self):
        """Sunset hour angle 2π (exceeds bounds) should raise ValueError."""
        with pytest.raises(ValueError):
            check_sunset_hour_angle_rad(2 * math.pi)


class TestCheckDayHours:
    """Test daylight hours validation."""

    def test_valid_hours_zero(self):
        """Zero daylight hours (pole, winter) should be valid."""
        check_day_hours(0, "daylight_hours")

    def test_valid_hours_12(self):
        """12 daylight hours (equinox) should be valid."""
        check_day_hours(12, "daylight_hours")

    def test_valid_hours_24(self):
        """24 daylight hours (pole, summer) should be valid."""
        check_day_hours(24, "daylight_hours")

    def test_valid_hours_fractional(self):
        """Fractional daylight hours should be valid."""
        check_day_hours(14.5, "daylight_hours")

    def test_valid_hours_small_positive(self):
        """Small positive daylight hours should be valid."""
        check_day_hours(0.5, "daylight_hours")

    def test_invalid_hours_negative(self):
        """Negative daylight hours should raise ValueError."""
        with pytest.raises(ValueError):
            check_day_hours(-1, "daylight_hours")

    def test_invalid_hours_exceeds_24(self):
        """Daylight hours greater than 24 should raise ValueError."""
        with pytest.raises(ValueError):
            check_day_hours(25, "daylight_hours")

    def test_invalid_hours_way_too_large(self):
        """Daylight hours 100 (exceeds bounds) should raise ValueError."""
        with pytest.raises(ValueError):
            check_day_hours(100, "daylight_hours")

    def test_check_day_hours_with_custom_name(self):
        """Custom argument name should be used in error message."""
        with pytest.raises(ValueError, match="sunshine_hours"):
            check_day_hours(-5, "sunshine_hours")


class TestCheckValidationEdgeCases:
    """Test boundary conditions for all validation functions."""

    def test_doy_boundary_lower(self):
        """DOY at lower boundary (1) should be valid."""
        check_doy(1)

    def test_doy_boundary_upper(self):
        """DOY at upper boundary (366) should be valid."""
        check_doy(366)

    def test_latitude_boundary_north(self):
        """Latitude at North Pole should be valid."""
        check_latitude_rad(math.pi / 2)

    def test_latitude_boundary_south(self):
        """Latitude at South Pole should be valid."""
        check_latitude_rad(-math.pi / 2)

    def test_sunset_angle_boundary_zero(self):
        """Sunset hour angle at 0 should be valid."""
        check_sunset_hour_angle_rad(0)

    def test_sunset_angle_boundary_pi(self):
        """Sunset hour angle at π should be valid."""
        check_sunset_hour_angle_rad(math.pi)

    def test_day_hours_boundary_zero(self):
        """Day hours at 0 should be valid."""
        check_day_hours(0, "hours")

    def test_day_hours_boundary_24(self):
        """Day hours at 24 should be valid."""
        check_day_hours(24, "hours")
