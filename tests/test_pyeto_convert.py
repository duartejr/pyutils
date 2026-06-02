"""Unit tests for pyeto unit conversion functions.

Tests cover temperature, angle, and unit conversion utilities
used in FAO-56 evapotranspiration calculations.
"""

import math

import pytest

from pyeto.convert import celsius2kelvin, deg2rad, kelvin2celsius, rad2deg


class TestTemperatureConversions:
    """Test Celsius ↔ Kelvin temperature conversions."""

    def test_celsius_to_kelvin_zero(self):
        """Zero Celsius equals 273.15 Kelvin."""
        result = celsius2kelvin(0)
        assert result == pytest.approx(273.15)

    def test_celsius_to_kelvin_positive(self):
        """Positive Celsius values convert correctly."""
        result = celsius2kelvin(25)
        assert result == pytest.approx(298.15)

    def test_celsius_to_kelvin_negative(self):
        """Negative Celsius values convert correctly."""
        result = celsius2kelvin(-10)
        assert result == pytest.approx(263.15)

    def test_kelvin_to_celsius_absolute_zero(self):
        """Zero Kelvin converts to -273.15 Celsius."""
        result = kelvin2celsius(0)
        assert result == pytest.approx(-273.15)

    def test_kelvin_to_celsius_positive(self):
        """Positive Kelvin values convert correctly."""
        result = kelvin2celsius(298.15)
        assert result == pytest.approx(25)

    def test_kelvin_to_celsius_negative(self):
        """Negative Kelvin values convert correctly (though physically invalid)."""
        result = kelvin2celsius(263.15)
        assert result == pytest.approx(-10)

    def test_celsius_kelvin_roundtrip(self):
        """Converting Celsius → Kelvin → Celsius returns original value."""
        original = 37.5
        kelvin = celsius2kelvin(original)
        result = kelvin2celsius(kelvin)
        assert result == pytest.approx(original)

    def test_kelvin_celsius_roundtrip(self):
        """Converting Kelvin → Celsius → Kelvin returns original value."""
        original = 310.0
        celsius = kelvin2celsius(original)
        result = celsius2kelvin(celsius)
        assert result == pytest.approx(original)


class TestAngleConversions:
    """Test degree ↔ radian angle conversions."""

    def test_deg_to_rad_right_angle(self):
        """90 degrees equals π/2 radians."""
        result = deg2rad(90)
        assert result == pytest.approx(math.pi / 2)

    def test_deg_to_rad_full_circle(self):
        """360 degrees equals 2π radians."""
        result = deg2rad(360)
        assert result == pytest.approx(2 * math.pi)

    def test_deg_to_rad_straight_angle(self):
        """180 degrees equals π radians."""
        result = deg2rad(180)
        assert result == pytest.approx(math.pi)

    def test_deg_to_rad_zero(self):
        """Zero degrees equals zero radians."""
        result = deg2rad(0)
        assert result == pytest.approx(0)

    def test_deg_to_rad_negative(self):
        """Negative degrees convert correctly."""
        result = deg2rad(-45)
        assert result == pytest.approx(-math.pi / 4)

    def test_rad_to_deg_half_circle(self):
        """π radians equals 180 degrees."""
        result = rad2deg(math.pi)
        assert result == pytest.approx(180)

    def test_rad_to_deg_quarter_circle(self):
        """π/2 radians equals 90 degrees."""
        result = rad2deg(math.pi / 2)
        assert result == pytest.approx(90)

    def test_rad_to_deg_full_circle(self):
        """2π radians equals 360 degrees."""
        result = rad2deg(2 * math.pi)
        assert result == pytest.approx(360)

    def test_rad_to_deg_zero(self):
        """Zero radians equals zero degrees."""
        result = rad2deg(0)
        assert result == pytest.approx(0)

    def test_rad_to_deg_negative(self):
        """Negative radians convert correctly."""
        result = rad2deg(-math.pi / 4)
        assert result == pytest.approx(-45)

    def test_deg_rad_roundtrip(self):
        """Converting degrees → radians → degrees returns original value."""
        original = 37.5
        radians = deg2rad(original)
        result = rad2deg(radians)
        assert result == pytest.approx(original)

    def test_rad_deg_roundtrip(self):
        """Converting radians → degrees → radians returns original value."""
        original = 0.75
        degrees = rad2deg(original)
        result = deg2rad(degrees)
        assert result == pytest.approx(original)


class TestConversionEdgeCases:
    """Test edge cases and special values."""

    def test_large_temperature_values(self):
        """Large temperature values should convert correctly."""
        result = celsius2kelvin(1000)
        assert result == pytest.approx(1273.15)

    def test_very_large_temperature_kelvin(self):
        """Very large Kelvin values convert correctly."""
        result = kelvin2celsius(5000)
        assert result == pytest.approx(4726.85)

    def test_large_degree_values(self):
        """Large degree values should convert correctly."""
        result = deg2rad(720)
        assert result == pytest.approx(4 * math.pi)

    def test_large_radian_values(self):
        """Large radian values convert correctly."""
        result = rad2deg(4 * math.pi)
        assert result == pytest.approx(720)

    def test_fractional_temperature(self):
        """Fractional temperature values convert correctly."""
        result = celsius2kelvin(37.25)
        assert result == pytest.approx(310.4)

    def test_fractional_degrees(self):
        """Fractional degree values convert correctly."""
        result = deg2rad(45.5)
        assert result == pytest.approx(45.5 * math.pi / 180)

    def test_fractional_radians(self):
        """Fractional radian values convert correctly."""
        result = rad2deg(0.5)
        assert result == pytest.approx(0.5 * 180 / math.pi)
