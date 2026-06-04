"""Unit tests for pyutils.utils module."""

import math
import numpy as np
import pytest
from pyutils.utils import DataValidator, UnitConverter
from pyutils.core import ValidationError


# ------------------------------------------------------------------ #
# DataValidator                                                        #
# ------------------------------------------------------------------ #

class TestDataValidatorPrecipitation:
    def test_valid_series_passes(self):
        assert DataValidator.validate_precipitation(np.array([0.0, 10.0, 50.0])) is True

    def test_negative_values_raise(self):
        with pytest.raises(ValidationError):
            DataValidator.validate_precipitation(np.array([-1.0, 10.0]))

    def test_nan_raises_by_default(self):
        with pytest.raises(ValidationError):
            DataValidator.validate_precipitation(np.array([np.nan, 10.0]))

    def test_nan_allowed_when_flag_set(self):
        assert DataValidator.validate_precipitation(
            np.array([np.nan, 10.0]), allow_missing=True
        ) is True

    def test_zero_disallowed(self):
        with pytest.raises(ValidationError):
            DataValidator.validate_precipitation(
                np.array([0.0, 10.0]), allow_zero=False
            )

    def test_excessive_value_raises(self):
        with pytest.raises(ValidationError):
            DataValidator.validate_precipitation(np.array([5000.0]))


class TestDataValidatorTemperature:
    def test_valid_range_passes(self):
        assert DataValidator.validate_temperature(np.array([-20.0, 0.0, 35.0])) is True

    def test_too_cold_raises(self):
        with pytest.raises(ValidationError):
            DataValidator.validate_temperature(np.array([-200.0]))

    def test_too_hot_raises(self):
        with pytest.raises(ValidationError):
            DataValidator.validate_temperature(np.array([70.0]))

    def test_nan_raises_by_default(self):
        with pytest.raises(ValidationError):
            DataValidator.validate_temperature(np.array([np.nan]))

    def test_nan_allowed_when_flag_set(self):
        assert DataValidator.validate_temperature(
            np.array([np.nan, 20.0]), allow_missing=True
        ) is True


class TestDataValidatorOther:
    def test_humidity_valid(self):
        assert DataValidator.validate_relative_humidity(np.array([0.0, 50.0, 100.0])) is True

    def test_humidity_negative_raises(self):
        with pytest.raises(ValidationError):
            DataValidator.validate_relative_humidity(np.array([-5.0]))

    def test_humidity_over_100_raises(self):
        with pytest.raises(ValidationError):
            DataValidator.validate_relative_humidity(np.array([101.0]))

    def test_wind_speed_valid(self):
        assert DataValidator.validate_wind_speed(np.array([0.0, 5.0, 20.0])) is True

    def test_negative_wind_raises(self):
        with pytest.raises(ValidationError):
            DataValidator.validate_wind_speed(np.array([-1.0]))

    def test_spatial_data_valid(self):
        assert DataValidator.validate_spatial_data(
            np.array([0.0, 1.0]), np.array([0.0, 1.0]), np.array([10.0, 20.0])
        ) is True

    def test_spatial_mismatch_raises(self):
        with pytest.raises(ValidationError):
            DataValidator.validate_spatial_data(
                np.array([0.0, 1.0]),
                np.array([0.0]),          # length mismatch
                np.array([10.0, 20.0]),
            )

    def test_spatial_single_point_raises(self):
        with pytest.raises(ValidationError):
            DataValidator.validate_spatial_data(
                np.array([0.0]),
                np.array([0.0]),
                np.array([10.0]),
            )

    def test_time_series_length_valid(self):
        assert DataValidator.validate_time_series_length(np.ones(30), min_length=12) is True

    def test_time_series_too_short_raises(self):
        with pytest.raises(ValidationError):
            DataValidator.validate_time_series_length(np.ones(5), min_length=12)


# ------------------------------------------------------------------ #
# UnitConverter — temperature                                          #
# ------------------------------------------------------------------ #

class TestUnitConverterTemperature:
    def test_celsius_to_kelvin_freezing(self):
        assert UnitConverter.celsius_to_kelvin(0.0) == pytest.approx(273.15)

    def test_celsius_to_kelvin_boiling(self):
        assert UnitConverter.celsius_to_kelvin(100.0) == pytest.approx(373.15)

    def test_kelvin_to_celsius_roundtrip(self):
        original = 25.0
        assert UnitConverter.kelvin_to_celsius(
            UnitConverter.celsius_to_kelvin(original)
        ) == pytest.approx(original)

    def test_fahrenheit_to_celsius_freezing(self):
        assert UnitConverter.fahrenheit_to_celsius(32.0) == pytest.approx(0.0)

    def test_fahrenheit_to_celsius_boiling(self):
        assert UnitConverter.fahrenheit_to_celsius(212.0) == pytest.approx(100.0)

    def test_celsius_to_fahrenheit_roundtrip(self):
        original = 37.0
        assert UnitConverter.fahrenheit_to_celsius(
            UnitConverter.celsius_to_fahrenheit(original)
        ) == pytest.approx(original)


# ------------------------------------------------------------------ #
# UnitConverter — angles                                               #
# ------------------------------------------------------------------ #

class TestUnitConverterAngles:
    def test_180_degrees_is_pi(self):
        assert UnitConverter.degrees_to_radians(180.0) == pytest.approx(math.pi)

    def test_pi_is_180_degrees(self):
        assert UnitConverter.radians_to_degrees(math.pi) == pytest.approx(180.0)

    def test_roundtrip(self):
        original = 45.0
        assert UnitConverter.radians_to_degrees(
            UnitConverter.degrees_to_radians(original)
        ) == pytest.approx(original)


# ------------------------------------------------------------------ #
# UnitConverter — depth / radiation                                    #
# ------------------------------------------------------------------ #

class TestUnitConverterOther:
    def test_mm_to_meters(self):
        assert UnitConverter.mm_to_meters(1000.0) == pytest.approx(1.0)

    def test_meters_to_mm(self):
        assert UnitConverter.meters_to_mm(1.0) == pytest.approx(1000.0)

    def test_depth_roundtrip(self):
        original = 250.0
        assert UnitConverter.meters_to_mm(
            UnitConverter.mm_to_meters(original)
        ) == pytest.approx(original)

    def test_radiation_roundtrip(self):
        original = 15.0
        assert UnitConverter.watts_per_m2_to_mj_per_m2(
            UnitConverter.mj_per_m2_to_watts_per_m2(original)
        ) == pytest.approx(original, rel=1e-4)

    def test_wind_speed_correction_reduces_at_height(self):
        corrected = UnitConverter.wind_speed_height_correction(5.0, 10.0, 2.0)
        assert corrected < 5.0
