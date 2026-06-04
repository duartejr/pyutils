"""Unit tests for pyutils.climate module."""

import numpy as np
import pytest
from pyutils.climate import BiasCorrection
from pyutils.core import ValidationError


@pytest.fixture
def climate_data():
    np.random.seed(99)
    return {
        "hindcast":    np.random.normal(loc=15, scale=3, size=120),
        "observation": np.random.normal(loc=18, scale=4, size=120),
        "forecast":    np.random.normal(loc=14, scale=3, size=24),
    }


@pytest.fixture
def precipitation_data():
    np.random.seed(55)
    return {
        "hindcast":    np.abs(np.random.gamma(shape=2, scale=10, size=120)),
        "observation": np.abs(np.random.gamma(shape=2, scale=14, size=120)),
        "forecast":    np.abs(np.random.gamma(shape=2, scale=10, size=12)),
    }


# ------------------------------------------------------------------ #
# BiasCorrection — linear_scaling                                      #
# ------------------------------------------------------------------ #

class TestLinearScaling:
    def test_additive_shifts_mean(self, climate_data):
        corrected = BiasCorrection.linear_scaling(
            climate_data["forecast"],
            climate_data["hindcast"],
            climate_data["observation"],
            method="additive",
        )
        obs_mean  = np.nanmean(climate_data["observation"])
        hind_mean = np.nanmean(climate_data["hindcast"])
        fcst_mean = np.nanmean(climate_data["forecast"])
        expected_shift = obs_mean - hind_mean
        assert abs(np.mean(corrected) - (fcst_mean + expected_shift)) < 0.5

    def test_multiplicative_scales_mean(self, climate_data):
        corrected = BiasCorrection.linear_scaling(
            climate_data["forecast"],
            climate_data["hindcast"],
            climate_data["observation"],
            method="multiplicative",
        )
        ratio = np.nanmean(climate_data["observation"]) / np.nanmean(climate_data["hindcast"])
        expected = np.nanmean(climate_data["forecast"]) * ratio
        assert abs(np.mean(corrected) - expected) < 0.5

    def test_returns_same_length(self, climate_data):
        corrected = BiasCorrection.linear_scaling(
            climate_data["forecast"],
            climate_data["hindcast"],
            climate_data["observation"],
        )
        assert len(corrected) == len(climate_data["forecast"])

    def test_invalid_method_raises(self, climate_data):
        with pytest.raises(ValidationError):
            BiasCorrection.linear_scaling(
                climate_data["forecast"],
                climate_data["hindcast"],
                climate_data["observation"],
                method="unknown",
            )

    def test_mismatched_arrays_raise(self, climate_data):
        with pytest.raises(ValidationError):
            BiasCorrection.linear_scaling(
                climate_data["forecast"],
                hindcast=np.ones(60),       # different length from obs
                observation=np.ones(120),
            )


# ------------------------------------------------------------------ #
# BiasCorrection — variance_scaling                                    #
# ------------------------------------------------------------------ #

class TestVarianceScaling:
    def test_returns_same_length(self, climate_data):
        corrected = BiasCorrection.variance_scaling(
            climate_data["forecast"],
            climate_data["hindcast"],
            climate_data["observation"],
        )
        assert len(corrected) == len(climate_data["forecast"])

    def test_no_negative_values_for_precipitation(self, precipitation_data):
        corrected = BiasCorrection.variance_scaling(
            precipitation_data["forecast"],
            precipitation_data["hindcast"],
            precipitation_data["observation"],
        )
        assert np.all(corrected >= 0)

    def test_corrected_mean_closer_to_obs(self, climate_data):
        raw_diff = abs(
            np.mean(climate_data["forecast"]) - np.mean(climate_data["observation"])
        )
        corrected = BiasCorrection.variance_scaling(
            climate_data["forecast"],
            climate_data["hindcast"],
            climate_data["observation"],
        )
        corr_diff = abs(np.mean(corrected) - np.mean(climate_data["observation"]))
        assert corr_diff < raw_diff


# ------------------------------------------------------------------ #
# BiasCorrection — quantile_mapping_gamma                              #
# ------------------------------------------------------------------ #

class TestQuantileMappingGamma:
    def test_single_value(self, precipitation_data):
        corrected = BiasCorrection.quantile_mapping_gamma(
            np.array([5.0]),
            precipitation_data["hindcast"],
            precipitation_data["observation"],
        )
        assert corrected.shape == (1,)
        assert corrected[0] >= 0

    def test_array_input(self, precipitation_data):
        corrected = BiasCorrection.quantile_mapping_gamma(
            precipitation_data["forecast"],
            precipitation_data["hindcast"],
            precipitation_data["observation"],
        )
        assert len(corrected) == len(precipitation_data["forecast"])

    def test_zero_input_returns_zero(self, precipitation_data):
        corrected = BiasCorrection.quantile_mapping_gamma(
            np.array([0.0]),
            precipitation_data["hindcast"],
            precipitation_data["observation"],
        )
        assert corrected[0] == 0.0

    def test_no_negative_output(self, precipitation_data):
        corrected = BiasCorrection.quantile_mapping_gamma(
            precipitation_data["forecast"],
            precipitation_data["hindcast"],
            precipitation_data["observation"],
        )
        assert np.all(corrected >= 0)

    def test_insufficient_obs_raises(self):
        with pytest.raises(ValidationError):
            BiasCorrection.quantile_mapping_gamma(
                np.array([5.0]),
                hindcast=np.array([1.0]),   # too little data
                observation=np.array([1.0]),
            )
