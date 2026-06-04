"""Unit tests for pyutils.hydrology module."""

import numpy as np
import pytest
from pyutils.hydrology import (
    Hargreaves,
    Thornthwaite,
    PenmanMonteith,
    ThornthwaiteMather,
    StandardizedPrecipitationIndex,
    FlowAnalyzer,
)
from pyutils.core import ValidationError


# ------------------------------------------------------------------ #
# Fixtures                                                             #
# ------------------------------------------------------------------ #

@pytest.fixture
def monthly_temperatures():
    return {
        "t_med": np.array([22, 23, 24, 25, 26, 27, 27, 26, 25, 24, 23, 22], dtype=float),
        "t_max": np.array([27, 28, 29, 30, 31, 32, 32, 31, 30, 29, 28, 27], dtype=float),
        "t_min": np.array([17, 18, 19, 20, 21, 22, 22, 21, 20, 19, 18, 17], dtype=float),
    }


@pytest.fixture
def water_balance_inputs():
    return {
        "precipitation": np.array([200, 150, 100, 50, 30, 20, 15, 20, 40, 80, 150, 200], dtype=float),
        "potential_et":  np.array([ 40,  45,  60, 80,100,110,120,110, 90, 70,  50,  40], dtype=float),
    }


@pytest.fixture
def long_precipitation():
    np.random.seed(42)
    return np.random.gamma(shape=2, scale=40, size=120)   # 10 years monthly


@pytest.fixture
def daily_flow():
    np.random.seed(7)
    return np.abs(np.random.gamma(shape=2, scale=5, size=365 * 10))


# ------------------------------------------------------------------ #
# Hargreaves                                                           #
# ------------------------------------------------------------------ #

class TestHargreaves:
    def test_returns_twelve_values(self, monthly_temperatures):
        et = Hargreaves().compute(**monthly_temperatures, y=-5.0, months=1)
        assert et.shape == (12,)

    def test_all_values_positive(self, monthly_temperatures):
        et = Hargreaves().compute(**monthly_temperatures, y=-5.0, months=1)
        assert np.all(et > 0)

    def test_higher_temperature_range_gives_higher_et(self):
        base = dict(
            t_med=np.full(12, 25.0),
            t_max=np.full(12, 30.0),
            t_min=np.full(12, 20.0),
        )
        wide = dict(
            t_med=np.full(12, 25.0),
            t_max=np.full(12, 35.0),
            t_min=np.full(12, 15.0),
        )
        et_base = Hargreaves().compute(**base, y=-5.0, months=1)
        et_wide = Hargreaves().compute(**wide, y=-5.0, months=1)
        assert np.all(et_wide > et_base)

    def test_invalid_temperatures_raise(self):
        with pytest.raises(ValidationError):
            Hargreaves().compute(
                t_med=np.full(12, 35.0),
                t_max=np.full(12, 30.0),   # t_max < t_med — invalid
                t_min=np.full(12, 20.0),
                y=-5.0, months=1,
            )

    def test_validate_returns_true(self):
        assert Hargreaves().validate() is True

    def test_repr(self):
        assert "Hargreaves" in repr(Hargreaves())


# ------------------------------------------------------------------ #
# ThornthwaiteMather                                                   #
# ------------------------------------------------------------------ #

class TestThornthwaiteMather:
    def test_output_keys(self, water_balance_inputs):
        result = ThornthwaiteMather().compute(**water_balance_inputs)
        assert set(result) == {"runoff", "actual_et", "soil_moisture", "water_deficit", "water_surplus"}

    def test_output_length_matches_input(self, water_balance_inputs):
        result = ThornthwaiteMather().compute(**water_balance_inputs)
        n = len(water_balance_inputs["precipitation"])
        for key, arr in result.items():
            assert len(arr) == n, f"{key} length mismatch"

    def test_soil_moisture_never_exceeds_capacity(self, water_balance_inputs):
        cap = 100.0
        result = ThornthwaiteMather(cap=cap).compute(**water_balance_inputs)
        assert np.all(result["soil_moisture"] <= cap + 1e-9)

    def test_no_negative_values(self, water_balance_inputs):
        result = ThornthwaiteMather().compute(**water_balance_inputs)
        for key, arr in result.items():
            assert np.all(arr >= -1e-9), f"{key} contains negative values"

    def test_all_three_methods_run(self, water_balance_inputs):
        for method in ("standard", "variant2", "ccs"):
            result = ThornthwaiteMather(method=method).compute(**water_balance_inputs)
            assert "actual_et" in result

    def test_invalid_capacity_raises(self):
        with pytest.raises(ValidationError):
            ThornthwaiteMather(cap=-10)

    def test_invalid_method_raises(self):
        with pytest.raises(ValidationError):
            ThornthwaiteMather(method="unknown")

    def test_mismatched_arrays_raise(self):
        with pytest.raises(ValidationError):
            ThornthwaiteMather().compute(
                precipitation=np.ones(12),
                potential_et=np.ones(6),
            )

    def test_negative_precipitation_raises(self):
        with pytest.raises(ValidationError):
            ThornthwaiteMather().compute(
                precipitation=np.array([-1.0] * 12),
                potential_et=np.ones(12),
            )


# ------------------------------------------------------------------ #
# StandardizedPrecipitationIndex                                       #
# ------------------------------------------------------------------ #

class TestSPI:
    def test_output_length(self, long_precipitation):
        spi, _, _ = StandardizedPrecipitationIndex(timescale=3).compute(long_precipitation)
        assert len(spi) == len(long_precipitation)

    def test_spi_roughly_standard_normal(self, long_precipitation):
        spi, _, _ = StandardizedPrecipitationIndex(timescale=3).compute(long_precipitation)
        valid = spi[~np.isnan(spi)]
        assert abs(np.mean(valid)) < 0.5
        assert 0.5 < np.std(valid) < 2.0

    def test_timescale_1_has_fewer_nans(self, long_precipitation):
        spi1, _, _ = StandardizedPrecipitationIndex(timescale=1).compute(long_precipitation)
        spi3, _, _ = StandardizedPrecipitationIndex(timescale=3).compute(long_precipitation)
        assert np.sum(np.isnan(spi1)) < np.sum(np.isnan(spi3))

    def test_lognormal_distribution(self, long_precipitation):
        spi, _, _ = StandardizedPrecipitationIndex(timescale=3, distribution="lognormal").compute(long_precipitation)
        valid = spi[~np.isnan(spi)]
        assert len(valid) > 0

    def test_short_series_raises(self):
        with pytest.raises(ValidationError):
            StandardizedPrecipitationIndex(timescale=12).compute(np.ones(5))

    def test_negative_precipitation_raises(self):
        with pytest.raises(ValidationError):
            StandardizedPrecipitationIndex().compute(np.array([-1.0] * 50))

    def test_invalid_timescale_raises(self):
        with pytest.raises(ValidationError):
            StandardizedPrecipitationIndex(timescale=0)

    def test_invalid_distribution_raises(self):
        with pytest.raises(ValidationError):
            StandardizedPrecipitationIndex(distribution="weibull")


# ------------------------------------------------------------------ #
# FlowAnalyzer                                                         #
# ------------------------------------------------------------------ #

class TestFlowAnalyzer:
    def test_duration_curve_length(self, daily_flow):
        prob, flows = FlowAnalyzer.flow_duration_curve(daily_flow)
        assert len(prob) == len(daily_flow)
        assert len(flows) == len(daily_flow)

    def test_duration_curve_descending(self, daily_flow):
        _, flows = FlowAnalyzer.flow_duration_curve(daily_flow)
        assert np.all(flows[:-1] >= flows[1:])

    def test_exceedance_probability_range(self, daily_flow):
        prob, _ = FlowAnalyzer.flow_duration_curve(daily_flow)
        assert prob[0] > 0
        assert prob[-1] <= 100.0

    def test_q90_less_than_q50(self, daily_flow):
        q50 = FlowAnalyzer.minimum_flow_quantile(daily_flow, 50)
        q90 = FlowAnalyzer.minimum_flow_quantile(daily_flow, 90)
        assert q90 < q50

    def test_q99_less_than_q95(self, daily_flow):
        q95 = FlowAnalyzer.minimum_flow_quantile(daily_flow, 95)
        q99 = FlowAnalyzer.minimum_flow_quantile(daily_flow, 99)
        assert q99 < q95

    def test_low_flow_statistics_keys(self, daily_flow):
        stats = FlowAnalyzer.low_flow_statistics(daily_flow)
        assert set(stats) == {"q50", "q90", "q95", "q99", "mean_annual", "cv"}

    def test_high_flow_statistics_keys(self, daily_flow):
        stats = FlowAnalyzer.high_flow_statistics(daily_flow)
        assert "q1" in stats and "max_flow" in stats

    def test_cv_positive(self, daily_flow):
        stats = FlowAnalyzer.low_flow_statistics(daily_flow)
        assert stats["cv"] > 0

    def test_q7_10_positive(self, daily_flow):
        q = FlowAnalyzer.q7_10_extreme_value(daily_flow)
        assert q > 0

    def test_mean_annual_flow(self, daily_flow):
        mean = FlowAnalyzer.mean_annual_flow(daily_flow)
        assert mean > 0
        assert abs(mean - np.mean(daily_flow)) < np.std(daily_flow)


# ------------------------------------------------------------------ #
# PenmanMonteith                                                       #
# ------------------------------------------------------------------ #

class TestPenmanMonteith:
    def test_returns_positive_et0(self):
        pm = PenmanMonteith()
        et0 = pm.compute(
            t_mean=25.0, t_max=32.0, t_min=18.0,
            rh_mean=60.0, u2=2.0,
            latitude=-15.0, day_of_year=180,
            elevation=500.0,
        )
        assert et0 > 0

    def test_higher_temperature_gives_higher_et0(self):
        pm = PenmanMonteith()
        et_cool = pm.compute(
            t_mean=15.0, t_max=20.0, t_min=10.0,
            rh_mean=60.0, u2=2.0,
            latitude=0.0, day_of_year=180,
        )
        et_warm = pm.compute(
            t_mean=30.0, t_max=35.0, t_min=25.0,
            rh_mean=60.0, u2=2.0,
            latitude=0.0, day_of_year=180,
        )
        assert et_warm > et_cool

    def test_with_provided_rs(self):
        pm = PenmanMonteith()
        et0 = pm.compute(
            t_mean=25.0, t_max=32.0, t_min=18.0,
            rh_mean=60.0, u2=2.0,
            latitude=-15.0, day_of_year=180,
            rs=18.0,
        )
        assert et0 > 0

    def test_invalid_temperatures_raises(self):
        with pytest.raises(ValidationError):
            PenmanMonteith().compute(
                t_mean=35.0, t_max=30.0, t_min=18.0,  # t_mean > t_max
                rh_mean=60.0, u2=2.0,
                latitude=-15.0, day_of_year=180,
            )

    def test_invalid_day_of_year_raises(self):
        with pytest.raises(ValidationError):
            PenmanMonteith().compute(
                t_mean=25.0, t_max=32.0, t_min=18.0,
                rh_mean=60.0, u2=2.0,
                latitude=-15.0, day_of_year=400,
            )

    def test_validate_returns_true(self):
        assert PenmanMonteith().validate() is True
