"""Stream flow analysis and hydrological indices."""

from typing import Tuple, List, Optional
import numpy as np
import pandas as pd
from datetime import datetime


class FlowAnalyzer:
    """Analyze stream flow characteristics and compute hydrological indices.

    Consolidates vazoes_minimas.py functionality for flow duration curves,
    minimum flow quantiles, and extreme value analysis.
    """

    @staticmethod
    def flow_duration_curve(
        flow: np.ndarray,
        return_sorted: bool = True,
    ) -> Tuple[np.ndarray, np.ndarray]:
        """Compute flow duration curve probabilities.

        Parameters
        ----------
        flow : np.ndarray
            Stream flow series (m³/s or mm/month).
        return_sorted : bool, optional
            If True, return sorted flow values. Default is True.

        Returns
        -------
        tuple
            - exceedance_probability: Exceedance probability (%)
            - flow_values: Flow values sorted in descending order
        """
        flow = np.asarray(flow, dtype=float)
        sorted_flow = np.sort(flow)[::-1]  # Sort descending

        # Compute exceedance probability
        n = len(flow)
        exceedance_prob = np.linspace(100 / n, 100, n)

        return exceedance_prob, sorted_flow

    @staticmethod
    def minimum_flow_quantile(
        flow: np.ndarray,
        quantile: float = 90.0,
    ) -> float:
        """Compute minimum flow at specified exceedance probability.

        Parameters
        ----------
        flow : np.ndarray
            Stream flow series.
        quantile : float, optional
            Exceedance probability (0-100). Default is 90% (Q90).

        Returns
        -------
        float
            Minimum flow at specified quantile.
        """
        flow = np.asarray(flow, dtype=float)
        sorted_flow = np.sort(flow)[::-1]

        exceedance_prob = np.linspace(100 / len(flow), 100, len(flow))
        idx = np.argmin(np.abs(exceedance_prob - quantile))

        return sorted_flow[idx]

    @staticmethod
    def q7_10_extreme_value(
        flow: np.ndarray,
        dates: Optional[np.ndarray] = None,
        return_period: float = 10.0,
    ) -> float:
        """Compute 7-day minimum flow with specified return period.

        The 7-day, 10-year minimum flow (Q7,10) is a critical low-flow index
        used for water availability assessment and environmental flow standards.

        Parameters
        ----------
        flow : np.ndarray
            Daily or sub-daily flow series (m³/s).
        dates : np.ndarray, optional
            Corresponding dates (datetime or string). If provided, groups by year.
        return_period : float, optional
            Return period in years. Default is 10 years.

        Returns
        -------
        float
            7-day minimum flow with specified return period.
        """
        flow = np.asarray(flow, dtype=float)

        # Compute 7-day moving average
        if len(flow) < 7:
            return np.min(flow)

        window = 7
        moving_avg = np.convolve(
            flow, np.ones(window) / window, mode="valid"
        )

        # Group by year if dates provided
        if dates is not None:
            dates = pd.to_datetime(dates)
            df = pd.DataFrame({"flow": moving_avg, "date": dates[window - 1 :]})
            df["year"] = df["date"].dt.year

            # Get minimum for each year
            annual_mins = df.groupby("year")["flow"].min().values
        else:
            # If no dates, use all data
            annual_mins = np.array([np.min(moving_avg)])

        # Fit Gumbel extreme value distribution
        mean = np.mean(annual_mins)
        std = np.std(annual_mins)

        # Gumbel parameters
        c = np.sqrt(6) / np.pi
        sigma = std / c
        mu = mean - 0.5772 * sigma

        # Quantile for return period
        t = return_period
        quantile = mu - sigma * np.log(np.log(t / (t - 1)))

        return quantile

    @staticmethod
    def mean_annual_flow(
        flow: np.ndarray,
        dates: Optional[np.ndarray] = None,
    ) -> float:
        """Compute mean annual flow.

        Parameters
        ----------
        flow : np.ndarray
            Flow series.
        dates : np.ndarray, optional
            Dates for grouping. If None, uses entire series.

        Returns
        -------
        float
            Mean annual flow.
        """
        flow = np.asarray(flow, dtype=float)

        if dates is not None:
            dates = pd.to_datetime(dates)
            df = pd.DataFrame({"flow": flow, "date": dates})
            df["year"] = df["date"].dt.year
            annual_means = df.groupby("year")["flow"].mean()
            return annual_means.mean()
        else:
            return np.mean(flow)

    @staticmethod
    def low_flow_statistics(
        flow: np.ndarray,
        dates: Optional[np.ndarray] = None,
    ) -> dict:
        """Compute comprehensive low flow statistics.

        Parameters
        ----------
        flow : np.ndarray
            Flow series.
        dates : np.ndarray, optional
            Dates for annual grouping.

        Returns
        -------
        dict
            Dictionary containing:
            - q50: Median flow (m³/s)
            - q90: 90-percentile low flow
            - q95: 95-percentile low flow
            - q99: 99-percentile low flow
            - mean_annual: Mean annual flow
            - cv: Coefficient of variation
        """
        flow = np.asarray(flow, dtype=float)

        return {
            "q50": FlowAnalyzer.minimum_flow_quantile(flow, 50),
            "q90": FlowAnalyzer.minimum_flow_quantile(flow, 90),
            "q95": FlowAnalyzer.minimum_flow_quantile(flow, 95),
            "q99": FlowAnalyzer.minimum_flow_quantile(flow, 99),
            "mean_annual": FlowAnalyzer.mean_annual_flow(flow, dates),
            "cv": np.std(flow) / np.mean(flow),
        }

    @staticmethod
    def high_flow_statistics(
        flow: np.ndarray,
        dates: Optional[np.ndarray] = None,
    ) -> dict:
        """Compute comprehensive high flow statistics.

        Parameters
        ----------
        flow : np.ndarray
            Flow series.
        dates : np.ndarray, optional
            Dates for annual grouping.

        Returns
        -------
        dict
            Dictionary containing:
            - q1: 1-percentile high flow
            - q5: 5-percentile high flow
            - q10: 10-percentile high flow
            - max_annual: Maximum annual flow
            - mean_of_maxes: Mean of annual maximum flows
        """
        flow = np.asarray(flow, dtype=float)

        max_flows = []
        if dates is not None:
            dates = pd.to_datetime(dates)
            df = pd.DataFrame({"flow": flow, "date": dates})
            df["year"] = df["date"].dt.year
            max_flows = df.groupby("year")["flow"].max().values

        return {
            "q1": np.percentile(flow, 99),
            "q5": np.percentile(flow, 95),
            "q10": np.percentile(flow, 90),
            "max_flow": np.max(flow),
            "mean_of_maxes": np.mean(max_flows) if max_flows else np.max(flow),
        }
