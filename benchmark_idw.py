#!/usr/bin/env python3
"""IDW vectorisation benchmark: loop-based vs vectorised implementation.

Runs both implementations on identical synthetic data and reports:
- Wall-clock time (best of N repeats)
- Speedup factor
- Maximum absolute difference between outputs (correctness check)
"""

import time

import numpy as np


# ── Reference implementation (original loop-based code) ──────────────────────


def _idw_loop_reference(
    station_rows: np.ndarray,
    station_cols: np.ndarray,
    station_values: np.ndarray,
    output_grid: np.ndarray,
    power: float,
) -> np.ndarray:
    """Original loop-based IDW — kept verbatim as the reference baseline."""
    for row in range(output_grid.shape[0]):
        for col in range(output_grid.shape[1]):
            distances = np.sqrt(
                (station_rows - row) ** 2 + (station_cols - col) ** 2
            )
            if (distances ** power).min() == 0:
                output_grid[row, col] = station_values[(distances ** power).argmin()]
            else:
                weight_total = np.sum(1 / (distances ** power))
                output_grid[row, col] = np.sum(
                    station_values / (distances ** power) / weight_total
                )
    return output_grid


# ── Timing helpers ────────────────────────────────────────────────────────────


def _time_idw_loop(
    station_rows: np.ndarray,
    station_cols: np.ndarray,
    station_values: np.ndarray,
    grid_shape: tuple,
    power: float,
    num_repeats: int,
) -> tuple[float, np.ndarray]:
    """Run the loop-based IDW *num_repeats* times and return best time + result.

    Args:
        station_rows: Station row coordinates.
        station_cols: Station column coordinates.
        station_values: Observed values at each station.
        grid_shape: ``(rows, cols)`` of the output grid.
        power: IDW distance-decay exponent.
        num_repeats: Number of timed runs; best (minimum) time is returned.

    Returns:
        Tuple of (best elapsed time in seconds, filled output grid).
    """
    elapsed_times = []
    output_grid = np.zeros(grid_shape)

    for _ in range(num_repeats):
        output_grid = np.zeros(grid_shape)
        start_time = time.perf_counter()
        _idw_loop_reference(station_rows, station_cols, station_values, output_grid, power)
        elapsed_times.append(time.perf_counter() - start_time)

    return min(elapsed_times), output_grid


def _time_idw_vectorised(
    station_rows: np.ndarray,
    station_cols: np.ndarray,
    station_values: np.ndarray,
    grid_shape: tuple,
    power: float,
    num_repeats: int,
) -> tuple[float, np.ndarray]:
    """Run the vectorised IDW *num_repeats* times and return best time + result.

    Args:
        station_rows: Station row coordinates.
        station_cols: Station column coordinates.
        station_values: Observed values at each station.
        grid_shape: ``(rows, cols)`` of the output grid.
        power: IDW distance-decay exponent.
        num_repeats: Number of timed runs; best (minimum) time is returned.

    Returns:
        Tuple of (best elapsed time in seconds, filled output grid).
    """
    from idw import idw

    elapsed_times = []
    output_grid = np.zeros(grid_shape)

    for _ in range(num_repeats):
        output_grid = np.zeros(grid_shape)
        start_time = time.perf_counter()
        idw(station_rows, station_cols, station_values, output_grid, power)
        elapsed_times.append(time.perf_counter() - start_time)

    return min(elapsed_times), output_grid


# ── Reporting ─────────────────────────────────────────────────────────────────


def _print_comparison(
    label: str,
    loop_time_seconds: float,
    vectorised_time_seconds: float,
    reference_grid: np.ndarray,
    comparison_grid: np.ndarray,
) -> None:
    """Print a formatted timing and correctness comparison.

    Args:
        label: Description of the grid size being benchmarked.
        loop_time_seconds: Best loop elapsed time in seconds.
        vectorised_time_seconds: Best vectorised elapsed time in seconds.
        reference_grid: Output from the loop-based implementation.
        comparison_grid: Output from the vectorised implementation.
    """
    max_absolute_difference = np.abs(reference_grid - comparison_grid).max()
    speedup = loop_time_seconds / vectorised_time_seconds
    outputs_match = max_absolute_difference < 1e-10

    print(f"\n{label}")
    print(f"  Loop-based :  {loop_time_seconds * 1000:.1f} ms")
    print(f"  Vectorised :  {vectorised_time_seconds * 1000:.1f} ms")
    print(f"  Speedup    :  {speedup:.1f}×")
    print(
        f"  Max |diff| :  {max_absolute_difference:.2e}  "
        f"{'✅ identical output' if outputs_match else '❌ outputs differ'}"
    )


# ── Main benchmark ────────────────────────────────────────────────────────────


def run(
    num_rows: int = 100,
    num_cols: int = 100,
    num_stations: int = 20,
    power: float = 2.0,
    num_repeats: int = 3,
) -> None:
    """Run the full benchmark suite and print results to stdout.

    Args:
        num_rows: Number of grid rows for the small benchmark.
        num_cols: Number of grid columns for the small benchmark.
        num_stations: Number of synthetic weather stations.
        power: IDW distance-decay exponent used in both implementations.
        num_repeats: Timed runs per configuration; best time is reported.
    """
    rng = np.random.default_rng(seed=42)

    station_rows = rng.uniform(0, num_rows - 1, num_stations)
    station_cols = rng.uniform(0, num_cols - 1, num_stations)
    station_values = rng.uniform(10, 100, num_stations)

    print(
        f"Config: grid={num_rows}×{num_cols}  "
        f"stations={num_stations}  power={power}  repeats={num_repeats}"
    )

    small_grid_shape = (num_rows, num_cols)
    loop_time, reference_grid = _time_idw_loop(
        station_rows, station_cols, station_values,
        small_grid_shape, power, num_repeats,
    )
    vectorised_time, vectorised_grid = _time_idw_vectorised(
        station_rows, station_cols, station_values,
        small_grid_shape, power, num_repeats,
    )
    _print_comparison(
        f"Small grid ({num_rows}×{num_cols})",
        loop_time, vectorised_time,
        reference_grid, vectorised_grid,
    )

    large_num_rows, large_num_cols = 300, 300
    large_station_rows = rng.uniform(0, large_num_rows - 1, num_stations)
    large_station_cols = rng.uniform(0, large_num_cols - 1, num_stations)

    large_loop_time, large_reference_grid = _time_idw_loop(
        large_station_rows, large_station_cols, station_values,
        (large_num_rows, large_num_cols), power, num_repeats,
    )
    large_vectorised_time, large_vectorised_grid = _time_idw_vectorised(
        large_station_rows, large_station_cols, station_values,
        (large_num_rows, large_num_cols), power, num_repeats,
    )
    _print_comparison(
        f"Large grid ({large_num_rows}×{large_num_cols})",
        large_loop_time, large_vectorised_time,
        large_reference_grid, large_vectorised_grid,
    )


if __name__ == "__main__":
    run()
