"""Inverse Distance Weighting (IDW) spatial interpolation.

IDW estimates values at unsampled locations as a weighted average of nearby
observations, where the influence of each observation decays with distance.

References:
    Shepard, D. (1968). A two-dimensional interpolation function for
    irregularly-spaced data. ACM '68 Proceedings, 517-524.
"""

import numpy as np


def _build_distance_matrix(
    station_rows: np.ndarray,
    station_cols: np.ndarray,
    num_rows: int,
    num_cols: int,
) -> np.ndarray:
    """Compute Euclidean distances from every station to every grid cell.

    Args:
        station_rows: Row coordinates of each station. Shape ``(n_stations,)``.
        station_cols: Column coordinates of each station. Shape ``(n_stations,)``.
        num_rows: Number of rows in the output grid.
        num_cols: Number of columns in the output grid.

    Returns:
        Distance array of shape ``(n_stations, num_rows, num_cols)``.
    """
    row_indices = np.arange(num_rows, dtype=float)[:, None]   # (num_rows, 1)
    col_indices = np.arange(num_cols, dtype=float)[None, :]   # (1, num_cols)

    station_rows_3d = station_rows[:, None, None]  # (n_stations, 1, 1)
    station_cols_3d = station_cols[:, None, None]  # (n_stations, 1, 1)

    return np.sqrt(
        (station_rows_3d - row_indices) ** 2
        + (station_cols_3d - col_indices) ** 2
    )  # (n_stations, num_rows, num_cols)


def _compute_weighted_mean(
    distances: np.ndarray,
    station_values: np.ndarray,
    power: float,
) -> np.ndarray:
    """Return the IDW weighted mean for every grid cell.

    Cells that coincide with a station (distance == 0) are excluded from
    weight computation and handled separately by the caller.

    Args:
        distances: Shape ``(n_stations, rows, cols)``. Must not contain zeros
            (replace them before calling, e.g. with ``np.finfo(float).tiny``).
        station_values: Observed values. Shape ``(n_stations,)``.
        power: Distance-decay exponent.

    Returns:
        Interpolated array of shape ``(rows, cols)``.
    """
    # Overflow/invalid warnings are expected for exact-hit cells, whose
    # safe_distances == np.finfo(float).tiny raises weights to inf, then
    # inf/inf → NaN. Those cells are overwritten by _apply_exact_hits so
    # the NaN never appears in the final result.
    with np.errstate(over="ignore", invalid="ignore"):
        weights = distances ** (-power)  # (n_stations, rows, cols)
        weight_totals = weights.sum(axis=0)  # (rows, cols)
        return (weights * station_values[:, None, None]).sum(axis=0) / weight_totals


def _apply_exact_hits(
    interpolated_values: np.ndarray,
    exact_hit_mask: np.ndarray,
    station_values: np.ndarray,
) -> None:
    """Overwrite cells that coincide exactly with a station.

    Modifies *interpolated_values* in-place.

    Args:
        interpolated_values: Grid array to update. Shape ``(rows, cols)``.
        exact_hit_mask: Boolean array, True where a station sits on a cell.
            Shape ``(n_stations, rows, cols)``.
        station_values: Observed values. Shape ``(n_stations,)``.
    """
    cells_with_exact_hit = exact_hit_mask.any(axis=0)  # (rows, cols)
    if not cells_with_exact_hit.any():
        return

    # For each hit cell, pick the first station that lands there.
    station_index = exact_hit_mask.argmax(axis=0)  # (rows, cols)
    interpolated_values[cells_with_exact_hit] = station_values[
        station_index[cells_with_exact_hit]
    ]


def idw(
    station_rows: np.ndarray,
    station_cols: np.ndarray,
    station_values: np.ndarray,
    output_grid: np.ndarray,
    power: float,
) -> np.ndarray:
    """Fill a regular grid by Inverse Distance Weighting interpolation.

    Each grid cell receives a weighted mean of all station observations, where
    the weight of each station is inversely proportional to its distance from
    the cell raised to *power*. Grid cells that coincide exactly with a station
    receive that station's value directly, regardless of *power*.

    Memory note:
        This implementation allocates an intermediate array of shape
        ``(n_stations, rows, cols)`` to vectorise distance computation.
        For very large grids or a high station count this may be a constraint;
        consider chunking the grid in that case.

    Args:
        station_rows: Row coordinates of each station in grid-index space.
            Shape ``(n_stations,)``.
        station_cols: Column coordinates of each station in grid-index space.
            Shape ``(n_stations,)``.
        station_values: Observed value at each station. Shape ``(n_stations,)``.
        output_grid: Pre-allocated output array filled in-place.
            Shape ``(rows, cols)``.
        power: Distance-decay exponent. Higher values give more weight to
            nearby stations. Typical range: 1 (linear) to 3 (localised).

    Returns:
        The filled *output_grid* (same object passed in).

    Example:
        >>> import numpy as np
        >>> from idw import idw
        >>> grid = np.zeros((5, 5))
        >>> idw(
        ...     station_rows=np.array([0.0, 4.0]),
        ...     station_cols=np.array([0.0, 4.0]),
        ...     station_values=np.array([10.0, 20.0]),
        ...     output_grid=grid,
        ...     power=2,
        ... )
        >>> grid[0, 0]  # exact hit → station value
        10.0
    """
    num_rows, num_cols = output_grid.shape

    station_rows = np.asarray(station_rows, dtype=float)
    station_cols = np.asarray(station_cols, dtype=float)
    station_values = np.asarray(station_values, dtype=float)

    distances = _build_distance_matrix(station_rows, station_cols, num_rows, num_cols)

    # Detect exact hits via distance == 0, not via np.isinf(weights).
    # This correctly handles any value of power, including power == 0.
    exact_hit_mask = distances == 0  # (n_stations, rows, cols)

    # Replace zero distances with a tiny positive value so division is safe.
    # The affected cells are overwritten by _apply_exact_hits below.
    safe_distances = np.where(exact_hit_mask, np.finfo(float).tiny, distances)

    interpolated_values = _compute_weighted_mean(safe_distances, station_values, power)
    _apply_exact_hits(interpolated_values, exact_hit_mask, station_values)

    output_grid[:] = interpolated_values
    return output_grid
