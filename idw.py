"""Inverse Distance Weighting (IDW) spatial interpolation."""

import numpy as np


def idw(
    x: np.ndarray,
    y: np.ndarray,
    v: np.ndarray,
    grid: np.ndarray,
    power: float,
) -> np.ndarray:
    """Fill a regular grid by Inverse Distance Weighting interpolation.

    Parameters
    ----------
    x, y : array_like
        Station row and column coordinates in grid-index space.
    v : array_like
        Observed values at each station. Same length as *x* and *y*.
    grid : np.ndarray
        Pre-allocated output array with shape ``(rows, cols)``. Filled in-place.
    power : float
        Distance-decay exponent. Common values: 1 (linear) to 3 (localized).

    Returns
    -------
    np.ndarray
        The filled *grid* array (same object passed in).

    Notes
    -----
    Cells where a station sits exactly (distance == 0) receive that station's
    value directly; all other cells use the standard weighted mean.

    The implementation is fully vectorised — no Python-level loops over grid
    cells — so it scales well to large grids.
    """
    rows, cols = grid.shape
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    v = np.asarray(v, dtype=float)

    # Grid index arrays, broadcast-ready: (rows, 1) and (1, cols)
    gi = np.arange(rows, dtype=float)[:, None]
    gj = np.arange(cols, dtype=float)[None, :]

    # Station coords reshaped for broadcasting: (n_stations, 1, 1)
    xs = x[:, None, None]
    ys = y[:, None, None]

    # Euclidean distances for every station × every cell: shape (n, rows, cols)
    dist = np.sqrt((xs - gi) ** 2 + (ys - gj) ** 2)

    # Weights: 1 / dist^power — becomes inf where a station is exactly on a cell.
    # inf/inf (exact-hit cells) raises an ignored invalid warning; those cells
    # are overwritten below by the exact-hit branch, so the NaN is harmless.
    with np.errstate(divide="ignore", invalid="ignore"):
        weights = dist ** (-power)                                       # (n, rows, cols)
        w_sum = weights.sum(axis=0)                                      # (rows, cols)
        result = (weights * v[:, None, None]).sum(axis=0) / w_sum        # (rows, cols)

    # Exact-hit cells: distance == 0 → assign that station's value directly
    exact = np.isinf(weights)          # shape (n, rows, cols)
    hit_cells = exact.any(axis=0)     # shape (rows, cols)
    if hit_cells.any():
        station_idx = exact.argmax(axis=0)
        result[hit_cells] = v[station_idx[hit_cells]]

    grid[:] = result
    return grid
