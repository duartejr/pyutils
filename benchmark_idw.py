#!/usr/bin/env python3
"""IDW vectorisation benchmark — before vs after.

Runs the old loop-based implementation alongside the new vectorised one on
the same synthetic data and prints timing and correctness results.
"""

import time
import numpy as np


# ── old implementation (loop-based, preserved here for comparison) ────────────

def _idw_loop(x, y, v, grid, power):
    for i in range(grid.shape[0]):
        for j in range(grid.shape[1]):
            distance = np.sqrt((x - i) ** 2 + (y - j) ** 2)
            if (distance ** power).min() == 0:
                grid[i, j] = v[(distance ** power).argmin()]
            else:
                total = np.sum(1 / (distance ** power))
                grid[i, j] = np.sum(v / (distance ** power) / total)
    return grid


# ── new implementation (vectorised) ──────────────────────────────────────────

def _idw_vec(x, y, v, grid, power):
    rows, cols = grid.shape
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    v = np.asarray(v, dtype=float)

    gi = np.arange(rows, dtype=float)[:, None]
    gj = np.arange(cols, dtype=float)[None, :]
    xs = x[:, None, None]
    ys = y[:, None, None]

    dist = np.sqrt((xs - gi) ** 2 + (ys - gj) ** 2)

    with np.errstate(divide="ignore", invalid="ignore"):
        weights = dist ** (-power)

    w_sum = weights.sum(axis=0)
    result = (weights * v[:, None, None]).sum(axis=0) / w_sum

    exact = np.isinf(weights)
    hit_cells = exact.any(axis=0)
    if hit_cells.any():
        station_idx = exact.argmax(axis=0)
        result[hit_cells] = v[station_idx[hit_cells]]

    grid[:] = result
    return grid


# ── benchmark ─────────────────────────────────────────────────────────────────

def run(rows=100, cols=100, n_stations=20, power=2, repeats=3):
    rng = np.random.default_rng(42)
    x = rng.uniform(0, rows - 1, n_stations)
    y = rng.uniform(0, cols - 1, n_stations)
    v = rng.uniform(10, 100, n_stations)

    print(f"Grid: {rows}×{cols}  |  Stations: {n_stations}  |  Power: {power}  |  Repeats: {repeats}\n")

    # --- loop ---
    times_loop = []
    for _ in range(repeats):
        g = np.zeros((rows, cols))
        t0 = time.perf_counter()
        _idw_loop(x, y, v, g, power)
        times_loop.append(time.perf_counter() - t0)
    grid_loop = g.copy()
    t_loop = min(times_loop)

    # --- vectorised ---
    times_vec = []
    for _ in range(repeats):
        g = np.zeros((rows, cols))
        t0 = time.perf_counter()
        _idw_vec(x, y, v, g, power)
        times_vec.append(time.perf_counter() - t0)
    grid_vec = g.copy()
    t_vec = min(times_vec)

    # --- results ---
    max_diff = np.abs(grid_loop - grid_vec).max()
    print(f"Loop-based:   {t_loop*1000:.1f} ms")
    print(f"Vectorised:   {t_vec*1000:.1f} ms")
    print(f"Speedup:      {t_loop/t_vec:.1f}×")
    print(f"Max abs diff: {max_diff:.2e}  {'✅ identical output' if max_diff < 1e-10 else '❌ outputs differ'}\n")

    # larger grid
    rows2, cols2 = 300, 300
    print(f"Larger grid: {rows2}×{cols2}")
    g_vec = np.zeros((rows2, cols2))
    x2 = rng.uniform(0, rows2 - 1, n_stations)
    y2 = rng.uniform(0, cols2 - 1, n_stations)
    t0 = time.perf_counter()
    _idw_vec(x2, y2, v, g_vec, power)
    print(f"Vectorised:   {(time.perf_counter()-t0)*1000:.1f} ms")

    g_loop = np.zeros((rows2, cols2))
    t0 = time.perf_counter()
    _idw_loop(x2, y2, v, g_loop, power)
    print(f"Loop-based:   {(time.perf_counter()-t0)*1000:.1f} ms")
    diff2 = np.abs(g_loop - g_vec).max()
    print(f"Max abs diff: {diff2:.2e}  {'✅ identical output' if diff2 < 1e-10 else '❌ outputs differ'}")


if __name__ == "__main__":
    run()
