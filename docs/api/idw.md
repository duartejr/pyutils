# idw — Inverse Distance Weighting

Vectorized IDW spatial interpolation. Fills a regular grid from station observations, where each cell receives a weighted mean of all stations with weights decaying as `distance ** (-power)`.

Cells coinciding exactly with a station receive that station's observed value directly.

!!! info "Performance"
    The vectorized implementation (NumPy broadcasting) is **30–45× faster** than a double-nested Python loop on typical grid sizes (100×100 to 300×300).

## Reference

Shepard, D. (1968). *A two-dimensional interpolation function for irregularly-spaced data.* ACM '68 Proceedings, 517-524.

---

## Main function

::: idw.idw
    options:
      show_signature_annotations: true

---

## Internal helpers

::: idw._build_distance_matrix

::: idw._compute_weighted_mean

::: idw._apply_exact_hits
